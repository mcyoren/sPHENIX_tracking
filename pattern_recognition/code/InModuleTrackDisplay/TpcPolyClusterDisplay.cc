#include "TpcPolyClusterDisplay.h"

#include "FinalTrack.h"
#include "FinalTrackContainer.h"
#include "FinalTrackVertex.h"
#include "FinalTrackVertexContainer.h"
#include "TpcPolyClusterTrack.h"
#include "TpcPolyClusterTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH3D.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <vector>

namespace
{
  int cluster_color(const unsigned int icluster)
  {
    static const int colors[] = {
      kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kCyan + 2,
      kOrange + 7, kViolet + 1, kAzure + 1, kPink + 7, kTeal + 3
    };
    return colors[icluster % (sizeof(colors) / sizeof(colors[0]))];
  }



  bool cluster_track_z_range(const TpcPolyClusterTrack* trk,
                             const double display_zmin,
                             const double display_zmax,
                             double& zmin,
                             double& zmax)
  {
    if (!trk || !trk->isValid()) return false;

    zmin = std::numeric_limits<double>::max();
    zmax = -std::numeric_limits<double>::max();

    auto update_range = [&](const double z)
    {
      if (!std::isfinite(z)) return;
      zmin = std::min(zmin, z);
      zmax = std::max(zmax, z);
    };

    for (unsigned int icluster = 0; icluster < trk->size_clusters(); ++icluster)
    {
      update_range(trk->get_cluster_z(icluster));
    }

    if (zmin == std::numeric_limits<double>::max() || zmax == -std::numeric_limits<double>::max()) return false;

    zmin = std::max(zmin, display_zmin);
    zmax = std::min(zmax, display_zmax);
    if (zmin > zmax) return false;

    if (zmin == zmax)
    {
      zmin = std::max(zmin - 0.1, display_zmin);
      zmax = std::min(zmax + 0.1, display_zmax);
    }

    return zmin < zmax;
  }

  bool track_vertex_z_selected(const FinalTrackVertex* vtx,
                               const double vertex_zmin,
                               const double vertex_zmax)
  {
    if (!vtx) return false;
    const double z0 = vtx->get_z0();
    return std::isfinite(z0) && z0 >= vertex_zmin && z0 <= vertex_zmax;
  }

  bool final_track_xy_at_z(const FinalTrack* trk,
                           const double z,
                           const double magnetic_field_tesla,
                           const double arc_direction,
                           double& x,
                           double& y)
  {
    if (!trk || trk->get_fit_status() == 0 || !std::isfinite(z)) return false;

    const double x0 = trk->get_x();
    const double y0 = trk->get_y();
    const double z0 = trk->get_z();
    const double px = trk->get_px();
    const double py = trk->get_py();
    const double pz = trk->get_pz();
    const double charge = trk->get_charge();
    if (!std::isfinite(x0) || !std::isfinite(y0) || !std::isfinite(z0) ||
        !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz) ||
        !std::isfinite(charge))
    {
      return false;
    }

    const double pt = std::hypot(px, py);
    if (pt <= 0.0 || std::fabs(pz) < 1.0e-12) return false;

    const double dz = z - z0;
    if (std::fabs(charge * magnetic_field_tesla) < 1.0e-12)
    {
      x = x0 + arc_direction * px / pz * dz;
      y = y0 + arc_direction * py / pz * dz;
      return std::isfinite(x) && std::isfinite(y);
    }

    const double signed_radius = pt / (0.003 * charge * magnetic_field_tesla);
    const double radius = std::fabs(signed_radius);
    if (!std::isfinite(radius) || radius <= 0.0) return false;

    const double tx = px / pt;
    const double ty = py / pt;
    const double sign = signed_radius > 0.0 ? 1.0 : -1.0;
    const double xc = x0 + sign * radius * ty;
    const double yc = y0 - sign * radius * tx;
    const double phi0 = std::atan2(y0 - yc, x0 - xc);
    const double dzds = pz / pt;
    if (std::fabs(dzds) < 1.0e-12) return false;

    const double arc = arc_direction * dz / dzds;
    const double phi = phi0 - sign * arc / radius;
    x = xc + radius * std::cos(phi);
    y = yc + radius * std::sin(phi);
    return std::isfinite(x) && std::isfinite(y);
  }

  double cluster_line_residual2(const FinalTrack* final_track,
                                const TpcPolyClusterTrack* cluster_track,
                                const double magnetic_field_tesla,
                                const double arc_direction)
  {
    if (!final_track || !cluster_track || !cluster_track->isValid()) return std::numeric_limits<double>::max();

    double sum = 0.0;
    unsigned int n = 0;
    for (unsigned int icluster = 0; icluster < cluster_track->size_clusters(); ++icluster)
    {
      const double cx = cluster_track->get_cluster_x(icluster);
      const double cy = cluster_track->get_cluster_y(icluster);
      const double cz = cluster_track->get_cluster_z(icluster);
      if (!std::isfinite(cx) || !std::isfinite(cy) || !std::isfinite(cz)) continue;

      double x = 0.0;
      double y = 0.0;
      if (!final_track_xy_at_z(final_track, cz, magnetic_field_tesla, arc_direction, x, y)) continue;

      const double dx = x - cx;
      const double dy = y - cy;
      sum += dx * dx + dy * dy;
      ++n;
    }

    return n > 0 ? sum / static_cast<double>(n) : std::numeric_limits<double>::max();
  }

  TPolyLine3D* make_final_track_line(const FinalTrack* trk,
                                     const double zmin,
                                     const double zmax,
                                     const double xymax,
                                     const double magnetic_field_tesla,
                                     const double arc_direction,
                                     const int color)
  {
    if (!trk || trk->get_fit_status() == 0) return nullptr;

    std::vector<double> zs;
    std::vector<double> xs;
    std::vector<double> ys;
    const unsigned int nsteps = 80;
    zs.reserve(nsteps + 1);
    xs.reserve(nsteps + 1);
    ys.reserve(nsteps + 1);

    for (unsigned int istep = 0; istep <= nsteps; ++istep)
    {
      const double f = static_cast<double>(istep) / static_cast<double>(nsteps);
      const double z = zmin + f * (zmax - zmin);
      double x = 0.0;
      double y = 0.0;
      if (!final_track_xy_at_z(trk, z, magnetic_field_tesla, arc_direction, x, y)) continue;
      if (std::fabs(x) > xymax || std::fabs(y) > xymax) continue;
      zs.push_back(z);
      xs.push_back(x);
      ys.push_back(y);
    }

    if (zs.size() < 2) return nullptr;
    TPolyLine3D* line = new TPolyLine3D(static_cast<int>(zs.size()));
    for (unsigned int i = 0; i < zs.size(); ++i) line->SetPoint(static_cast<int>(i), zs[i], xs[i], ys[i]);
    line->SetLineColor(color);
    line->SetLineWidth(3);
    return line;
  }

  TPolyMarker3D* make_cluster_marker(const double z,
                                     const double x,
                                     const double y,
                                     const int color)
  {
    TPolyMarker3D* marker = new TPolyMarker3D(1);
    marker->SetPoint(0, z, x, y);
    marker->SetMarkerColor(color);
    marker->SetMarkerStyle(20);
    marker->SetMarkerSize(1.4);
    return marker;
  }

  TPolyMarker3D* make_pca_marker(const FinalTrackVertex* vtx,
                                 const double zmin,
                                 const double zmax,
                                 const double xymax,
                                 const int color)
  {
    if (!vtx || !vtx->get_pca_valid()) return nullptr;

    const double x = vtx->get_pca_x();
    const double y = vtx->get_pca_y();
    const double z = vtx->get_pca_z();
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) return nullptr;
    if (z < zmin || z > zmax) return nullptr;
    if (std::fabs(x) > xymax || std::fabs(y) > xymax) return nullptr;

    TPolyMarker3D* marker = new TPolyMarker3D(1);
    marker->SetPoint(0, z, x, y);
    marker->SetMarkerColor(color);
    marker->SetMarkerStyle(20);
    marker->SetMarkerSize(1.2);
    return marker;
  }

  TPolyMarker3D* make_collision_vertex_marker(const double x,
                                              const double y,
                                              const double z,
                                              const double zmin,
                                              const double zmax,
                                              const double xymax)
  {
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) return nullptr;
    if (z < zmin || z > zmax) return nullptr;
    if (std::fabs(x) > xymax || std::fabs(y) > xymax) return nullptr;

    TPolyMarker3D* marker = new TPolyMarker3D(1);
    marker->SetPoint(0, z, x, y);
    marker->SetMarkerColor(kBlack);
    marker->SetMarkerStyle(29);
    marker->SetMarkerSize(2.0);
    return marker;
  }
}

TpcPolyClusterDisplay::TpcPolyClusterDisplay(const std::string& name,
                                             const std::string& outfilename,
                                             const std::string& clusterNodeName,
                                             const unsigned int maxEventDisplays)
  : SubsysReco(name)
  , m_outfilename(outfilename)
  , m_clusterNodeName(clusterNodeName)
  , m_finalTrackNodeName("FINALTRACKS")
  , m_finalTrackVertexNodeName("FINALTRACKVERTICES")
  , m_maxEventDisplays(maxEventDisplays)
  , m_evt(0)
  , m_eventsSaved(0)
  , m_zmin(-102.0)
  , m_zmax(102.0)
  , m_trackVertexZMin(-20.0)
  , m_trackVertexZMax(20.0)
  , m_xymax(85.0)
  , m_magneticFieldTesla(1.4)
  , m_outfile(nullptr)
  , m_clusterTracks(nullptr)
  , m_finalTracks(nullptr)
  , m_finalTrackVertices(nullptr)
{
}

TpcPolyClusterDisplay::~TpcPolyClusterDisplay()
{
  if (m_outfile)
  {
    delete m_outfile;
    m_outfile = nullptr;
  }
}

int TpcPolyClusterDisplay::Init(PHCompositeNode*)
{
  m_outfile = new TFile(m_outfilename.c_str(), "RECREATE");
  if (!m_outfile || m_outfile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot open output file " << m_outfilename << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }
  m_outfile->mkdir("events");
  return Fun4AllReturnCodes::EVENT_OK;
}

bool TpcPolyClusterDisplay::get_nodes(PHCompositeNode* topNode)
{
  m_clusterTracks = findNode::getClass<TpcPolyClusterTrackContainer>(topNode, m_clusterNodeName.c_str());
  if (!m_clusterTracks)
  {
    std::cerr << Name() << " - missing " << m_clusterNodeName << std::endl;
    return false;
  }
  m_finalTracks = findNode::getClass<FinalTrackContainer>(topNode, m_finalTrackNodeName.c_str());
  if (!m_finalTracks)
  {
    std::cerr << Name() << " - missing " << m_finalTrackNodeName << ", drawing clusters without fitted lines" << std::endl;
  }
  m_finalTrackVertices = findNode::getClass<FinalTrackVertexContainer>(topNode, m_finalTrackVertexNodeName.c_str());
  if (!m_finalTrackVertices)
  {
    std::cerr << Name() << " - missing " << m_finalTrackVertexNodeName
              << ", drawing without final-track PCA/collision vertices" << std::endl;
  }
  return true;
}

int TpcPolyClusterDisplay::process_event(PHCompositeNode* topNode)
{
  ++m_evt;
  if (!get_nodes(topNode)) return Fun4AllReturnCodes::EVENT_OK;
  if (!m_outfile || m_eventsSaved >= m_maxEventDisplays) return Fun4AllReturnCodes::EVENT_OK;

  TDirectory* eventsTop = m_outfile->GetDirectory("events");
  if (!eventsTop) eventsTop = m_outfile->mkdir("events");
  eventsTop->cd();

  TDirectory* eventDir = eventsTop->mkdir(Form("event_%06u", m_evt));
  if (!eventDir) return Fun4AllReturnCodes::EVENT_OK;
  eventDir->cd();

  TH3D* h3 = new TH3D(Form("h3_evt%06u_tpcpolycluster_centroids_bothsides", m_evt),
                      Form("event %u both sides TpcPolyCluster centroids;z [cm];x [cm];y [cm]", m_evt),
                      204, m_zmin, m_zmax,
                      170, -m_xymax, m_xymax,
                      170, -m_xymax, m_xymax);
  h3->SetStats(0);
  h3->SetDirectory(nullptr);

  
  std::vector<TPolyMarker3D*> markers;
  std::vector<TPolyMarker3D*> pca_markers;
  std::vector<TPolyMarker3D*> collision_vertex_markers;
  std::vector<TPolyLine3D*> lines;
  std::map<unsigned int, const FinalTrackVertex*> final_track_vertices_by_full_track_id;
  std::map<unsigned int, const TpcPolyClusterTrack*> cluster_tracks_by_full_track_id;
  const unsigned int nvertices = m_finalTrackVertices ? m_finalTrackVertices->size() : 0;
  for (unsigned int ivtx = 0; ivtx < nvertices; ++ivtx)
  {
    const FinalTrackVertex* vtx = m_finalTrackVertices->get_vertex(ivtx);
    if (!vtx) continue;
    final_track_vertices_by_full_track_id[vtx->get_source_full_track_id()] = vtx;
  }

  const bool applyTrackVertexZRange = m_finalTrackVertices != nullptr;
  const unsigned int ncluster_tracks = m_clusterTracks->size();
  unsigned int nclusters = 0;
  unsigned int nclustersPlotted = 0;
  for (unsigned int itrack = 0; itrack < ncluster_tracks; ++itrack)
  {
    const TpcPolyClusterTrack* trk = m_clusterTracks->get_track(itrack);
    if (!trk || !trk->isValid()) continue;
    const auto vertex_iter = final_track_vertices_by_full_track_id.find(trk->get_source_full_track_id());
    if (applyTrackVertexZRange &&
        (vertex_iter == final_track_vertices_by_full_track_id.end() ||
         !track_vertex_z_selected(vertex_iter->second, m_trackVertexZMin, m_trackVertexZMax)))
    {
      continue;
    }
    cluster_tracks_by_full_track_id[trk->get_source_full_track_id()] = trk;

    for (unsigned int icluster = 0; icluster < trk->size_clusters(); ++icluster)
    {
      ++nclusters;
      const double x = trk->get_cluster_x(icluster);
      const double y = trk->get_cluster_y(icluster);
      const double z = trk->get_cluster_z(icluster);
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) continue;
      if (z < m_zmin || z > m_zmax) continue;
      markers.push_back(make_cluster_marker(z, x, y, cluster_color(itrack)));

      ++nclustersPlotted;
    }
  }

  const unsigned int nfinal_tracks = m_finalTracks ? m_finalTracks->size() : 0;
  for (unsigned int ifinal = 0; ifinal < nfinal_tracks; ++ifinal)
  {
    const FinalTrack* trk = m_finalTracks->get_track(ifinal);
    if (!trk) continue;

    const auto cluster_iter = cluster_tracks_by_full_track_id.find(trk->get_source_full_track_id());
    if (cluster_iter == cluster_tracks_by_full_track_id.end()) continue;

    double line_zmin = m_zmin;
    double line_zmax = m_zmax;
    if (!cluster_track_z_range(cluster_iter->second, m_zmin, m_zmax, line_zmin, line_zmax)) continue;

    const double forward_residual2 = cluster_line_residual2(trk, cluster_iter->second, m_magneticFieldTesla, 1.0);
    const double reverse_residual2 = cluster_line_residual2(trk, cluster_iter->second, m_magneticFieldTesla, -1.0);
    const double arc_direction = forward_residual2 <= reverse_residual2 ? 1.0 : -1.0;

    TPolyLine3D* line = make_final_track_line(trk, line_zmin, line_zmax, m_xymax,
                                               m_magneticFieldTesla, arc_direction,
                                               cluster_color(cluster_iter->second->get_track_id()));
    if (line) lines.push_back(line);
  }

  for (unsigned int ivtx = 0; ivtx < nvertices; ++ivtx)
  {
    const FinalTrackVertex* vtx = m_finalTrackVertices->get_vertex(ivtx);
    if (!vtx) continue;
    if (!track_vertex_z_selected(vtx, m_trackVertexZMin, m_trackVertexZMax)) continue;

    int color = cluster_color(vtx->get_track_id());
    const auto cluster_iter = cluster_tracks_by_full_track_id.find(vtx->get_source_full_track_id());
    if (cluster_iter != cluster_tracks_by_full_track_id.end())
    {
      color = cluster_color(cluster_iter->second->get_track_id());
    }

    TPolyMarker3D* marker = make_pca_marker(vtx, m_zmin, m_zmax, m_xymax, color);
    if (marker) pca_markers.push_back(marker);
  }

  const unsigned int ncollision_vertices = (m_finalTrackVertices && m_finalTrackVertices->get_collision_vertex_valid())
                                           ? m_finalTrackVertices->get_collision_vertex_count()
                                           : 0;
  for (unsigned int ivtx = 0; ivtx < ncollision_vertices; ++ivtx)
  {
    TPolyMarker3D* marker = make_collision_vertex_marker(m_finalTrackVertices->get_collision_x(ivtx),
                                                         m_finalTrackVertices->get_collision_y(ivtx),
                                                         m_finalTrackVertices->get_collision_z(ivtx),
                                                         m_zmin, m_zmax, m_xymax);
    if (marker) collision_vertex_markers.push_back(marker);
  }

  TCanvas* c3 = new TCanvas(Form("c3_evt%06u_tpcpolycluster_z_x_y_centroids_bothsides", m_evt),
                            Form("event %u both sides TpcPolyCluster centroids", m_evt),
                            1200, 900);
  h3->Draw();
  for (TPolyLine3D* line : lines)
  {
    if (line) line->Draw("same");
  }
  for (TPolyMarker3D* marker : markers)
  {
    if (marker) marker->Draw("same");
  }
  for (TPolyMarker3D* marker : pca_markers)
  {
    if (marker) marker->Draw("same");
  }
  for (TPolyMarker3D* marker : collision_vertex_markers)
  {
    if (marker) marker->Draw("same");
  }
  c3->Modified();
  c3->Update();
  c3->Write();

  std::cout << Name() << " - saved event " << m_evt
            << " cluster_tracks=" << ncluster_tracks
            << " layer_clusters=" << nclusters
            << " plotted=" << nclustersPlotted
            << " final_tracks=" << nfinal_tracks
            << " lines=" << lines.size()
            << " final_track_vertices=" << nvertices
            << " pca_markers=" << pca_markers.size()
            << " collision_vertices=" << ncollision_vertices
            << " collision_markers=" << collision_vertex_markers.size() << std::endl;

  ++m_eventsSaved;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPolyClusterDisplay::End(PHCompositeNode*)
{
  if (m_outfile)
  {
    m_outfile->Close();
    delete m_outfile;
    m_outfile = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
