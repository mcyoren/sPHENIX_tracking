#include "TpcPolyTrackDisplay.h"

#include "TpcPolyTrack.h"
#include "TpcPolyTrackContainer.h"

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
#include <vector>


namespace
{
  int track_color(const unsigned int itrk)
  {
    static const int colors[] = {
      kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kCyan + 2,
      kOrange + 7, kViolet + 1, kAzure + 1, kPink + 7, kTeal + 3
    };
    return colors[itrk % (sizeof(colors) / sizeof(colors[0]))];
  }


  bool choose_perigee_branch(const TpcPolyTrack* trk,
                             const double curvature,
                             const double d0,
                             const double phi0,
                             double& px,
                             double& py,
                             double& xc,
                             double& yc)
  {
    if (!trk) return false;
    if (!std::isfinite(curvature) || std::fabs(curvature) < 1.0e-12) return false;
    if (!std::isfinite(d0) || !std::isfinite(phi0)) return false;

    const double sign = curvature > 0.0 ? 1.0 : -1.0;
    const double R = 1.0 / std::fabs(curvature);

    double best_score = std::numeric_limits<double>::max();
    for (int iphi = 0; iphi < 2; ++iphi)
    {
      const double phi_try = phi0 + static_cast<double>(iphi) * M_PI;
      const double tx = std::cos(phi_try);
      const double ty = std::sin(phi_try);
      for (int id0 = 0; id0 < 2; ++id0)
      {
        const double d0_try = (id0 == 0) ? d0 : -d0;
        const double px_try = -d0_try * ty;
        const double py_try =  d0_try * tx;
        const double rhatx = sign * ty;
        const double rhaty = -sign * tx;
        const double xc_try = px_try - R * rhatx;
        const double yc_try = py_try - R * rhaty;

        double score = 0.0;
        unsigned int nscore = 0;
        for (unsigned int ih = 0; ih < trk->size_hits(); ++ih)
        {
          const double x = trk->get_hit_x(ih);
          const double y = trk->get_hit_y(ih);
          if (!std::isfinite(x) || !std::isfinite(y)) continue;
          score += std::fabs(std::hypot(x - xc_try, y - yc_try) - R);
          ++nscore;
        }
        if (nscore == 0) continue;
        score /= static_cast<double>(nscore);
        if (score < best_score)
        {
          best_score = score;
          px = px_try;
          py = py_try;
          xc = xc_try;
          yc = yc_try;
        }
      }
    }

    return best_score != std::numeric_limits<double>::max();
  }

  bool update_s_range_from_hit(const TpcPolyTrack* trk,
                               const double xc,
                               const double yc,
                               const double R,
                               const double phi_perigee,
                               const double z0,
                               const double dzds,
                               bool& have,
                               double& smin,
                               double& smax)
  {
    if (!trk) return false;

    for (unsigned int ih = 0; ih < trk->size_hits(); ++ih)
    {
      const double x = trk->get_hit_x(ih);
      const double y = trk->get_hit_y(ih);
      const double z = trk->get_hit_z(ih);
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) continue;

      const double angle0 = std::atan2(y - yc, x - xc);
      double best_s = 0.0;
      double best_resid = std::numeric_limits<double>::max();
      for (int k = -4; k <= 4; ++k)
      {
        const double dangle = angle0 - phi_perigee + 2.0 * M_PI * static_cast<double>(k);
        const double s = R * dangle;
        const double zfit = z0 + dzds * s;
        const double resid = std::fabs(zfit - z);
        if (std::isfinite(s) && std::isfinite(resid) && resid < best_resid)
        {
          best_resid = resid;
          best_s = s;
        }
      }

      if (best_resid == std::numeric_limits<double>::max()) continue;
      if (!have)
      {
        smin = best_s;
        smax = best_s;
        have = true;
      }
      else
      {
        smin = std::min(smin, best_s);
        smax = std::max(smax, best_s);
      }
    }

    return have && smax > smin;
  }

  TPolyLine3D* make_trajectory(const TpcPolyTrack* trk,
                               const double zmin,
                               const double zmax,
                               const int color)
  {
    if (!trk) return nullptr;

    const double curvature = trk->get_curvature();
    const double d0 = trk->get_d0();
    const double z0 = trk->get_z0();
    const double phi0 = trk->get_phi0();
    const double theta = trk->get_theta();
    if (!std::isfinite(curvature) || !std::isfinite(d0) ||
        !std::isfinite(z0) || !std::isfinite(phi0) || !std::isfinite(theta))
    {
      return nullptr;
    }
    if (std::fabs(curvature) < 1.0e-12) return nullptr;

    const double R = 1.0 / std::fabs(curvature);
    double px = 0.0;
    double py = 0.0;
    double xc = 0.0;
    double yc = 0.0;
    if (!choose_perigee_branch(trk, curvature, d0, phi0, px, py, xc, yc)) return nullptr;
    const double phi_perigee = std::atan2(py - yc, px - xc);
    const double dzds = 1.0 / std::tan(theta);

    double smin = 0.0;
    double smax = 0.0;
    bool have_s_range = false;
    if (!update_s_range_from_hit(trk, xc, yc, R, phi_perigee, z0, dzds, have_s_range, smin, smax)) return nullptr;

    std::vector<double> zs;
    std::vector<double> xs;
    std::vector<double> ys;
    const unsigned int nsteps = 40;
    zs.reserve(nsteps + 1);
    xs.reserve(nsteps + 1);
    ys.reserve(nsteps + 1);

    for (unsigned int istep = 0; istep <= nsteps; ++istep)
    {
      const double f = static_cast<double>(istep) / static_cast<double>(nsteps);
      const double s = smin + f * (smax - smin);
      const double angle = phi_perigee + s / R;
      const double x = xc + R * std::cos(angle);
      const double y = yc + R * std::sin(angle);
      const double z = z0 + dzds * s;
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) continue;
      if (z < zmin || z > zmax) continue;
      zs.push_back(z);
      xs.push_back(x);
      ys.push_back(y);
    }

    if (zs.size() < 2) return nullptr;

    TPolyLine3D* line = new TPolyLine3D(static_cast<int>(zs.size()));
    for (unsigned int i = 0; i < zs.size(); ++i)
    {
      line->SetPoint(static_cast<int>(i), zs[i], xs[i], ys[i]);
    }
    line->SetLineColor(color);
    line->SetLineWidth(3);
    line->SetLineStyle(1);
    return line;
  }

  TPolyMarker3D* make_perigee_marker(const TpcPolyTrack* trk,
                                     const double zmin,
                                     const double zmax,
                                     const int color)
  {
    if (!trk) return nullptr;

    const double curvature = trk->get_curvature();
    const double d0 = trk->get_d0();
    const double z0 = trk->get_z0();
    const double phi0 = trk->get_phi0();
    if (!std::isfinite(curvature) || !std::isfinite(d0) ||
        !std::isfinite(z0) || !std::isfinite(phi0))
    {
      return nullptr;
    }
    if (z0 < zmin || z0 > zmax) return nullptr;

    double px = 0.0;
    double py = 0.0;
    double xc = 0.0;
    double yc = 0.0;
    if (!choose_perigee_branch(trk, curvature, d0, phi0, px, py, xc, yc)) return nullptr;

    TPolyMarker3D* marker = new TPolyMarker3D(1);
    marker->SetPoint(0, z0, px, py);
    marker->SetMarkerColor(color);
    marker->SetMarkerStyle(20);
    marker->SetMarkerSize(1.2);
    return marker;
  }

}

TpcPolyTrackDisplay::TpcPolyTrackDisplay(const std::string& name,
                                         const std::string& outfilename,
                                         const std::string& trackNodeName,
                                         const unsigned int maxEventDisplays)
  : SubsysReco(name)
  , m_outfilename(outfilename)
  , m_trackNodeName(trackNodeName)
  , m_maxEventDisplays(maxEventDisplays)
  , m_evt(0)
  , m_eventsSaved(0)
  , m_zmin(-102.0)
  , m_zmax(102.0)
  , m_xymax(85.0)
  , m_lowPtMax(0.1)
  , m_highPtMin(0.3)
  , m_highPtMaxAbsD0(3.0)
  , m_magneticFieldTesla(1.4)
  , m_outfile(nullptr)
  , m_tracks(nullptr)
{
}

TpcPolyTrackDisplay::~TpcPolyTrackDisplay()
{
  if (m_outfile)
  {
    delete m_outfile;
    m_outfile = nullptr;
  }
}

int TpcPolyTrackDisplay::Init(PHCompositeNode*)
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

bool TpcPolyTrackDisplay::get_nodes(PHCompositeNode* topNode)
{
  m_tracks = findNode::getClass<TpcPolyTrackContainer>(topNode, m_trackNodeName.c_str());
  if (!m_tracks)
  {
    std::cerr << Name() << " - missing " << m_trackNodeName << std::endl;
    return false;
  }
  return true;
}

int TpcPolyTrackDisplay::process_event(PHCompositeNode* topNode)
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

  TH3D* h3 = new TH3D(Form("h3_evt%06u_tpcpolytrack_hits_bothsides", m_evt),
                      Form("event %u both sides TpcPolyTrack hits;z [cm];x [cm];y [cm]", m_evt),
                      204, m_zmin, m_zmax,
                      170, -m_xymax, m_xymax,
                      170, -m_xymax, m_xymax);
  h3->SetStats(0);
  h3->SetDirectory(nullptr);

  TH3D* h3_lowpt = new TH3D(Form("h3_evt%06u_tpcpolytrack_lowpt_hits_bothsides", m_evt),
                            Form("event %u both sides TpcPolyTrack p_{T}<%.3g GeV hits;z [cm];x [cm];y [cm]", m_evt, m_lowPtMax),
                            204, m_zmin, m_zmax,
                            170, -m_xymax, m_xymax,
                            170, -m_xymax, m_xymax);
  h3_lowpt->SetStats(0);
  h3_lowpt->SetDirectory(nullptr);

  TH3D* h3_highpt = new TH3D(Form("h3_evt%06u_tpcpolytrack_highpt_d0_hits_bothsides", m_evt),
                             Form("event %u both sides TpcPolyTrack p_{T}>%.3g GeV |d_{0}|<%.3g cm hits;z [cm];x [cm];y [cm]", m_evt, m_highPtMin, m_highPtMaxAbsD0),
                             204, m_zmin, m_zmax,
                             170, -m_xymax, m_xymax,
                             170, -m_xymax, m_xymax);
  h3_highpt->SetStats(0);
  h3_highpt->SetDirectory(nullptr);

  std::vector<TPolyLine3D*> trajectories;
  std::vector<TPolyMarker3D*> perigee_markers;
  std::vector<TPolyLine3D*> lowpt_trajectories;
  std::vector<TPolyMarker3D*> lowpt_perigee_markers;
  std::vector<TPolyLine3D*> highpt_trajectories;
  std::vector<TPolyMarker3D*> highpt_perigee_markers;
  const unsigned int ntracks = m_tracks->size();
  unsigned int ntracksPlotted = 0;
  unsigned int ntracksLowPt = 0;
  unsigned int ntracksHighPt = 0;
  for (unsigned int itrk = 0; itrk < ntracks; ++itrk)
  {
    const TpcPolyTrack* trk = m_tracks->get_track(itrk);
    if (!trk) continue;

    const bool fit_ok = trk->get_fit_status() != 0;
    const double curvature = trk->get_curvature();
    const double pt = (fit_ok && std::isfinite(curvature) && std::fabs(curvature) > 1.0e-12)
                      ? 0.003 * m_magneticFieldTesla / std::fabs(curvature)
                      : -1.0;
    const double d0 = trk->get_d0();
    const bool is_lowpt = fit_ok && std::isfinite(pt) && pt >= 0.0 && pt < m_lowPtMax;
    const bool is_highpt = fit_ok && std::isfinite(pt) && pt > m_highPtMin &&
                           std::isfinite(d0) && std::fabs(d0) < m_highPtMaxAbsD0;

    bool has_plotted_hit = false;
    for (unsigned int ih = 0; ih < trk->size_hits(); ++ih)
    {
      const double x = trk->get_hit_x(ih);
      const double y = trk->get_hit_y(ih);
      const double z = trk->get_hit_z(ih);
      if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) continue;
      if (z < m_zmin || z > m_zmax) continue;
      if (fit_ok) h3->Fill(z, x, y);
      if (is_lowpt) h3_lowpt->Fill(z, x, y);
      if (is_highpt) h3_highpt->Fill(z, x, y);
      has_plotted_hit = true;
    }

    if (fit_ok)
    {
      TPolyLine3D* trajectory = make_trajectory(trk, m_zmin, m_zmax, track_color(itrk));
      if (trajectory)
      {
        trajectories.push_back(trajectory);
        TPolyMarker3D* marker = make_perigee_marker(trk, m_zmin, m_zmax, track_color(itrk));
        if (marker) perigee_markers.push_back(marker);
        if (is_lowpt)
        {
          TPolyLine3D* lowpt_trajectory = make_trajectory(trk, m_zmin, m_zmax, track_color(itrk));
          if (lowpt_trajectory) lowpt_trajectories.push_back(lowpt_trajectory);
          TPolyMarker3D* lowpt_marker = make_perigee_marker(trk, m_zmin, m_zmax, track_color(itrk));
          if (lowpt_marker) lowpt_perigee_markers.push_back(lowpt_marker);
        }
        if (is_highpt)
        {
          TPolyLine3D* highpt_trajectory = make_trajectory(trk, m_zmin, m_zmax, track_color(itrk));
          if (highpt_trajectory) highpt_trajectories.push_back(highpt_trajectory);
          TPolyMarker3D* highpt_marker = make_perigee_marker(trk, m_zmin, m_zmax, track_color(itrk));
          if (highpt_marker) highpt_perigee_markers.push_back(highpt_marker);
        }
      }
      if (has_plotted_hit) ++ntracksPlotted;
      if (is_lowpt && has_plotted_hit) ++ntracksLowPt;
      if (is_highpt && has_plotted_hit) ++ntracksHighPt;
    }
  }

  TCanvas* c3 = new TCanvas(Form("c3_evt%06u_tpcpolytrack_z_x_y_hits_trajectory_bothsides", m_evt),
                            Form("event %u both sides TpcPolyTrack hits and trajectory", m_evt),
                            1200, 900);
  h3->Draw("BOX2Z");
  for (TPolyLine3D* trajectory : trajectories)
  {
    if (trajectory) trajectory->Draw("same");
  }
  for (TPolyMarker3D* marker : perigee_markers)
  {
    if (marker) marker->Draw("same");
  }
  c3->Modified();
  c3->Update();
  c3->Write();

  TCanvas* c3_lowpt = new TCanvas(Form("c3_evt%06u_tpcpolytrack_z_x_y_lowpt_hits_trajectory_bothsides", m_evt),
                                  Form("event %u both sides TpcPolyTrack p_{T}<%.3g GeV", m_evt, m_lowPtMax),
                                  1200, 900);
  h3_lowpt->Draw("BOX2Z");
  for (TPolyLine3D* trajectory : lowpt_trajectories)
  {
    if (trajectory) trajectory->Draw("same");
  }
  for (TPolyMarker3D* marker : lowpt_perigee_markers)
  {
    if (marker) marker->Draw("same");
  }
  c3_lowpt->Modified();
  c3_lowpt->Update();
  c3_lowpt->Write();

  TCanvas* c3_highpt = new TCanvas(Form("c3_evt%06u_tpcpolytrack_z_x_y_highpt_d0_hits_trajectory_bothsides", m_evt),
                                   Form("event %u both sides TpcPolyTrack p_{T}>%.3g GeV |d_{0}|<%.3g cm", m_evt, m_highPtMin, m_highPtMaxAbsD0),
                                   1200, 900);
  h3_highpt->Draw("BOX2Z");
  for (TPolyLine3D* trajectory : highpt_trajectories)
  {
    if (trajectory) trajectory->Draw("same");
  }
  for (TPolyMarker3D* marker : highpt_perigee_markers)
  {
    if (marker) marker->Draw("same");
  }
  c3_highpt->Modified();
  c3_highpt->Update();
  c3_highpt->Write();

  std::cout << Name() << " - saved event " << m_evt
            << " tracks=" << ntracks
            << " plotted=" << ntracksPlotted
            << " lowpt=" << ntracksLowPt
            << " highpt=" << ntracksHighPt << std::endl;

  ++m_eventsSaved;
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcPolyTrackDisplay::End(PHCompositeNode*)
{
  if (m_outfile)
  {
    m_outfile->Close();
    delete m_outfile;
    m_outfile = nullptr;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
