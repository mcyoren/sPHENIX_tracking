#include "TpcDstV0Finder.h"

#include <inmoduletracks/FinalTrack.h>
#include <inmoduletracks/FinalTrackContainer.h>
#include <inmoduletracks/FinalTrackVertexContainer.h>

#include <ffaobjects/EventHeader.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <utility>

namespace
{
  constexpr double kPionMass = 0.13957039;
  constexpr double kProtonMass = 0.938272088;
  constexpr double kTiny = 1.0e-12;

  double clamp_value(double x, double lo, double hi)
  {
    return std::max(lo, std::min(x, hi));
  }
}

TpcDstV0Finder::TpcDstV0Finder(const std::string& name, const std::string& output)
  : SubsysReco(name)
  , m_output_file(output)
{
}

void TpcDstV0Finder::set_primary_vertex(double x, double y, double z)
{
  m_fixed_primary = {x, y, z};
}

void TpcDstV0Finder::set_search_range(double backward_cm, double forward_cm)
{
  m_search_backward_cm = std::max(0.0, backward_cm);
  m_search_forward_cm = std::max(0.0, forward_cm);
}

int TpcDstV0Finder::Init(PHCompositeNode*)
{
  m_file = TFile::Open(m_output_file.c_str(), "RECREATE");
  if (!m_file || m_file->IsZombie())
  {
    std::cerr << Name() << ": cannot create " << m_output_file << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_pair_tree = new TTree("pairTree", "V0 candidates from stored FINALTRACKS states; no track refit");
  if (m_write_track_tree)
  {
    m_track_tree = new TTree("trackTree", "Single FINALTRACKS states used by the V0 finder");
  }
  create_branches();

  if (m_refit_after_vertex)
  {
    std::cout << Name() << ": vertex-constrained refit was requested, but this version intentionally "
              << "keeps the stored DST track state. refit_used will remain 0." << std::endl;
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcDstV0Finder::process_event(PHCompositeNode* topNode)
{
  ++m_events;
  const int run = get_run(topNode);
  const int evt = get_event(topNode);

  auto* track_object = findNode::getClass<PHObject>(topNode, m_track_node);
  auto* tracks = dynamic_cast<FinalTrackContainer*>(track_object);
  if (!tracks)
  {
    std::cerr << Name() << ": missing or wrong-type node " << m_track_node << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  FinalTrackVertexContainer* vertices = nullptr;
  if (m_use_vertex_node)
  {
    auto* vertex_object = findNode::getClass<PHObject>(topNode, m_vertex_node);
    vertices = dynamic_cast<FinalTrackVertexContainer*>(vertex_object);
  }
  const Vec3 primary = choose_primary_vertex(vertices);

  std::vector<State> selected;
  selected.reserve(tracks->size());
  for (unsigned int i = 0; i < tracks->size(); ++i)
  {
    State state;
    if (!make_state(tracks->get_track(i), state)) continue;
    selected.push_back(state);
    ++m_tracks;
    if (m_track_tree) fill_track(state, primary, run, evt);
  }

  for (std::size_t i = 0; i < selected.size(); ++i)
  {
    for (std::size_t j = i + 1; j < selected.size(); ++j)
    {
      ++m_pairs_tested;
      fill_pair(selected[i], selected[j], primary, run, evt);
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int TpcDstV0Finder::End(PHCompositeNode*)
{
  if (m_file)
  {
    m_file->cd();
    if (m_pair_tree) m_pair_tree->Write();
    if (m_track_tree) m_track_tree->Write();
    m_file->Close();
    delete m_file;
    m_file = nullptr;
  }

  std::cout << Name() << ": events=" << m_events
            << " selected_tracks=" << m_tracks
            << " tested_pairs=" << m_pairs_tested
            << " written_pairs=" << m_pairs_written << std::endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

bool TpcDstV0Finder::make_state(const FinalTrack* track, State& out) const
{
  if (!track || !track->isValid() || track->get_fit_status() == 0) return false;

  out.id = track->get_track_id();
  out.source_id = track->get_source_full_track_id();
  out.nclusters = track->get_nclusters();
  out.fit_status = track->get_fit_status();
  out.charge = (track->get_charge() > 0.0) - (track->get_charge() < 0.0);
  out.chi2 = track->get_chi2();
  out.ndf = track->get_ndf();
  out.dedx = track->get_dedx();
  out.r0 = {track->get_x(), track->get_y(), track->get_z()};
  out.p0 = {track->get_px(), track->get_py(), track->get_pz()};
  out.pt = std::hypot(out.p0.x, out.p0.y);
  out.p = norm(out.p0);

  if (out.charge == 0 || out.nclusters < m_min_clusters || out.pt < m_min_pt || out.p <= 0.0) return false;
  if (!std::isfinite(out.r0.x) || !std::isfinite(out.r0.y) || !std::isfinite(out.r0.z) ||
      !std::isfinite(out.p0.x) || !std::isfinite(out.p0.y) || !std::isfinite(out.p0.z)) return false;
  if (m_max_track_chi2_ndf >= 0.0 && out.ndf > 0.0 && out.chi2 / out.ndf > m_max_track_chi2_ndf) return false;

  out.phi0 = std::atan2(out.p0.y, out.p0.x);
  out.tanl = out.p0.z / out.pt;
  out.omega = 0.003 * static_cast<double>(out.charge) * m_bfield_t / out.pt;
  return std::isfinite(out.omega);
}

TpcDstV0Finder::Vec3 TpcDstV0Finder::position(const State& t, double s) const
{
  if (std::abs(t.omega) < 1.0e-10)
  {
    return {t.r0.x + s * std::cos(t.phi0),
            t.r0.y + s * std::sin(t.phi0),
            t.r0.z + s * t.tanl};
  }
  const double phi = t.phi0 - t.omega * s;
  return {t.r0.x + (std::sin(t.phi0) - std::sin(phi)) / t.omega,
          t.r0.y + (std::cos(phi) - std::cos(t.phi0)) / t.omega,
          t.r0.z + s * t.tanl};
}

TpcDstV0Finder::Vec3 TpcDstV0Finder::momentum(const State& t, double s) const
{
  const double phi = t.phi0 - t.omega * s;
  return {t.pt * std::cos(phi), t.pt * std::sin(phi), t.p0.z};
}

TpcDstV0Finder::Vec3 TpcDstV0Finder::tangent(const State& t, double s) const
{
  const double phi = t.phi0 - t.omega * s;
  return {std::cos(phi), std::sin(phi), t.tanl};
}

TpcDstV0Finder::PcaResult TpcDstV0Finder::refine_pair(
    const State& a, const State& b, double s1, double s2) const
{
  const double lo = -m_search_backward_cm;
  const double hi = m_search_forward_cm;
  PcaResult best;

  for (unsigned int iteration = 0; iteration < 30; ++iteration)
  {
    const Vec3 x1 = position(a, s1);
    const Vec3 x2 = position(b, s2);
    const Vec3 u = tangent(a, s1);
    const Vec3 v = tangent(b, s2);
    const Vec3 w = sub(x1, x2);

    const double uu = dot(u, u);
    const double vv = dot(v, v);
    const double uv = dot(u, v);
    const double uw = dot(u, w);
    const double vw = dot(v, w);
    const double det = uu * vv - uv * uv;
    if (std::abs(det) < 1.0e-12) break;

    double ds1 = (uv * vw - vv * uw) / det;
    double ds2 = (uu * vw - uv * uw) / det;
    ds1 = clamp_value(ds1, -10.0, 10.0);
    ds2 = clamp_value(ds2, -10.0, 10.0);

    const double next1 = clamp_value(s1 + ds1, lo, hi);
    const double next2 = clamp_value(s2 + ds2, lo, hi);
    if (std::abs(next1 - s1) + std::abs(next2 - s2) < 1.0e-5)
    {
      s1 = next1;
      s2 = next2;
      break;
    }
    s1 = next1;
    s2 = next2;
  }

  best.s1 = s1;
  best.s2 = s2;
  best.x1 = position(a, s1);
  best.x2 = position(b, s2);
  best.dca = distance(best.x1, best.x2);
  best.valid = std::isfinite(best.dca);
  return best;
}

TpcDstV0Finder::PcaResult TpcDstV0Finder::pair_pca(const State& a, const State& b) const
{
  const unsigned int steps = std::max(4u, m_coarse_steps);
  const double lo = -m_search_backward_cm;
  const double hi = m_search_forward_cm;
  PcaResult best;
  best.dca = std::numeric_limits<double>::infinity();

  // Coarse global scan avoids choosing the wrong turn/branch of a curved trajectory.
  for (unsigned int i = 0; i <= steps; ++i)
  {
    const double s1 = lo + (hi - lo) * static_cast<double>(i) / steps;
    const Vec3 x1 = position(a, s1);
    for (unsigned int j = 0; j <= steps; ++j)
    {
      const double s2 = lo + (hi - lo) * static_cast<double>(j) / steps;
      const double d = distance(x1, position(b, s2));
      if (d < best.dca)
      {
        best.dca = d;
        best.s1 = s1;
        best.s2 = s2;
      }
    }
  }

  // Refine the best seed plus nearby coarse seeds. This is still propagation, not a refit.
  const double ds = (hi - lo) / steps;
  for (int di = -1; di <= 1; ++di)
  {
    for (int dj = -1; dj <= 1; ++dj)
    {
      const PcaResult candidate = refine_pair(
          a, b, clamp_value(best.s1 + di * ds, lo, hi), clamp_value(best.s2 + dj * ds, lo, hi));
      if (candidate.valid && candidate.dca < best.dca) best = candidate;
    }
  }
  best.valid = std::isfinite(best.dca);
  return best;
}

std::pair<double, double> TpcDstV0Finder::dca_to_vertex(const State& t, const Vec3& vertex) const
{
  // One-dimensional minimization of 3D distance to the chosen primary vertex.
  double best_s = 0.0;
  double best_d2 = std::numeric_limits<double>::infinity();
  const unsigned int steps = std::max(20u, m_coarse_steps);
  const double lo = -m_search_backward_cm;
  const double hi = m_search_forward_cm;
  for (unsigned int i = 0; i <= steps; ++i)
  {
    const double s = lo + (hi - lo) * static_cast<double>(i) / steps;
    const Vec3 d = sub(position(t, s), vertex);
    const double d2 = dot(d, d);
    if (d2 < best_d2) { best_d2 = d2; best_s = s; }
  }
  double step = (hi - lo) / steps;
  for (unsigned int iteration = 0; iteration < 30; ++iteration)
  {
    const Vec3 x = position(t, best_s);
    const Vec3 u = tangent(t, best_s);
    const double uu = dot(u, u);
    if (uu < kTiny) break;
    const double ds = clamp_value(-dot(sub(x, vertex), u) / uu, -step, step);
    const double next = clamp_value(best_s + ds, lo, hi);
    if (std::abs(next - best_s) < 1.0e-6) break;
    best_s = next;
    step *= 0.8;
  }

  const Vec3 x = position(t, best_s);
  const Vec3 p = momentum(t, best_s);
  const double dxy = std::hypot(x.x - vertex.x, x.y - vertex.y);
  const double sign = ((x.x - vertex.x) * p.y - (x.y - vertex.y) * p.x) >= 0.0 ? 1.0 : -1.0;
  return {sign * dxy, x.z - vertex.z};
}

TpcDstV0Finder::Vec3 TpcDstV0Finder::choose_primary_vertex(const FinalTrackVertexContainer* vertices) const
{
  if (vertices && vertices->get_collision_vertex_valid() && vertices->get_collision_vertex_count() > 0)
  {
    return {vertices->get_collision_x(0), vertices->get_collision_y(0), vertices->get_collision_z(0)};
  }
  return m_fixed_primary;
}

bool TpcDstV0Finder::fill_pair(
    const State& a, const State& b, const Vec3& primary, int run, int evt)
{
  if (!m_write_same_sign && a.charge == b.charge) return false;

  const PcaResult pca = pair_pca(a, b);
  if (!pca.valid || pca.dca > m_max_pair_dca) return false;

  const Vec3 vertex = scale(add(pca.x1, pca.x2), 0.5);
  const Vec3 mom1 = momentum(a, pca.s1);
  const Vec3 mom2 = momentum(b, pca.s2);
  const Vec3 v0mom = add(mom1, mom2);
  const Vec3 flight = sub(vertex, primary);
  const double lproj = norm(flight);
  const double dira = cosine(flight, v0mom);
  const double radius = std::hypot(vertex.x - primary.x, vertex.y - primary.y);

  if (!std::isfinite(dira) || lproj < m_min_lproj || dira < m_min_dira ||
      radius < m_min_decay_radius || std::abs(vertex.z) > m_max_abs_decay_z) return false;

  const Vec3& pplus = a.charge > 0 ? mom1 : mom2;
  const Vec3& pminus = a.charge > 0 ? mom2 : mom1;
  double alpha = 0.0;
  double qt = 0.0;
  if (!armenteros(pplus, pminus, alpha, qt) || std::abs(alpha) > m_max_abs_alpha) return false;

  const auto dca1 = dca_to_vertex(a, primary);
  const auto dca2 = dca_to_vertex(b, primary);

  m_pair = PairRow{};
  m_pair.run = run;
  m_pair.evt = evt;
  m_pair.track_id1 = a.id;
  m_pair.track_id2 = b.id;
  m_pair.source_id1 = a.source_id;
  m_pair.source_id2 = b.source_id;
  m_pair.charge1 = static_cast<short>(a.charge);
  m_pair.charge2 = static_cast<short>(b.charge);
  m_pair.nclusters1 = static_cast<unsigned short>(std::min(a.nclusters, 65535u));
  m_pair.nclusters2 = static_cast<unsigned short>(std::min(b.nclusters, 65535u));
  m_pair.quality1 = f(a.ndf > 0.0 ? a.chi2 / a.ndf : -1.0);
  m_pair.quality2 = f(b.ndf > 0.0 ? b.chi2 / b.ndf : -1.0);
  m_pair.dedx_1 = f(a.dedx);
  m_pair.dedx_2 = f(b.dedx);
  m_pair.ref_x1=f(a.r0.x); m_pair.ref_y1=f(a.r0.y); m_pair.ref_z1=f(a.r0.z);
  m_pair.ref_x2=f(b.r0.x); m_pair.ref_y2=f(b.r0.y); m_pair.ref_z2=f(b.r0.z);
  m_pair.ref_px1=f(a.p0.x); m_pair.ref_py1=f(a.p0.y); m_pair.ref_pz1=f(a.p0.z);
  m_pair.ref_px2=f(b.p0.x); m_pair.ref_py2=f(b.p0.y); m_pair.ref_pz2=f(b.p0.z);
  m_pair.pca1_x=f(pca.x1.x); m_pair.pca1_y=f(pca.x1.y); m_pair.pca1_z=f(pca.x1.z);
  m_pair.pca2_x=f(pca.x2.x); m_pair.pca2_y=f(pca.x2.y); m_pair.pca2_z=f(pca.x2.z);
  m_pair.pca_x=f(vertex.x); m_pair.pca_y=f(vertex.y); m_pair.pca_z=f(vertex.z);
  m_pair.px1=f(mom1.x); m_pair.py1=f(mom1.y); m_pair.pz1=f(mom1.z);
  m_pair.px2=f(mom2.x); m_pair.py2=f(mom2.y); m_pair.pz2=f(mom2.z);
  m_pair.v0_px=f(v0mom.x); m_pair.v0_py=f(v0mom.y); m_pair.v0_pz=f(v0mom.z);
  m_pair.v0_pt=f(std::hypot(v0mom.x, v0mom.y));
  m_pair.pairDCA=f(pca.dca);
  m_pair.dca_xy1=f(dca1.first); m_pair.dca_z1=f(dca1.second);
  m_pair.dca_xy2=f(dca2.first); m_pair.dca_z2=f(dca2.second);
  m_pair.decay_radius=f(radius);
  m_pair.Lproj=f(lproj);
  m_pair.cosThetaReco=f(dira);
  m_pair.alpha=f(alpha);
  m_pair.qT=f(qt);
  m_pair.mass_Kshort=f(invariant_mass(mom1, kPionMass, mom2, kPionMass));
  m_pair.mass_Lambda=f(invariant_mass(pplus, kProtonMass, pminus, kPionMass));
  m_pair.mass_AntiLambda=f(invariant_mass(pplus, kPionMass, pminus, kProtonMass));
  m_pair.primary_x=f(primary.x); m_pair.primary_y=f(primary.y); m_pair.primary_z=f(primary.z);
  m_pair.s1=f(pca.s1); m_pair.s2=f(pca.s2);
  m_pair.refit_used=0;

  m_pair_tree->Fill();
  ++m_pairs_written;
  return true;
}

void TpcDstV0Finder::fill_track(const State& t, const Vec3& primary, int run, int evt)
{
  const auto dca = dca_to_vertex(t, primary);
  m_track = TrackRow{};
  m_track.run=run; m_track.evt=evt; m_track.track_id=t.id; m_track.source_id=t.source_id;
  m_track.nclusters=t.nclusters; m_track.fit_status=t.fit_status; m_track.charge=static_cast<short>(t.charge);
  m_track.x=f(t.r0.x); m_track.y=f(t.r0.y); m_track.z=f(t.r0.z);
  m_track.px=f(t.p0.x); m_track.py=f(t.p0.y); m_track.pz=f(t.p0.z);
  m_track.pt=f(t.pt); m_track.p=f(t.p);
  m_track.eta=f(std::asinh(t.p0.z/t.pt));
  m_track.chi2=f(t.chi2); m_track.ndf=f(t.ndf);
  m_track.quality=f(t.ndf>0.0 ? t.chi2/t.ndf : -1.0); m_track.dedx=f(t.dedx);
  m_track.dca_xy=f(dca.first); m_track.dca_z=f(dca.second);
  m_track.primary_x=f(primary.x); m_track.primary_y=f(primary.y); m_track.primary_z=f(primary.z);
  m_track_tree->Fill();
}

void TpcDstV0Finder::create_branches()
{
#define B(tree, row, field) tree->Branch(#field, &row.field)
  B(m_pair_tree,m_pair,run); B(m_pair_tree,m_pair,evt);
  B(m_pair_tree,m_pair,track_id1); B(m_pair_tree,m_pair,track_id2);
  B(m_pair_tree,m_pair,source_id1); B(m_pair_tree,m_pair,source_id2);
  B(m_pair_tree,m_pair,charge1); B(m_pair_tree,m_pair,charge2);
  B(m_pair_tree,m_pair,nclusters1); B(m_pair_tree,m_pair,nclusters2);
  B(m_pair_tree,m_pair,quality1); B(m_pair_tree,m_pair,quality2);
  B(m_pair_tree,m_pair,dedx_1); B(m_pair_tree,m_pair,dedx_2);
  B(m_pair_tree,m_pair,ref_x1); B(m_pair_tree,m_pair,ref_y1); B(m_pair_tree,m_pair,ref_z1);
  B(m_pair_tree,m_pair,ref_x2); B(m_pair_tree,m_pair,ref_y2); B(m_pair_tree,m_pair,ref_z2);
  B(m_pair_tree,m_pair,ref_px1); B(m_pair_tree,m_pair,ref_py1); B(m_pair_tree,m_pair,ref_pz1);
  B(m_pair_tree,m_pair,ref_px2); B(m_pair_tree,m_pair,ref_py2); B(m_pair_tree,m_pair,ref_pz2);
  B(m_pair_tree,m_pair,pca1_x); B(m_pair_tree,m_pair,pca1_y); B(m_pair_tree,m_pair,pca1_z);
  B(m_pair_tree,m_pair,pca2_x); B(m_pair_tree,m_pair,pca2_y); B(m_pair_tree,m_pair,pca2_z);
  B(m_pair_tree,m_pair,pca_x); B(m_pair_tree,m_pair,pca_y); B(m_pair_tree,m_pair,pca_z);
  B(m_pair_tree,m_pair,px1); B(m_pair_tree,m_pair,py1); B(m_pair_tree,m_pair,pz1);
  B(m_pair_tree,m_pair,px2); B(m_pair_tree,m_pair,py2); B(m_pair_tree,m_pair,pz2);
  B(m_pair_tree,m_pair,v0_px); B(m_pair_tree,m_pair,v0_py); B(m_pair_tree,m_pair,v0_pz); B(m_pair_tree,m_pair,v0_pt);
  B(m_pair_tree,m_pair,pairDCA); B(m_pair_tree,m_pair,dca_xy1); B(m_pair_tree,m_pair,dca_z1);
  B(m_pair_tree,m_pair,dca_xy2); B(m_pair_tree,m_pair,dca_z2);
  B(m_pair_tree,m_pair,decay_radius); B(m_pair_tree,m_pair,Lproj); B(m_pair_tree,m_pair,cosThetaReco);
  B(m_pair_tree,m_pair,alpha); B(m_pair_tree,m_pair,qT);
  B(m_pair_tree,m_pair,mass_Kshort); B(m_pair_tree,m_pair,mass_Lambda); B(m_pair_tree,m_pair,mass_AntiLambda);
  B(m_pair_tree,m_pair,primary_x); B(m_pair_tree,m_pair,primary_y); B(m_pair_tree,m_pair,primary_z);
  B(m_pair_tree,m_pair,s1); B(m_pair_tree,m_pair,s2); B(m_pair_tree,m_pair,refit_used);

  if (m_track_tree)
  {
    B(m_track_tree,m_track,run); B(m_track_tree,m_track,evt); B(m_track_tree,m_track,track_id);
    B(m_track_tree,m_track,source_id); B(m_track_tree,m_track,nclusters); B(m_track_tree,m_track,fit_status);
    B(m_track_tree,m_track,charge); B(m_track_tree,m_track,x); B(m_track_tree,m_track,y); B(m_track_tree,m_track,z);
    B(m_track_tree,m_track,px); B(m_track_tree,m_track,py); B(m_track_tree,m_track,pz);
    B(m_track_tree,m_track,pt); B(m_track_tree,m_track,p); B(m_track_tree,m_track,eta);
    B(m_track_tree,m_track,chi2); B(m_track_tree,m_track,ndf); B(m_track_tree,m_track,quality); B(m_track_tree,m_track,dedx);
    B(m_track_tree,m_track,dca_xy); B(m_track_tree,m_track,dca_z);
    B(m_track_tree,m_track,primary_x); B(m_track_tree,m_track,primary_y); B(m_track_tree,m_track,primary_z);
  }
#undef B
}

int TpcDstV0Finder::get_run(PHCompositeNode* topNode) const
{
  if (auto* h = findNode::getClass<EventHeader>(topNode, "EventHeader")) return h->get_RunNumber();
  return 0;
}

int TpcDstV0Finder::get_event(PHCompositeNode* topNode) const
{
  if (auto* h = findNode::getClass<EventHeader>(topNode, "EventHeader")) return h->get_EvtSequence();
  return static_cast<int>(m_events);
}

TpcDstV0Finder::Vec3 TpcDstV0Finder::add(const Vec3& a,const Vec3& b){return {a.x+b.x,a.y+b.y,a.z+b.z};}
TpcDstV0Finder::Vec3 TpcDstV0Finder::sub(const Vec3& a,const Vec3& b){return {a.x-b.x,a.y-b.y,a.z-b.z};}
TpcDstV0Finder::Vec3 TpcDstV0Finder::scale(const Vec3& a,double s){return {s*a.x,s*a.y,s*a.z};}
double TpcDstV0Finder::dot(const Vec3& a,const Vec3& b){return a.x*b.x+a.y*b.y+a.z*b.z;}
double TpcDstV0Finder::norm(const Vec3& a){return std::sqrt(dot(a,a));}
double TpcDstV0Finder::distance(const Vec3& a,const Vec3& b){return norm(sub(a,b));}
double TpcDstV0Finder::cosine(const Vec3& a,const Vec3& b){const double n=norm(a)*norm(b);return n>kTiny?dot(a,b)/n:std::numeric_limits<double>::quiet_NaN();}

double TpcDstV0Finder::invariant_mass(const Vec3& p1,double m1,const Vec3& p2,double m2)
{
  const double e1=std::sqrt(dot(p1,p1)+m1*m1); const double e2=std::sqrt(dot(p2,p2)+m2*m2);
  const Vec3 p=add(p1,p2); const double m2pair=(e1+e2)*(e1+e2)-dot(p,p);
  return std::sqrt(std::max(0.0,m2pair));
}

bool TpcDstV0Finder::armenteros(const Vec3& pp,const Vec3& pm,double& alpha,double& qt)
{
  const Vec3 v=add(pp,pm); const double pv=norm(v); if(pv<kTiny) return false;
  const Vec3 u=scale(v,1.0/pv); const double lp=dot(pp,u); const double lm=dot(pm,u);
  const double den=lp+lm; if(std::abs(den)<kTiny) return false;
  alpha=(lp-lm)/den; const Vec3 transverse=sub(pp,scale(u,lp)); qt=norm(transverse);
  return std::isfinite(alpha)&&std::isfinite(qt);
}

float TpcDstV0Finder::f(double value)
{
  return std::isfinite(value) ? static_cast<float>(value) : std::numeric_limits<float>::quiet_NaN();
}
