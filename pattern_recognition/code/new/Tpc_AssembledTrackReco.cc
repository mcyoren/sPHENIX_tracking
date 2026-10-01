#include "Tpc_AssembledTrackReco.h"

#include "IdealPadMap.h"
#include "Tpc_AssembledTrack.h"
#include "Tpc_AssembledTrackContainer.h"
#include "Tpc_AssembledTrackContainerv1.h"
#include "Tpc_AssembledTrackv1.h"
#include "Tpc_FittingTools.h"

#include "Tpc_ModuleTrack.h"
#include "Tpc_ModuleTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/getClass.h>

#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <iterator>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

namespace
{
  inline double sqr(const double value)
  {
    return value * value;
  }

  double piece_match_strength(unsigned int nblobs)
  {
    constexpr double kReliabilityCap = 12.0;
    return std::max(1.0, std::min(static_cast<double>(nblobs), kReliabilityCap));
  }

  double piece_slope_reliability(unsigned int nblobs_a, unsigned int nblobs_b)
  {
    const double ra = piece_match_strength(nblobs_a);
    const double rb = piece_match_strength(nblobs_b);
    return 2.0 * std::min(ra, rb) / std::max(ra + rb, 1.0e-9);
  }

  double reliability_weighted_coordinate(double xa, double xb,
                                         unsigned int nblobs_a,
                                         unsigned int nblobs_b)
  {
    const double ra = piece_match_strength(nblobs_a);
    const double rb = piece_match_strength(nblobs_b);
    return (rb * xa + ra * xb) / (ra + rb);
  }

  unsigned int minimum_layer_mask_gap(uint64_t a, uint64_t b)
  {
    if ((a & b) != 0U)
    {
      return std::numeric_limits<unsigned int>::max();
    }

    unsigned int best = std::numeric_limits<unsigned int>::max();
    for (unsigned int ia = 0; ia < 48U; ++ia)
    {
      if ((a & (uint64_t{1} << ia)) == 0U) continue;
      for (unsigned int ib = 0; ib < 48U; ++ib)
      {
        if ((b & (uint64_t{1} << ib)) == 0U) continue;
        const unsigned int diff = ia > ib ? ia - ib : ib - ia;
        best = std::min(best, diff > 0U ? diff - 1U : 0U);
      }
    }
    return best;
  }

  // NOLINTBEGIN(misc-non-private-member-variables-in-classes)
  struct PieceStartSort
  {
    const std::vector<Tpc_AssembledTrackReco::Piece>* pieces;
    explicit PieceStartSort(const std::vector<Tpc_AssembledTrackReco::Piece>* p)
      : pieces(p)
    {
    }
    bool operator()(unsigned int a, unsigned int b) const
    {
      const Tpc_AssembledTrackReco::Piece& pa = (*pieces)[a];
      const Tpc_AssembledTrackReco::Piece& pb = (*pieces)[b];
      if (pa.nblobs != pb.nblobs)
      {
        return pa.nblobs > pb.nblobs;
      }
      const unsigned int span_a = pa.last_layer >= pa.first_layer
                                      ? pa.last_layer - pa.first_layer + 1U
                                      : 0U;
      const unsigned int span_b = pb.last_layer >= pb.first_layer
                                      ? pb.last_layer - pb.first_layer + 1U
                                      : 0U;
      if (span_a != span_b)
      {
        return span_a > span_b;
      }
      if (pa.first_layer != pb.first_layer)
      {
        return pa.first_layer < pb.first_layer;
      }
      if (pa.sector != pb.sector)
      {
        return pa.sector < pb.sector;
      }
      return pa.source_track_id < pb.source_track_id;
    }
  };

  struct CandidateStartSort
  {
    const std::vector<Tpc_AssembledTrackReco::Candidate>* candidates;
    explicit CandidateStartSort(const std::vector<Tpc_AssembledTrackReco::Candidate>* c)
      : candidates(c)
    {
    }
    bool operator()(unsigned int a, unsigned int b) const
    {
      const Tpc_AssembledTrackReco::Candidate& ca = (*candidates)[a];
      const Tpc_AssembledTrackReco::Candidate& cb = (*candidates)[b];
      if (ca.nblobs != cb.nblobs)
      {
        return ca.nblobs > cb.nblobs;
      }
      const unsigned int span_a = ca.last_layer >= ca.first_layer
                                      ? ca.last_layer - ca.first_layer + 1U
                                      : 0U;
      const unsigned int span_b = cb.last_layer >= cb.first_layer
                                      ? cb.last_layer - cb.first_layer + 1U
                                      : 0U;
      if (span_a != span_b)
      {
        return span_a > span_b;
      }
      if (ca.first_layer != cb.first_layer)
      {
        return ca.first_layer < cb.first_layer;
      }
      if (ca.first_sector != cb.first_sector)
      {
        return ca.first_sector < cb.first_sector;
      }
      return ca.nsegments > cb.nsegments;
    }
  };

  // NOLINTEND(misc-non-private-member-variables-in-classes)

  double unwrap_phi_to_reference(double phi, const double ref)
  {
    while (phi - ref > M_PI)
    {
      phi -= 2.0 * M_PI;
    }
    while (phi - ref < -M_PI)
    {
      phi += 2.0 * M_PI;
    }
    return phi;
  }

  double wrap_to_pi(double phi)
  {
    while (phi > M_PI)
    {
      phi -= 2.0 * M_PI;
    }
    while (phi <= -M_PI)
    {
      phi += 2.0 * M_PI;
    }
    return phi;
  }

  int wrapped_sector_delta(const unsigned int sector_a, const unsigned int sector_b)
  {
    int d = static_cast<int>(sector_b) - static_cast<int>(sector_a);
    while (d > 6)
    {
      d -= 12;
    }
    while (d < -6)
    {
      d += 12;
    }
    return d;
  }

  double sagitta_model_derivative(double xrot, double x0, double invR)
  {
    const double dx = xrot - x0;
    const double dx2 = dx * dx;
    const double invR2 = invR * invR;
    const double invR3 = invR2 * invR;
    const double invR5 = invR3 * invR2;
    return -invR * dx - 0.5 * invR3 * dx2 * dx - 0.375 * invR5 * dx2 * dx2 * dx;
  }

  std::pair<double, double> predict_sagitta_phi_and_slope(double radius,
                                                          double S,
                                                          double x0,
                                                          double invR,
                                                          double theta,
                                                          double bline)
  {
    const double c = std::cos(theta);
    const double s = std::sin(theta);
    double yy = std::tan(theta) * radius;

    for (unsigned int iter = 0; iter < 25; ++iter)
    {
      const double xrot = c * radius + s * yy;
      const double yrot = -s * radius + c * yy;
      const double f = Tpc_FittingTools::sagittaModel(xrot, S, x0, invR);
      const double g = yrot - f;
      const double df = sagitta_model_derivative(xrot, x0, invR);
      const double dg = c - df * s;
      if (std::fabs(dg) < 1.0e-12)
      {
        break;
      }
      const double step = g / dg;
      yy -= step;
      if (std::fabs(step) < 1.0e-10)
      {
        break;
      }
    }

    const double xrot = c * radius + s * yy;
    const double df = sagitta_model_derivative(xrot, x0, invR);
    const double denom = c - df * s;
    const double slope = std::fabs(denom) > 1.0e-12
                             ? (s + df * c) / denom
                             : std::tan(theta);
    return {bline + yy, slope};
  }

  bool fit_points(const std::vector<Tpc_AssembledTrackReco::PiecePoint>& points,
                  bool use_sagitta,
                  double& phi_slope,
                  double& phi_intercept,
                  double& phi_S,
                  double& phi_x0,
                  double& phi_invR,
                  double& phi_theta,
                  double& phi_bline,
                  bool& phi_sagitta_ok,
                  double& tbin_slope,
                  double& tbin_intercept,
                  double& chi2_phi,
                  double& chi2_tbin,
                  int& ndof_phi,
                  int& ndof_tbin)
  {
    if (points.size() < 2)
    {
      return false;
    }

    // Tpc_FittingTools currently consumes FitPoint vectors.
    // Build only the two arrays actually required by the fit instead
    // of four parallel radius/phi/tbin/weight arrays.
    std::vector<Tpc_FittingTools::FitPoint> phi_points;
    std::vector<Tpc_FittingTools::FitPoint> tbin_points;
    phi_points.reserve(points.size());
    tbin_points.reserve(points.size());

    for (const auto& point : points)
    {
      phi_points.emplace_back(point.radius, point.phi, point.weight);
      tbin_points.emplace_back(point.radius, point.tbin, point.weight);
    }

    const Tpc_FittingTools::LineFit phi_line = Tpc_FittingTools::fitLine(phi_points);
    if (!phi_line.ok)
    {
      return false;
    }

    phi_slope = phi_line.slope;
    phi_intercept = phi_line.intercept;
    chi2_phi = phi_line.chi2;
    ndof_phi = phi_line.ndof;

    phi_S = 0.0;
    phi_x0 = 0.0;
    phi_invR = 0.0;
    phi_theta = std::atan(phi_slope);
    phi_bline = phi_intercept;
    phi_sagitta_ok = false;

    if (use_sagitta && points.size() >= 3)
    {
      const Tpc_FittingTools::SagittaFit phi_sagitta =
          Tpc_FittingTools::fitSagitta(phi_points);

      if (phi_sagitta.ok)
      {
        phi_S = phi_sagitta.S;
        phi_x0 = phi_sagitta.x0;
        phi_invR = phi_sagitta.invR;
        phi_theta = phi_sagitta.theta;
        phi_bline = phi_sagitta.b;
        chi2_phi = phi_sagitta.chi2;
        ndof_phi = phi_sagitta.ndof;
        phi_sagitta_ok = true;
      }
    }

    const Tpc_FittingTools::LineFit tbin_line =
        Tpc_FittingTools::fitLine(tbin_points);

    if (!tbin_line.ok)
    {
      return false;
    }

    tbin_slope = tbin_line.slope;
    tbin_intercept = tbin_line.intercept;
    chi2_tbin = tbin_line.chi2;
    ndof_tbin = tbin_line.ndof;

    return true;
  }

}  // namespace

Tpc_AssembledTrackReco::Piece::Piece()
  : source_index(0)
  , source_track_id(0)
  , event(0)
  , region(0)
  , sector(0)
  , side(0)
  , first_layer(0)
  , last_layer(0)
  , nblobs(0)
  , nrawhits(0)
  , layer_mask(0)
  , phi_slope(0.0)
  , phi_intercept(0.0)
  , phi_S(0.0)
  , phi_x0(0.0)
  , phi_invR(0.0)
  , phi_theta(0.0)
  , phi_bline(0.0)
  , phi_sagitta_ok(false)
  , phi_sagitta_evaluated(false)
  , tbin_slope(0.0)
  , tbin_intercept(0.0)
{
}

Tpc_AssembledTrackReco::Candidate::Candidate()
  : event(0)
  , side(0)
  , first_layer(0)
  , last_layer(0)
  , first_sector(0)
  , last_sector(0)
  , first_region(0)
  , last_region(0)
  , nsegments(0)
  , nblobs(0)
  , nrawhits(0)
  , layer_mask(0)
  , phi_slope(0.0)
  , phi_intercept(0.0)
  , phi_S(0.0)
  , phi_x0(0.0)
  , phi_invR(0.0)
  , phi_theta(0.0)
  , phi_bline(0.0)
  , phi_sagitta_ok(false)
  , tbin_slope_r(0.0)
  , tbin_intercept_r(0.0)
  , chi2_phi(0.0)
  , chi2_tbin(0.0)
  , ndof_phi(0)
  , ndof_tbin(0)
{
}

Tpc_AssembledTrackReco::Tpc_AssembledTrackReco(const std::string& name, const std::string& filename)
  : SubsysReco(name)
  , m_outputFileName(filename)
  , m_debugOutputFileName("Tpc_AssembledTrackRecoDebug.root")
  , m_inputNodeName("TPC_MODULETRACKS")
  , m_outputNodeName("TPC_ASSEMBLEDTRACKS")
  , m_outputFile(nullptr)
  , m_debugOutputFile(nullptr)
  , m_tree(nullptr)
  , m_tpcModuleTrackContainer(nullptr)
  , m_assembledTrackContainer(nullptr)
  , m_hits(nullptr)
  , m_event(0)
  , m_idealPadMap(nullptr)
  , m_connectMaxLayerGap{{16, 16, 16}}
  , m_connect_dphi{{0.040, 0.040, 0.035}}
  , m_connect_dtbin{{8.0, 8.0, 8.0}}
  , m_connect_dphi_slope{{0.018, 0.015, 0.012}}
  , m_connect_dtbin_slope{{2.0, 2.0, 2.0}}
  , m_useSagittaPhiFit(true)
  , m_useAnalyticSagittaSlope(true)
  , m_doDebugHistograms(false)
  , m_minClustersPerTrack(0)
  , m_linearPhiPrecutScale(8.0)
  , m_minSagittaLayers(8)
  , m_robustPieceMinBlobs(8)
  , m_score_dphi{{0.018, 0.010, 0.012}}
  , m_score_dtbin{{2.0, 1.5, 1.6}}
  , m_score_dphi_slope{{0.0070, 0.0045, 0.0045}}
  , m_score_dtbin_slope{{0.90, 0.50, 0.50}}
  , m_weight_phi{{0.75, 0.75, 0.75}}
  , m_weight_tbin{{1.50, 1.50, 1.50}}
  , m_weight_phi_slope{{0.75, 0.75, 0.75}}
  , m_weight_tbin_slope{{2.00, 2.00, 2.00}}
  , m_enableFastSameSectorPass(true)
  , m_fastBroadScale(1.0)
  , m_fastTightScale(0.35)
  , m_enableNormalUniqueTight(true)
  , m_normalTightScale(0.5)
  , m_enableRescuePass(true)
  , m_rescueIterations(2)
  , m_rescueMaxWindowScale(2.0)
  , m_rescueMaxSegments(2)
  , m_enableFinalReassociation(false)
  , m_finalReassociationIterations(2)
  , m_finalReassociationImprovementRatio(0.8)
  , m_seedSigmaX(5.0)
  , m_seedSigmaY(5.0)
  , m_seedSigmaZ(10.0)
  , m_seedSigmaPx(1.0)
  , m_seedSigmaPy(1.0)
  , m_seedSigmaPz(1.0)
  , m_h_track_tbin_slope_vs_tbin_span_3modules(nullptr)
  , m_h_track_tbin_slope_vs_first_tbin_3modules(nullptr)
  , m_h_track_tbin_slope_vs_last_tbin_3modules(nullptr)
  , m_h_nsegments(nullptr)
{
}

Tpc_AssembledTrackReco::~Tpc_AssembledTrackReco()
{
  delete m_idealPadMap;
  m_idealPadMap = nullptr;

  if (m_debugOutputFile)
  {
    delete m_debugOutputFile;
    m_debugOutputFile = nullptr;
  }

  if (m_outputFile)
  {
    delete m_outputFile;
    m_outputFile = nullptr;
  }
}

int Tpc_AssembledTrackReco::Init(PHCompositeNode* /*unused*/)
{
  if (Verbosity() <= 0)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  m_outputFile = new TFile(m_outputFileName.c_str(), "RECREATE");
  if (!m_outputFile || m_outputFile->IsZombie())
  {
    std::cerr << Name() << "::Init - cannot create " << m_outputFileName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  if (m_doDebugHistograms)
  {
    m_debugOutputFile = new TFile(m_debugOutputFileName.c_str(), "RECREATE");
    if (!m_debugOutputFile || m_debugOutputFile->IsZombie())
    {
      std::cerr << Name() << "::Init - cannot create debug file "
                << m_debugOutputFileName << std::endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
    create_debug_histograms();
  }

  m_tree = new TTree("Tpc_AssembledTracks", "Assembled tracks connected from Tpc_ModuleTrackReco");
  m_tree->Branch("event", &m_tree_event, "event/I");
  m_tree->Branch("track_id", &m_tree_track_id);
  m_tree->Branch("side", &m_tree_side);
  m_tree->Branch("nsegments", &m_tree_nsegments);
  m_tree->Branch("nblobs", &m_tree_nblobs);
  m_tree->Branch("nrawhits", &m_tree_nrawhits);
  m_tree->Branch("first_layer", &m_tree_first_layer);
  m_tree->Branch("last_layer", &m_tree_last_layer);
  m_tree->Branch("first_sector", &m_tree_first_sector);
  m_tree->Branch("last_sector", &m_tree_last_sector);
  m_tree->Branch("first_region", &m_tree_first_region);
  m_tree->Branch("last_region", &m_tree_last_region);

  m_tree->Branch("source_assembled_track_id", &m_tree_source_assembled_track_id);
  m_tree->Branch("source_inmodule_track_id", &m_tree_source_inmodule_track_id);
  m_tree->Branch("source_region", &m_tree_source_region);
  m_tree->Branch("source_sector", &m_tree_source_sector);
  m_tree->Branch("source_side", &m_tree_source_side);

  m_tree->Branch("hit_assembled_track_id", &m_tree_hit_assembled_track_id);
  m_tree->Branch("hit_hitsetkey", &m_tree_hit_hitsetkey);
  m_tree->Branch("hit_hitkey", &m_tree_hit_hitkey);

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_AssembledTrackReco::InitRun(PHCompositeNode* topNode)
{
  if (getNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }
  if (createNodes(topNode) != Fun4AllReturnCodes::EVENT_OK)
  {
    return Fun4AllReturnCodes::ABORTRUN;
  }

  delete m_idealPadMap;
  m_idealPadMap = new IdealPadMap();
  if (m_idealPadMap->load_from_cdb(Verbosity()) != 0 || !m_idealPadMap->is_loaded())
  {
    std::cerr << Name() << "::InitRun - cannot load IdealPadMap from CDB" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_event = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_AssembledTrackReco::End(PHCompositeNode* /*unused*/)
{
  if (m_outputFile)
  {
    m_outputFile->cd();
    if (m_tree)
    {
      m_tree->Write();
    }
    m_outputFile->Close();
    delete m_outputFile;
    m_outputFile = nullptr;
  }

  write_debug_histograms();
  if (m_debugOutputFile)
  {
    m_debugOutputFile->Close();
    delete m_debugOutputFile;
    m_debugOutputFile = nullptr;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_AssembledTrackReco::getNodes(PHCompositeNode* topNode)
{
  m_tpcModuleTrackContainer = findNode::getClass<Tpc_ModuleTrackContainer>(topNode, m_inputNodeName);
  if (!m_tpcModuleTrackContainer)
  {
    std::cerr << Name() << "::getNodes - missing " << m_inputNodeName << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  m_hits = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!m_hits)
  {
    std::cerr << Name() << "::getNodes - missing TRKR_HITSET" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

int Tpc_AssembledTrackReco::createNodes(PHCompositeNode* topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode* dstNode = dynamic_cast<PHCompositeNode*>(iter.findFirst("PHCompositeNode", "DST"));

  if (!dstNode)
  {
    dstNode = new PHCompositeNode("DST");
    topNode->addNode(dstNode);
  }

  m_assembledTrackContainer = findNode::getClass<Tpc_AssembledTrackContainer>(topNode, m_outputNodeName);
  if (!m_assembledTrackContainer)
  {
    m_assembledTrackContainer = new Tpc_AssembledTrackContainerv1();
    PHIODataNode<PHObject>* node = new PHIODataNode<PHObject>(m_assembledTrackContainer, m_outputNodeName, "PHObject");
    dstNode->addNode(node);
    std::cout << Name() << "::createNodes - created " << m_outputNodeName << " node" << std::endl;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

void Tpc_AssembledTrackReco::create_debug_histograms()
{
  if (!m_debugOutputFile)
  {
    return;
  }
  m_debugOutputFile->cd();

  for (unsigned int region = 0; region < kRegionCount; ++region)
  {
    const std::string rtag = "r" + std::to_string(region + 1);
    const std::string prefix = "R" + std::to_string(region + 1) + " ";
    const double dphi_cut = m_connect_dphi[region];
    const double dtbin_cut = m_connect_dtbin[region];
    const double dmphi_cut = m_connect_dphi_slope[region];
    const double dmtbin_cut = m_connect_dtbin_slope[region];

    m_h_dphi[region] = new TH1D(("h_dphi_" + rtag).c_str(), std::format("{}#Delta#phi at match point, cut={:.4g};|#Delta#phi| [rad];tested pairs", prefix, dphi_cut).c_str(), 200, 0.0, std::max(0.2, 5.0 * dphi_cut));
    m_h_dtbin[region] = new TH1D(("h_dtbin_" + rtag).c_str(), std::format("{}#Deltatbin at match point, cut={:.4g};|#Deltatbin|;tested pairs", prefix, dtbin_cut).c_str(), 200, 0.0, std::max(50.0, 5.0 * dtbin_cut));
    m_h_dmphi[region] = new TH1D(("h_dmphi_" + rtag).c_str(), std::format("{}#Delta(d#phi/dr), cut={:.4g};|#Delta(d#phi/dr)| [rad/cm];tested pairs", prefix, dmphi_cut).c_str(), 200, 0.0, std::max(0.08, 5.0 * dmphi_cut));
    m_h_dmtbin[region] = new TH1D(("h_dmtbin_" + rtag).c_str(), std::format("{}#Delta(dtbin/dr), cut={:.4g};|#Delta(dtbin/dr)| [tbin/cm];tested pairs", prefix, dmtbin_cut).c_str(), 200, 0.0, std::max(20.0, 5.0 * dmtbin_cut));
    m_h_score[region] = new TH1D(("h_score_" + rtag).c_str(), (prefix + "accepted connection score;score;accepted connections").c_str(), 200, 0.0, 20.0);

    m_h_dphi_vs_dtbin[region] = new TH2D(("h_dphi_vs_dtbin_" + rtag).c_str(), std::format("{}#Delta#phi vs #Deltatbin, cuts #Delta#phi<{:.4g} #Deltatbin<{:.4g};|#Delta#phi| [rad];|#Deltatbin|", prefix, dphi_cut, dtbin_cut).c_str(), 160, 0.0, std::max(0.2, 5.0 * dphi_cut), 160, 0.0, std::max(50.0, 5.0 * dtbin_cut));
    m_h_dmphi_vs_dmtbin[region] = new TH2D(("h_dmphi_vs_dmtbin_" + rtag).c_str(), std::format("{}slope residuals, cuts #Delta(d#phi/dr)<{:.4g} #Delta(dtbin/dr)<{:.4g};|#Delta(d#phi/dr)| [rad/cm];|#Delta(dtbin/dr)| [tbin/cm]", prefix, dmphi_cut, dmtbin_cut).c_str(), 160, 0.0, std::max(0.08, 5.0 * dmphi_cut), 160, 0.0, std::max(20.0, 5.0 * dmtbin_cut));
    m_h_dphi_vs_dmphi[region] = new TH2D(("h_dphi_vs_dmphi_" + rtag).c_str(), std::format("{}#phi position vs slope residual, cuts #Delta#phi<{:.4g} #Delta(d#phi/dr)<{:.4g};|#Delta#phi| [rad];|#Delta(d#phi/dr)| [rad/cm]", prefix, dphi_cut, dmphi_cut).c_str(), 160, 0.0, std::max(0.2, 5.0 * dphi_cut), 160, 0.0, std::max(0.08, 5.0 * dmphi_cut));

    m_h_tbin_slope_vs_first_tbin[region] = new TH2D(("h_tbin_slope_vs_first_tbin_" + rtag).c_str(), (prefix + "tested connection dtbin/dr vs first timebin;first timebin;dtbin/dr [tbin/cm]").c_str(), 200, 0.0, 600.0, 200, -20.0, 20.0);
    m_h_tbin_slope_vs_last_tbin[region] = new TH2D(("h_tbin_slope_vs_last_tbin_" + rtag).c_str(), (prefix + "tested connection dtbin/dr vs last timebin;last timebin;dtbin/dr [tbin/cm]").c_str(), 200, 0.0, 600.0, 200, -20.0, 20.0);

    m_h_layer_gap[region] = new TH1D(("h_layer_gap_" + rtag).c_str(), (prefix + "accepted connection layer gap;b.first_layer - a.last_layer - 1;accepted connections").c_str(), 16, -0.5, 15.5);
    m_h_matched_sector_delta[region] = new TH1D(("h_matched_sector_delta_" + rtag).c_str(), (prefix + "accepted matched sector difference;wrapped #Delta sector;accepted connections").c_str(), 25, -12.5, 12.5);
  }

  m_h_track_tbin_slope_vs_tbin_span_3modules = new TH2D("h_track_tbin_slope_vs_tbin_span_3modules", "3-module tracks dtbin/dr vs last-first timebin;last timebin - first timebin;dtbin/dr [tbin/cm]", 200, -600.0, 600.0, 200, -20.0, 20.0);
  m_h_track_tbin_slope_vs_first_tbin_3modules = new TH2D("h_track_tbin_slope_vs_first_tbin_3modules", "3-module tracks dtbin/dr vs first timebin;first timebin;dtbin/dr [tbin/cm]", 200, 0.0, 600.0, 200, -20.0, 20.0);
  m_h_track_tbin_slope_vs_last_tbin_3modules = new TH2D("h_track_tbin_slope_vs_last_tbin_3modules", "3-module tracks dtbin/dr vs last timebin;last timebin;dtbin/dr [tbin/cm]", 200, 0.0, 600.0, 200, -20.0, 20.0);
  m_h_nsegments = new TH1D("h_nsegments", "pieces per assembled track;nsegments;assembled tracks", 16, -0.5, 15.5);
}

void Tpc_AssembledTrackReco::write_debug_histograms()
{
  if (!m_debugOutputFile)
  {
    return;
  }
  m_debugOutputFile->cd();

  for (unsigned int region = 0; region < kRegionCount; ++region)
  {
    if (m_h_dphi[region]) m_h_dphi[region]->Write();
    if (m_h_dtbin[region]) m_h_dtbin[region]->Write();
    if (m_h_dmphi[region]) m_h_dmphi[region]->Write();
    if (m_h_dmtbin[region]) m_h_dmtbin[region]->Write();
    if (m_h_score[region]) m_h_score[region]->Write();
    if (m_h_dphi_vs_dtbin[region]) m_h_dphi_vs_dtbin[region]->Write();
    if (m_h_dmphi_vs_dmtbin[region]) m_h_dmphi_vs_dmtbin[region]->Write();
    if (m_h_dphi_vs_dmphi[region]) m_h_dphi_vs_dmphi[region]->Write();
    if (m_h_tbin_slope_vs_first_tbin[region]) m_h_tbin_slope_vs_first_tbin[region]->Write();
    if (m_h_tbin_slope_vs_last_tbin[region]) m_h_tbin_slope_vs_last_tbin[region]->Write();
    if (m_h_layer_gap[region]) m_h_layer_gap[region]->Write();
    if (m_h_matched_sector_delta[region]) m_h_matched_sector_delta[region]->Write();
  }

  if (m_h_track_tbin_slope_vs_tbin_span_3modules) m_h_track_tbin_slope_vs_tbin_span_3modules->Write();
  if (m_h_track_tbin_slope_vs_first_tbin_3modules) m_h_track_tbin_slope_vs_first_tbin_3modules->Write();
  if (m_h_track_tbin_slope_vs_last_tbin_3modules) m_h_track_tbin_slope_vs_last_tbin_3modules->Write();
  if (m_h_nsegments) m_h_nsegments->Write();
}

void Tpc_AssembledTrackReco::reset_tree_vars()
{
  m_tree_event = m_event;
  m_tree_track_id.clear();
  m_tree_side.clear();
  m_tree_nsegments.clear();
  m_tree_nblobs.clear();
  m_tree_nrawhits.clear();
  m_tree_first_layer.clear();
  m_tree_last_layer.clear();
  m_tree_first_sector.clear();
  m_tree_last_sector.clear();
  m_tree_first_region.clear();
  m_tree_last_region.clear();
  m_tree_source_assembled_track_id.clear();
  m_tree_source_inmodule_track_id.clear();
  m_tree_source_region.clear();
  m_tree_source_sector.clear();
  m_tree_source_side.clear();
  m_tree_hit_assembled_track_id.clear();
  m_tree_hit_hitsetkey.clear();
  m_tree_hit_hitkey.clear();
}

bool Tpc_AssembledTrackReco::make_piece(unsigned int source_index, Piece& p) const
{
  const Tpc_ModuleTrack* trk = m_tpcModuleTrackContainer->get_track(source_index);
  if (!trk || !trk->isValid())
  {
    return false;
  }
  if (trk->get_last_layer() < trk->get_first_layer())
  {
    return false;
  }

  if (!m_idealPadMap || !m_idealPadMap->is_loaded() || !m_hits)
  {
    return false;
  }

  p.source_index = source_index;
  p.source_track_id = trk->get_track_id();
  p.event = trk->get_event();
  p.region = trk->get_region();
  p.sector = trk->get_sector();
  p.side = trk->get_side();
  p.first_layer = trk->get_first_layer();
  p.last_layer = trk->get_last_layer();
  p.nblobs = trk->get_nblobs();
  p.nrawhits = trk->get_nrawhits();
  p.layer_mask = 0;

  p.points.clear();
  p.hit_indices.clear();
  p.hit_indices.reserve(trk->size_hit_indices());

  struct LayerAccumulator
  {
    bool used{false};
    double radius{0.0};
    double adc_sum{0.0};
    double tbin_adc_sum{0.0};
    double cos_adc_sum{0.0};
    double sin_adc_sum{0.0};
    TrkrDefs::hitsetkey first_hitsetkey{0};
    TrkrDefs::hitkey first_hitkey{0};
  };

  std::array<LayerAccumulator, 48> layer_acc{};

  for (unsigned int ih = 0; ih < trk->size_hit_indices(); ++ih)
  {
    const Tpc_ModuleTrack::HitIndex hi = trk->get_hit_index(ih);
    TrkrHitSet* hitset = m_hits->findHitSet(hi.first);
    if (!hitset)
    {
      continue;
    }

    TrkrHit* hit = hitset->getHit(hi.second);
    if (!hit)
    {
      continue;
    }

    const unsigned int layer = TrkrDefs::getLayer(hi.first);
    if (layer < 7U || layer > 54U)
    {
      continue;
    }

    const unsigned int pad = TpcDefs::getPad(hi.second);
    const unsigned int tbin = TpcDefs::getTBin(hi.second);
    const double radius = m_idealPadMap->get_radius(layer);
    const double phi =
        wrap_to_pi(m_idealPadMap->get_phi(static_cast<unsigned int>(p.side),
                                         layer,
                                         pad));
    if (!std::isfinite(radius) || !std::isfinite(phi))
    {
      continue;
    }

    p.hit_indices.emplace_back(hi.first, hi.second);
    p.layer_mask |= (uint64_t{1} << (layer - 7U));

    auto& acc = layer_acc[layer - 7U];
    const double adc = std::max(1.0, static_cast<double>(hit->getAdc()));
    if (!acc.used)
    {
      acc.used = true;
      acc.radius = radius;
      acc.first_hitsetkey = hi.first;
      acc.first_hitkey = hi.second;
    }
    acc.adc_sum += adc;
    acc.tbin_adc_sum += adc * static_cast<double>(tbin);
    acc.cos_adc_sum += adc * std::cos(phi);
    acc.sin_adc_sum += adc * std::sin(phi);
  }

  double max_layer_adc = 0.0;
  for (const auto& acc : layer_acc)
  {
    if (acc.used)
    {
      max_layer_adc = std::max(max_layer_adc, acc.adc_sum);
    }
  }

  p.points.reserve(p.nblobs);
  for (unsigned int ilayer = 0; ilayer < layer_acc.size(); ++ilayer)
  {
    const auto& acc = layer_acc[ilayer];
    if (!acc.used || acc.adc_sum <= 0.0)
    {
      continue;
    }

    PiecePoint point;
    point.radius = acc.radius;
    point.phi = std::atan2(acc.sin_adc_sum, acc.cos_adc_sum);
    point.tbin = acc.tbin_adc_sum / acc.adc_sum;
    point.weight = Tpc_FittingTools::adcWeight(acc.adc_sum,
                                                max_layer_adc,
                                                0.5,
                                                0.15);
    point.hitsetkey = acc.first_hitsetkey;
    point.hitkey = acc.first_hitkey;
    p.points.push_back(point);
  }

  if (p.points.size() < 2)
  {
    return false;
  }

  // The assembler fit now has one statistically independent centroid per
  // layer/blob rather than one point per raw hit.  This keeps a broad blob
  // from receiving tens of times more fit weight than a narrow blob.
  std::sort(p.points.begin(), p.points.end(),
            [](const PiecePoint& a, const PiecePoint& b)
            {
              return a.radius < b.radius;
            });

  for (unsigned int i = 1; i < p.points.size(); ++i)
  {
    p.points[i].phi =
        unwrap_phi_to_reference(p.points[i].phi, p.points[i - 1].phi);
  }

  double chi2_phi = 0.0;
  double chi2_tbin = 0.0;
  int ndof_phi = 0;
  int ndof_tbin = 0;

  const bool fit_ok =
      fit_points(p.points,
                 false,
                 p.phi_slope,
                 p.phi_intercept,
                 p.phi_S,
                 p.phi_x0,
                 p.phi_invR,
                 p.phi_theta,
                 p.phi_bline,
                 p.phi_sagitta_ok,
                 p.tbin_slope,
                 p.tbin_intercept,
                 chi2_phi,
                 chi2_tbin,
                 ndof_phi,
                 ndof_tbin);

  p.phi_sagitta_evaluated = !m_useSagittaPhiFit ||
                            p.points.size() < m_minSagittaLayers;
  return fit_ok;
}

void Tpc_AssembledTrackReco::ensure_piece_sagitta(const Piece& p) const
{
  if (!m_useSagittaPhiFit || p.phi_sagitta_evaluated)
  {
    return;
  }

  // Mark evaluated even if the fit fails so a bad Piece is not refitted
  // repeatedly every time it reaches exact matching.
  p.phi_sagitta_evaluated = true;
  p.phi_sagitta_ok = false;

  if (p.points.size() < m_minSagittaLayers)
  {
    return;
  }

  std::vector<Tpc_FittingTools::FitPoint> phi_points;
  phi_points.reserve(p.points.size());

  for (const PiecePoint& point : p.points)
  {
    phi_points.emplace_back(point.radius, point.phi, point.weight);
  }

  const Tpc_FittingTools::SagittaFit phi_sagitta =
      Tpc_FittingTools::fitSagitta(phi_points);

  if (!phi_sagitta.ok)
  {
    return;
  }

  p.phi_S = phi_sagitta.S;
  p.phi_x0 = phi_sagitta.x0;
  p.phi_invR = phi_sagitta.invR;
  p.phi_theta = phi_sagitta.theta;
  p.phi_bline = phi_sagitta.b;
  p.phi_sagitta_ok = true;
}

double Tpc_AssembledTrackReco::predict_phi(const Piece& p, double radius) const
{
  if (m_useSagittaPhiFit)
  {
    ensure_piece_sagitta(p);
  }
  if (m_useSagittaPhiFit && p.phi_sagitta_ok)
  {
    return predict_sagitta_phi_and_slope(radius,
                                         p.phi_S,
                                         p.phi_x0,
                                         p.phi_invR,
                                         p.phi_theta,
                                         p.phi_bline).first;
  }
  return p.phi_slope * radius + p.phi_intercept;
}

double Tpc_AssembledTrackReco::predict_phi(const Candidate& c, double radius) const
{
  if (m_useSagittaPhiFit && c.phi_sagitta_ok)
  {
    return predict_sagitta_phi_and_slope(radius,
                                         c.phi_S,
                                         c.phi_x0,
                                         c.phi_invR,
                                         c.phi_theta,
                                         c.phi_bline).first;
  }
  return c.phi_slope * radius + c.phi_intercept;
}

double Tpc_AssembledTrackReco::predict_phi_slope(const Piece& p, double radius) const
{
  if (m_useSagittaPhiFit)
  {
    ensure_piece_sagitta(p);
  }
  if (!(m_useSagittaPhiFit && p.phi_sagitta_ok))
  {
    return p.phi_slope;
  }
  if (m_useAnalyticSagittaSlope)
  {
    return predict_sagitta_phi_and_slope(radius,
                                         p.phi_S,
                                         p.phi_x0,
                                         p.phi_invR,
                                         p.phi_theta,
                                         p.phi_bline).second;
  }
  const double eps = 1.0e-3;
  return (predict_phi(p, radius + eps) - predict_phi(p, radius - eps)) / (2.0 * eps);
}

double Tpc_AssembledTrackReco::predict_phi_slope(const Candidate& c, double radius) const
{
  if (!(m_useSagittaPhiFit && c.phi_sagitta_ok))
  {
    return c.phi_slope;
  }
  if (m_useAnalyticSagittaSlope)
  {
    return predict_sagitta_phi_and_slope(radius,
                                         c.phi_S,
                                         c.phi_x0,
                                         c.phi_invR,
                                         c.phi_theta,
                                         c.phi_bline).second;
  }
  const double eps = 1.0e-3;
  return (predict_phi(c, radius + eps) - predict_phi(c, radius - eps)) / (2.0 * eps);
}

std::pair<double, double> Tpc_AssembledTrackReco::predict_phi_and_slope(const Piece& p, double radius) const
{
  if (m_useSagittaPhiFit)
  {
    ensure_piece_sagitta(p);
  }
  if (m_useSagittaPhiFit && p.phi_sagitta_ok && m_useAnalyticSagittaSlope)
  {
    return predict_sagitta_phi_and_slope(radius,
                                         p.phi_S,
                                         p.phi_x0,
                                         p.phi_invR,
                                         p.phi_theta,
                                         p.phi_bline);
  }
  return {predict_phi(p, radius), predict_phi_slope(p, radius)};
}

std::pair<double, double> Tpc_AssembledTrackReco::predict_phi_and_slope(const Candidate& c, double radius) const
{
  if (m_useSagittaPhiFit && c.phi_sagitta_ok && m_useAnalyticSagittaSlope)
  {
    return predict_sagitta_phi_and_slope(radius,
                                         c.phi_S,
                                         c.phi_x0,
                                         c.phi_invR,
                                         c.phi_theta,
                                         c.phi_bline);
  }
  return {predict_phi(c, radius), predict_phi_slope(c, radius)};
}

bool Tpc_AssembledTrackReco::refit_candidate(const std::vector<Piece>& pieces,
                                               const std::vector<unsigned int>& piece_indices,
                                               Candidate& c) const
{
  if (piece_indices.empty())
  {
    return false;
  }

  std::vector<unsigned int> ordered_indices = piece_indices;
  std::sort(ordered_indices.begin(), ordered_indices.end(),
            [&pieces](unsigned int a, unsigned int b)
            {
              const Piece& pa = pieces[a];
              const Piece& pb = pieces[b];
              if (pa.first_layer != pb.first_layer)
              {
                return pa.first_layer < pb.first_layer;
              }
              if (pa.last_layer != pb.last_layer)
              {
                return pa.last_layer < pb.last_layer;
              }
              if (pa.sector != pb.sector)
              {
                return pa.sector < pb.sector;
              }
              return pa.source_track_id < pb.source_track_id;
            });

  c = Candidate();
  c.piece_indices = ordered_indices;

  std::size_t total_points = 0;
  for (unsigned int piece_index : ordered_indices)
  {
    total_points += pieces[piece_index].points.size();
  }

  std::vector<PiecePoint> points;
  points.reserve(total_points);

  for (unsigned int ii = 0; ii < ordered_indices.size(); ++ii)
  {
    const Piece& p = pieces[ordered_indices[ii]];
    c.layer_mask |= p.layer_mask;

    for (const PiecePoint& source_point : p.points)
    {
      points.push_back(source_point);
    }

    if (ii == 0)
    {
      c.event = p.event;
      c.side = p.side;
      c.first_layer = p.first_layer;
      c.last_layer = p.last_layer;
      c.first_sector = p.sector;
      c.last_sector = p.sector;
      c.first_region = p.region;
      c.last_region = p.region;
    }
    else
    {
      if (p.first_layer < c.first_layer)
      {
        c.first_layer = p.first_layer;
        c.first_sector = p.sector;
        c.first_region = p.region;
      }

      if (p.last_layer > c.last_layer)
      {
        c.last_layer = p.last_layer;
        c.last_sector = p.sector;
        c.last_region = p.region;
      }
    }

    c.nblobs += p.nblobs;
    c.nrawhits += p.nrawhits;
  }

  if (points.size() < 2)
  {
    return false;
  }

  // Always fit in radial order, independent of the order in which pieces were
  // attached.  This is essential once a long middle piece can collect a short
  // inner piece after it has already collected an outer continuation.
  std::sort(points.begin(), points.end(),
            [](const PiecePoint& a, const PiecePoint& b)
            {
              return a.radius < b.radius;
            });
  for (unsigned int i = 1; i < points.size(); ++i)
  {
    points[i].phi = unwrap_phi_to_reference(points[i].phi, points[i - 1].phi);
  }

  c.nsegments = static_cast<unsigned int>(ordered_indices.size());

  return fit_points(points,
                    m_useSagittaPhiFit && points.size() >= m_minSagittaLayers,
                    c.phi_slope,
                    c.phi_intercept,
                    c.phi_S,
                    c.phi_x0,
                    c.phi_invR,
                    c.phi_theta,
                    c.phi_bline,
                    c.phi_sagitta_ok,
                    c.tbin_slope_r,
                    c.tbin_intercept_r,
                    c.chi2_phi,
                    c.chi2_tbin,
                    c.ndof_phi,
                    c.ndof_tbin);
}

bool Tpc_AssembledTrackReco::candidates_can_connect(const Candidate& a,
                                                      const Piece& b,
                                                      double& score,
                                                      double& b_phi_intercept_shifted) const
{
  score = std::numeric_limits<double>::max();
  b_phi_intercept_shifted = b.phi_intercept;

  const unsigned int region = region_index(b.region);
  if (a.side != b.side || a.last_layer >= b.first_layer)
  {
    return false;
  }

  const unsigned int gap = b.first_layer - a.last_layer - 1;
  if (gap > m_connectMaxLayerGap[region])
  {
    return false;
  }

  const double ra = m_idealPadMap->get_radius(a.last_layer);
  const double rb = m_idealPadMap->get_radius(b.first_layer);
  if (!std::isfinite(ra) || !std::isfinite(rb) || ra <= 0.0 || rb <= 0.0)
  {
    return false;
  }

  // Reliability-weighted comparison point.  A longer candidate is propagated
  // farther across the gap, so the match point moves closer to the short piece.
  const double rmatch = reliability_weighted_coordinate(
      ra, rb, a.nblobs, b.nblobs);
  const double slope_rel = piece_slope_reliability(a.nblobs, b.nblobs);
  const double slope_cut_scale = 1.0 / std::max(slope_rel, 0.20);

  const double tbin_a = a.tbin_slope_r * rmatch + a.tbin_intercept_r;
  const double tbin_b = b.tbin_slope * rmatch + b.tbin_intercept;
  const double dtbin = std::fabs(tbin_a - tbin_b);
  if (dtbin > m_connect_dtbin[region])
  {
    return false;
  }

  const double dmtbin = std::fabs(a.tbin_slope_r - b.tbin_slope);
  if (dmtbin > slope_cut_scale * m_connect_dtbin_slope[region])
  {
    return false;
  }

  if (m_linearPhiPrecutScale > 0.0)
  {
    const double phi_a_linear = a.phi_slope * rmatch + a.phi_intercept;
    const double phi_b_linear_raw = b.phi_slope * rmatch + b.phi_intercept;
    const double phi_b_linear = unwrap_phi_to_reference(phi_b_linear_raw, phi_a_linear);
    const double dphi_linear = std::fabs(phi_a_linear - phi_b_linear);
    const double dmphi_linear = std::fabs(a.phi_slope - b.phi_slope);

    if (dphi_linear > m_linearPhiPrecutScale * m_connect_dphi[region] ||
        dmphi_linear > m_linearPhiPrecutScale * slope_cut_scale * m_connect_dphi_slope[region])
    {
      return false;
    }
  }

  const auto phi_a_pair = predict_phi_and_slope(a, rmatch);
  const auto phi_b_pair = predict_phi_and_slope(b, rmatch);
  const double phi_b = unwrap_phi_to_reference(phi_b_pair.first, phi_a_pair.first);

  b_phi_intercept_shifted = b.phi_intercept + (phi_b - phi_b_pair.first);

  const double dphi = std::fabs(phi_a_pair.first - phi_b);
  const double dmphi = std::fabs(phi_a_pair.second - phi_b_pair.second);
  if (dphi > m_connect_dphi[region] ||
      dmphi > slope_cut_scale * m_connect_dphi_slope[region])
  {
    return false;
  }

  if (m_doDebugHistograms)
  {
    std::lock_guard<std::mutex> lock(m_debugMutex);
    if (m_h_dphi[region]) m_h_dphi[region]->Fill(dphi);
    if (m_h_dtbin[region]) m_h_dtbin[region]->Fill(dtbin);
    if (m_h_dmphi[region]) m_h_dmphi[region]->Fill(dmphi);
    if (m_h_dmtbin[region]) m_h_dmtbin[region]->Fill(dmtbin);
    if (m_h_dphi_vs_dtbin[region]) m_h_dphi_vs_dtbin[region]->Fill(dphi, dtbin);
    if (m_h_dmphi_vs_dmtbin[region]) m_h_dmphi_vs_dmtbin[region]->Fill(dmphi, dmtbin);
    if (m_h_dphi_vs_dmphi[region]) m_h_dphi_vs_dmphi[region]->Fill(dphi, dmphi);
    if (m_h_tbin_slope_vs_last_tbin[region])
    {
      m_h_tbin_slope_vs_last_tbin[region]->Fill(
          a.tbin_slope_r * ra + a.tbin_intercept_r, a.tbin_slope_r);
    }
    if (m_h_tbin_slope_vs_first_tbin[region] && !b.points.empty())
    {
      m_h_tbin_slope_vs_first_tbin[region]->Fill(b.points.front().tbin, b.tbin_slope);
    }
  }

  score =
      m_weight_phi[region] * sqr(dphi / std::max(m_score_dphi[region], 1.0e-9)) +
      m_weight_tbin[region] * sqr(dtbin / std::max(m_score_dtbin[region], 1.0e-9)) +
      slope_rel * m_weight_phi_slope[region] *
          sqr(dmphi / std::max(m_score_dphi_slope[region], 1.0e-9)) +
      slope_rel * m_weight_tbin_slope[region] *
          sqr(dmtbin / std::max(m_score_dtbin_slope[region], 1.0e-9)) +
      0.05 * static_cast<double>(gap);

  return true;
}

void Tpc_AssembledTrackReco::connect_sector_pieces(const std::vector<Piece>& pieces, int side, unsigned int sector, std::vector<Candidate>& output) const
{
  std::vector<unsigned int> order;
  for (unsigned int i = 0; i < pieces.size(); ++i)
  {
    if (pieces[i].side == side && pieces[i].sector == sector &&
        pieces[i].nblobs >= m_robustPieceMinBlobs)
    {
      order.push_back(i);
    }
  }

  // Only robust pieces are allowed to pre-assemble inside a sector.  Small
  // 2/3/4-layer fragments are emitted as single candidates below and are
  // collected later by the global long-first pass.
  std::sort(order.begin(), order.end(), PieceStartSort(&pieces));
  std::vector<uint8_t> used(pieces.size(), 0);

  for (unsigned int io = 0; io < order.size(); ++io)
  {
    const unsigned int iseed = order[io];
    if (used[iseed])
    {
      continue;
    }

    std::vector<unsigned int> current_indices;
    current_indices.push_back(iseed);
    used[iseed] = 1;

    Candidate current;
    if (!refit_candidate(pieces, current_indices, current))
    {
      continue;
    }

    bool merged_any = true;
    while (merged_any)
    {
      merged_any = false;
      int best_j = -1;
      double best_score = std::numeric_limits<double>::max();

      for (unsigned int j : order)
      {
        if (used[j])
        {
          continue;
        }

        Candidate piece_candidate;
        if (!refit_candidate(pieces, std::vector<unsigned int>{j}, piece_candidate))
        {
          continue;
        }

        double score = 0.0;
        if (!candidates_can_connect(current, piece_candidate, score, true, 1.0))
        {
          continue;
        }

        if (score < best_score)
        {
          best_score = score;
          best_j = static_cast<int>(j);
        }
      }

      if (best_j >= 0)
      {
        std::vector<unsigned int> trial_indices = current_indices;
        trial_indices.push_back(static_cast<unsigned int>(best_j));

        Candidate refit;
        if (refit_candidate(pieces, trial_indices, refit))
        {
          const Piece& accepted_piece = pieces[static_cast<unsigned int>(best_j)];
          const unsigned int accepted_gap =
              minimum_layer_mask_gap(current.layer_mask, accepted_piece.layer_mask);
          const unsigned int qa_region = region_index(
              current.first_layer <= accepted_piece.first_layer
                  ? accepted_piece.region
                  : current.first_region);
          if (m_doDebugHistograms)
          {
            std::lock_guard<std::mutex> lock(m_debugMutex);
            if (m_h_score[qa_region])
            {
              m_h_score[qa_region]->Fill(best_score);
            }
            if (m_h_layer_gap[qa_region])
            {
              m_h_layer_gap[qa_region]->Fill(static_cast<double>(accepted_gap));
            }
            if (m_h_matched_sector_delta[qa_region])
            {
              m_h_matched_sector_delta[qa_region]->Fill(static_cast<double>(wrapped_sector_delta(current.last_sector, accepted_piece.sector)));
            }
          }

          current = refit;
          current_indices.swap(trial_indices);
          used[best_j] = 1;
          merged_any = true;
        }
      }
    }

    output.push_back(current);
  }

  for (unsigned int i = 0; i < pieces.size(); ++i)
  {
    if (pieces[i].side != side || pieces[i].sector != sector ||
        pieces[i].nblobs >= m_robustPieceMinBlobs)
    {
      continue;
    }
    Candidate single;
    if (refit_candidate(pieces, std::vector<unsigned int>{i}, single))
    {
      output.push_back(std::move(single));
    }
  }
}

bool Tpc_AssembledTrackReco::fast_piece_relation(const Candidate& a,
                                                    const Piece& b,
                                                    bool& broad,
                                                    bool& tight,
                                                    double& score) const
{
  broad = false;
  tight = false;
  score = std::numeric_limits<double>::max();

  const unsigned int region = region_index(b.region);
  if (a.side != b.side || a.last_layer >= b.first_layer)
  {
    return false;
  }

  const unsigned int gap = b.first_layer - a.last_layer - 1;
  if (gap > m_connectMaxLayerGap[region])
  {
    return false;
  }

  const double ra = m_idealPadMap->get_radius(a.last_layer);
  const double rb = m_idealPadMap->get_radius(b.first_layer);
  if (!std::isfinite(ra) || !std::isfinite(rb) || ra <= 0.0 || rb <= 0.0)
  {
    return false;
  }

  const double rmatch = reliability_weighted_coordinate(
      ra, rb, a.nblobs, b.nblobs);
  const double slope_rel = piece_slope_reliability(a.nblobs, b.nblobs);
  const double slope_cut_scale = 1.0 / std::max(slope_rel, 0.20);

  const double tbin_a = a.tbin_slope_r * rmatch + a.tbin_intercept_r;
  const double tbin_b = b.tbin_slope * rmatch + b.tbin_intercept;
  const double dtbin = std::fabs(tbin_a - tbin_b);
  const double dmtbin = std::fabs(a.tbin_slope_r - b.tbin_slope);

  const double phi_a = a.phi_slope * rmatch + a.phi_intercept;
  const double phi_b_raw = b.phi_slope * rmatch + b.phi_intercept;
  const double phi_b = unwrap_phi_to_reference(phi_b_raw, phi_a);
  const double dphi = std::fabs(phi_a - phi_b);
  const double dmphi = std::fabs(a.phi_slope - b.phi_slope);

  broad = dtbin <= m_fastBroadScale * m_connect_dtbin[region] &&
          dmtbin <= m_fastBroadScale * slope_cut_scale * m_connect_dtbin_slope[region] &&
          dphi <= m_fastBroadScale * m_connect_dphi[region] &&
          dmphi <= m_fastBroadScale * slope_cut_scale * m_connect_dphi_slope[region];

  if (!broad)
  {
    return true;
  }

  tight = dtbin <= m_fastTightScale * m_connect_dtbin[region] &&
          dmtbin <= m_fastTightScale * slope_cut_scale * m_connect_dtbin_slope[region] &&
          dphi <= m_fastTightScale * m_connect_dphi[region] &&
          dmphi <= m_fastTightScale * slope_cut_scale * m_connect_dphi_slope[region];

  score =
      sqr(dphi / std::max(m_fastBroadScale * m_connect_dphi[region], 1.0e-9)) +
      sqr(dtbin / std::max(m_fastBroadScale * m_connect_dtbin[region], 1.0e-9)) +
      slope_rel * sqr(dmphi / std::max(m_fastBroadScale * m_connect_dphi_slope[region], 1.0e-9)) +
      slope_rel * sqr(dmtbin / std::max(m_fastBroadScale * m_connect_dtbin_slope[region], 1.0e-9)) +
      0.05 * static_cast<double>(gap);
  return true;
}

void Tpc_AssembledTrackReco::connect_sector_pieces_fast(const std::vector<Piece>& pieces,
                                                         int side,
                                                         unsigned int sector,
                                                         std::vector<Candidate>& output) const
{
  std::vector<unsigned int> order;
  for (unsigned int i = 0; i < pieces.size(); ++i)
  {
    if (pieces[i].side == side && pieces[i].sector == sector &&
        pieces[i].nblobs >= m_robustPieceMinBlobs)
    {
      order.push_back(i);
    }
  }

  // Only robust pieces are allowed to pre-assemble inside a sector.  Small
  // 2/3/4-layer fragments are emitted as single candidates below and are
  // collected later by the global long-first pass.
  std::sort(order.begin(), order.end(), PieceStartSort(&pieces));
  std::vector<uint8_t> used(pieces.size(), 0);

  for (unsigned int iseed : order)
  {
    if (used[iseed])
    {
      continue;
    }

    std::vector<unsigned int> current_indices{iseed};
    used[iseed] = 1;

    Candidate current;
    if (!refit_candidate(pieces, current_indices, current))
    {
      continue;
    }

    while (true)
    {
      unsigned int broad_count = 0;
      int unique_piece = -1;
      bool unique_tight = false;
      double unique_score = std::numeric_limits<double>::max();

      for (unsigned int j = 0; j < pieces.size(); ++j)
      {
        const Piece& candidate_piece = pieces[j];
        if (candidate_piece.side != side ||
            candidate_piece.nblobs < m_robustPieceMinBlobs)
        {
          continue;
        }
        if (wrapped_sector_delta(sector, candidate_piece.sector) < -1 ||
            wrapped_sector_delta(sector, candidate_piece.sector) > 1)
        {
          continue;
        }
        if (candidate_piece.sector == sector && used[j])
        {
          continue;
        }

        bool broad = false;
        bool tight = false;
        double score = 0.0;
        if (!fast_piece_relation(current, candidate_piece, broad, tight, score) || !broad)
        {
          continue;
        }

        ++broad_count;
        if (broad_count == 1)
        {
          unique_piece = static_cast<int>(j);
          unique_tight = tight;
          unique_score = score;
        }
      }

      if (broad_count != 1 || unique_piece < 0 || !unique_tight)
      {
        break;
      }

      const unsigned int accepted_index = static_cast<unsigned int>(unique_piece);
      if (pieces[accepted_index].sector != sector || used[accepted_index])
      {
        // A neighboring-sector continuation is plausible, so defer the whole
        // decision to the global side stage where both options can compete.
        break;
      }

      Candidate accepted_candidate;
      double exact_score = 0.0;
      if (!refit_candidate(pieces, std::vector<unsigned int>{accepted_index}, accepted_candidate) ||
          !candidates_can_connect(current, accepted_candidate, exact_score, true, 1.0))
      {
        break;
      }
      unique_score = exact_score;

      std::vector<unsigned int> trial_indices = current_indices;
      trial_indices.push_back(accepted_index);
      Candidate refit;
      if (!refit_candidate(pieces, trial_indices, refit))
      {
        break;
      }

      if (m_doDebugHistograms)
      {
        const unsigned int qa_region = region_index(pieces[accepted_index].region);
        std::lock_guard<std::mutex> lock(m_debugMutex);
        if (m_h_score[qa_region]) m_h_score[qa_region]->Fill(unique_score);
        if (m_h_layer_gap[qa_region])
        {
          const unsigned int gap = pieces[accepted_index].first_layer - current.last_layer - 1;
          m_h_layer_gap[qa_region]->Fill(static_cast<double>(gap));
        }
        if (m_h_matched_sector_delta[qa_region]) m_h_matched_sector_delta[qa_region]->Fill(0.0);
      }

      current = std::move(refit);
      current_indices.swap(trial_indices);
      used[accepted_index] = 1;
    }

    output.push_back(std::move(current));
  }

  for (unsigned int i = 0; i < pieces.size(); ++i)
  {
    if (pieces[i].side != side || pieces[i].sector != sector ||
        pieces[i].nblobs >= m_robustPieceMinBlobs)
    {
      continue;
    }
    Candidate single;
    if (refit_candidate(pieces, std::vector<unsigned int>{i}, single))
    {
      output.push_back(std::move(single));
    }
  }
}

bool Tpc_AssembledTrackReco::cheap_candidate_relation(const Candidate& a,
                                                         const Candidate& b,
                                                         bool allow_same_sector,
                                                         double window_scale,
                                                         bool& broad,
                                                         bool& tight,
                                                         double& score) const
{
  broad = false;
  tight = false;
  score = std::numeric_limits<double>::max();

  if (a.side != b.side)
  {
    return false;
  }
  const bool same_single_sector =
      a.first_sector == a.last_sector &&
      b.first_sector == b.last_sector &&
      a.first_sector == b.first_sector;
  if (!allow_same_sector && same_single_sector)
  {
    return false;
  }

  const unsigned int gap = minimum_layer_mask_gap(a.layer_mask, b.layer_mask);
  if (gap == std::numeric_limits<unsigned int>::max())
  {
    return false;
  }

  const unsigned int region = region_index(
      a.first_layer <= b.first_layer ? b.first_region : a.first_region);
  if (gap > m_connectMaxLayerGap[region])
  {
    return false;
  }

  const double slope_rel = piece_slope_reliability(a.nblobs, b.nblobs);

  auto phi_linear = [](const Candidate& c, double r)
  {
    return c.phi_slope * r + c.phi_intercept;
  };
  auto tbin_at = [](const Candidate& c, double r)
  {
    return c.tbin_slope_r * r + c.tbin_intercept_r;
  };

  double dphi0 = 0.0;
  double dt0 = 0.0;
  double dphi1 = 0.0;
  double dt1 = 0.0;
  double w0 = 1.0;
  double w1 = 1.0;
  double cut_scale0 = 1.0;
  double cut_scale1 = 1.0;

  if (a.last_layer < b.first_layer || b.last_layer < a.first_layer)
  {
    const Candidate& inner = a.last_layer < b.first_layer ? a : b;
    const Candidate& outer = a.last_layer < b.first_layer ? b : a;

    const double r_inner = m_idealPadMap->get_radius(inner.last_layer);
    const double r_outer = m_idealPadMap->get_radius(outer.first_layer);
    if (!std::isfinite(r_outer) || !std::isfinite(r_inner) ||
        r_outer <= 0.0 || r_inner <= 0.0)
    {
      return false;
    }

    const double rmatch = reliability_weighted_coordinate(
        r_inner, r_outer, inner.nblobs, outer.nblobs);
    const double phi_inner = phi_linear(inner, rmatch);
    const double phi_outer =
        unwrap_phi_to_reference(phi_linear(outer, rmatch), phi_inner);

    dphi0 = std::fabs(phi_inner - phi_outer);
    dt0 = std::fabs(tbin_at(inner, rmatch) - tbin_at(outer, rmatch));
    dphi1 = dphi0;
    dt1 = dt0;
    w0 = w1 = 1.0;
  }
  else
  {
    // A fragment can sit inside a radial gap of an already assembled long
    // candidate.  Use the more reliable candidate as the road and test the
    // fragment at both of its measured edges.
    const Candidate& reference =
        (a.nblobs > b.nblobs ||
         (a.nblobs == b.nblobs &&
          (a.last_layer - a.first_layer) >= (b.last_layer - b.first_layer)))
            ? a
            : b;
    const Candidate& fragment = &reference == &a ? b : a;
    const double r0 = m_idealPadMap->get_radius(fragment.first_layer);
    const double r1 = m_idealPadMap->get_radius(fragment.last_layer);
    if (!std::isfinite(r0) || !std::isfinite(r1) || r0 <= 0.0 || r1 <= 0.0)
    {
      return false;
    }

    const double phi_ref0 = phi_linear(reference, r0);
    const double phi_frag0 = unwrap_phi_to_reference(phi_linear(fragment, r0), phi_ref0);
    const double phi_ref1 = phi_linear(reference, r1);
    const double phi_frag1 = unwrap_phi_to_reference(phi_linear(fragment, r1), phi_ref1);
    dphi0 = std::fabs(phi_ref0 - phi_frag0);
    dt0 = std::fabs(tbin_at(reference, r0) - tbin_at(fragment, r0));
    dphi1 = std::fabs(phi_ref1 - phi_frag1);
    dt1 = std::fabs(tbin_at(reference, r1) - tbin_at(fragment, r1));
    w0 = w1 = 1.0;
    cut_scale0 = cut_scale1 = 1.25;
  }

  const double dmphi = std::fabs(a.phi_slope - b.phi_slope);
  const double dmtbin = std::fabs(a.tbin_slope_r - b.tbin_slope_r);
  const double slope_cut_scale = 1.0 / std::max(slope_rel, 0.20);
  const double phi_pre_scale = m_linearPhiPrecutScale > 0.0
                                   ? m_linearPhiPrecutScale
                                   : 1.0;

  const bool phi_broad = m_linearPhiPrecutScale <= 0.0 ||
      (dphi0 <= window_scale * phi_pre_scale * cut_scale0 * m_connect_dphi[region] &&
       dphi1 <= window_scale * phi_pre_scale * cut_scale1 * m_connect_dphi[region] &&
       dmphi <= window_scale * phi_pre_scale * slope_cut_scale * m_connect_dphi_slope[region]);

  broad = phi_broad &&
          dt0 <= window_scale * cut_scale0 * m_connect_dtbin[region] &&
          dt1 <= window_scale * cut_scale1 * m_connect_dtbin[region] &&
          dmtbin <= window_scale * slope_cut_scale * m_connect_dtbin_slope[region];
  if (!broad)
  {
    return true;
  }

  const double wsum = std::max(w0 + w1, 1.0e-9);
  const double eff_dphi = std::sqrt((w0 * dphi0 * dphi0 + w1 * dphi1 * dphi1) / wsum);
  const double eff_dt = std::sqrt((w0 * dt0 * dt0 + w1 * dt1 * dt1) / wsum);

  tight = eff_dphi <= window_scale * m_normalTightScale * m_connect_dphi[region] &&
          eff_dt <= window_scale * m_normalTightScale * m_connect_dtbin[region] &&
          dmphi <= window_scale * m_normalTightScale * slope_cut_scale * m_connect_dphi_slope[region] &&
          dmtbin <= window_scale * m_normalTightScale * slope_cut_scale * m_connect_dtbin_slope[region];

  score =
      m_weight_phi[region] * sqr(eff_dphi / std::max(window_scale * m_score_dphi[region], 1.0e-9)) +
      m_weight_tbin[region] * sqr(eff_dt / std::max(window_scale * m_score_dtbin[region], 1.0e-9)) +
      slope_rel * m_weight_phi_slope[region] *
          sqr(dmphi / std::max(window_scale * m_score_dphi_slope[region], 1.0e-9)) +
      slope_rel * m_weight_tbin_slope[region] *
          sqr(dmtbin / std::max(window_scale * m_score_dtbin_slope[region], 1.0e-9)) +
      0.05 * static_cast<double>(gap);
  return true;
}

bool Tpc_AssembledTrackReco::candidates_can_connect(const Candidate& a,
                                                      const Candidate& b,
                                                      double& score,
                                                      bool allow_same_sector,
                                                      double window_scale) const
{
  score = std::numeric_limits<double>::max();

  bool broad = false;
  bool tight = false;
  double cheap_score = 0.0;
  if (!cheap_candidate_relation(a, b, allow_same_sector, window_scale,
                                broad, tight, cheap_score) || !broad)
  {
    return false;
  }

  const unsigned int gap = minimum_layer_mask_gap(a.layer_mask, b.layer_mask);
  if (gap == std::numeric_limits<unsigned int>::max())
  {
    return false;
  }
  const unsigned int region = region_index(
      a.first_layer <= b.first_layer ? b.first_region : a.first_region);

  const double slope_rel = piece_slope_reliability(a.nblobs, b.nblobs);

  auto tbin_at = [](const Candidate& c, double r)
  {
    return c.tbin_slope_r * r + c.tbin_intercept_r;
  };

  double dphi0 = 0.0;
  double dt0 = 0.0;
  double dphi1 = 0.0;
  double dt1 = 0.0;
  double w0 = 1.0;
  double w1 = 1.0;
  double cut_scale0 = 1.0;
  double cut_scale1 = 1.0;
  double dmphi = 0.0;

  if (a.last_layer < b.first_layer || b.last_layer < a.first_layer)
  {
    const Candidate& inner = a.last_layer < b.first_layer ? a : b;
    const Candidate& outer = a.last_layer < b.first_layer ? b : a;

    const double r_inner = m_idealPadMap->get_radius(inner.last_layer);
    const double r_outer = m_idealPadMap->get_radius(outer.first_layer);
    if (!std::isfinite(r_outer) || !std::isfinite(r_inner) ||
        r_outer <= 0.0 || r_inner <= 0.0)
    {
      return false;
    }

    const double rmatch = reliability_weighted_coordinate(
        r_inner, r_outer, inner.nblobs, outer.nblobs);
    const auto inner_at_match = predict_phi_and_slope(inner, rmatch);
    const auto outer_at_match_raw = predict_phi_and_slope(outer, rmatch);
    const double outer_phi_match =
        unwrap_phi_to_reference(outer_at_match_raw.first, inner_at_match.first);

    dphi0 = std::fabs(inner_at_match.first - outer_phi_match);
    dt0 = std::fabs(tbin_at(inner, rmatch) - tbin_at(outer, rmatch));
    dphi1 = dphi0;
    dt1 = dt0;
    dmphi = std::fabs(inner_at_match.second - outer_at_match_raw.second);
    w0 = w1 = 1.0;
  }
  else
  {
    const Candidate& reference =
        (a.nblobs > b.nblobs ||
         (a.nblobs == b.nblobs &&
          (a.last_layer - a.first_layer) >= (b.last_layer - b.first_layer)))
            ? a
            : b;
    const Candidate& fragment = &reference == &a ? b : a;
    const double r0 = m_idealPadMap->get_radius(fragment.first_layer);
    const double r1 = m_idealPadMap->get_radius(fragment.last_layer);
    if (!std::isfinite(r0) || !std::isfinite(r1) || r0 <= 0.0 || r1 <= 0.0)
    {
      return false;
    }

    const auto ref0 = predict_phi_and_slope(reference, r0);
    const auto frag0raw = predict_phi_and_slope(fragment, r0);
    const double fragphi0 = unwrap_phi_to_reference(frag0raw.first, ref0.first);
    const auto ref1 = predict_phi_and_slope(reference, r1);
    const auto frag1raw = predict_phi_and_slope(fragment, r1);
    const double fragphi1 = unwrap_phi_to_reference(frag1raw.first, ref1.first);
    dphi0 = std::fabs(ref0.first - fragphi0);
    dt0 = std::fabs(tbin_at(reference, r0) - tbin_at(fragment, r0));
    dphi1 = std::fabs(ref1.first - fragphi1);
    dt1 = std::fabs(tbin_at(reference, r1) - tbin_at(fragment, r1));
    dmphi = 0.5 * (std::fabs(ref0.second - frag0raw.second) +
                   std::fabs(ref1.second - frag1raw.second));
    w0 = w1 = 1.0;
    cut_scale0 = cut_scale1 = 1.25;
  }

  const double dmtbin = std::fabs(a.tbin_slope_r - b.tbin_slope_r);
  const double slope_cut_scale = 1.0 / std::max(slope_rel, 0.20);
  if (dphi0 > window_scale * cut_scale0 * m_connect_dphi[region] ||
      dt0 > window_scale * cut_scale0 * m_connect_dtbin[region] ||
      dphi1 > window_scale * cut_scale1 * m_connect_dphi[region] ||
      dt1 > window_scale * cut_scale1 * m_connect_dtbin[region] ||
      dmphi > window_scale * slope_cut_scale * m_connect_dphi_slope[region] ||
      dmtbin > window_scale * slope_cut_scale * m_connect_dtbin_slope[region])
  {
    return false;
  }

  const double wsum = std::max(w0 + w1, 1.0e-9);
  const double eff_dphi = std::sqrt((w0 * dphi0 * dphi0 + w1 * dphi1 * dphi1) / wsum);
  const double eff_dt = std::sqrt((w0 * dt0 * dt0 + w1 * dt1 * dt1) / wsum);

  if (m_doDebugHistograms)
  {
    std::lock_guard<std::mutex> lock(m_debugMutex);
    if (m_h_dphi[region]) m_h_dphi[region]->Fill(eff_dphi);
    if (m_h_dtbin[region]) m_h_dtbin[region]->Fill(eff_dt);
    if (m_h_dmphi[region]) m_h_dmphi[region]->Fill(dmphi);
    if (m_h_dmtbin[region]) m_h_dmtbin[region]->Fill(dmtbin);
    if (m_h_dphi_vs_dtbin[region]) m_h_dphi_vs_dtbin[region]->Fill(eff_dphi, eff_dt);
    if (m_h_dmphi_vs_dmtbin[region]) m_h_dmphi_vs_dmtbin[region]->Fill(dmphi, dmtbin);
    if (m_h_dphi_vs_dmphi[region]) m_h_dphi_vs_dmphi[region]->Fill(eff_dphi, dmphi);

    const double ra_last = m_idealPadMap->get_radius(a.last_layer);
    const double rb_first = m_idealPadMap->get_radius(b.first_layer);
    if (m_h_tbin_slope_vs_last_tbin[region] && std::isfinite(ra_last))
    {
      m_h_tbin_slope_vs_last_tbin[region]->Fill(
          a.tbin_slope_r * ra_last + a.tbin_intercept_r, a.tbin_slope_r);
    }
    if (m_h_tbin_slope_vs_first_tbin[region] && std::isfinite(rb_first))
    {
      m_h_tbin_slope_vs_first_tbin[region]->Fill(
          b.tbin_slope_r * rb_first + b.tbin_intercept_r, b.tbin_slope_r);
    }
  }

  score =
      m_weight_phi[region] * sqr(eff_dphi / std::max(window_scale * m_score_dphi[region], 1.0e-9)) +
      m_weight_tbin[region] * sqr(eff_dt / std::max(window_scale * m_score_dtbin[region], 1.0e-9)) +
      slope_rel * m_weight_phi_slope[region] *
          sqr(dmphi / std::max(window_scale * m_score_dphi_slope[region], 1.0e-9)) +
      slope_rel * m_weight_tbin_slope[region] *
          sqr(dmtbin / std::max(window_scale * m_score_dtbin_slope[region], 1.0e-9)) +
      0.05 * static_cast<double>(gap);

  return true;
}

void Tpc_AssembledTrackReco::connect_side_candidates(const std::vector<Piece>& pieces,
                                                     const std::vector<Candidate>& seeds,
                                                     int side,
                                                     std::vector<Candidate>& output,
                                                     bool allow_same_sector,
                                                     double window_scale,
                                                     bool enable_unique_tight) const
{
  std::vector<unsigned int> order;
  for (unsigned int i = 0; i < seeds.size(); ++i)
  {
    if (seeds[i].side == side)
    {
      order.push_back(i);
    }
  }
  if (order.empty())
  {
    return;
  }

  std::sort(order.begin(), order.end(), CandidateStartSort(&seeds));
  std::vector<uint8_t> used(seeds.size(), 0);

  for (unsigned int iseed : order)
  {
    if (used[iseed])
    {
      continue;
    }

    std::vector<unsigned int> current_indices = seeds[iseed].piece_indices;
    used[iseed] = 1;

    Candidate current;
    if (!refit_candidate(pieces, current_indices, current))
    {
      continue;
    }

    while (true)
    {
      int best_j = -1;
      double best_score = std::numeric_limits<double>::max();

      // First connect established/robust candidates to one another.  Only
      // after no >=m_robustPieceMinBlobs continuation survives do we let the
      // current track collect short dead-layer fragments.
      for (unsigned int j : order)
      {
        if (used[j] || seeds[j].nblobs < m_robustPieceMinBlobs)
        {
          continue;
        }
        double score = 0.0;
        if (!candidates_can_connect(current, seeds[j], score,
                                    allow_same_sector, window_scale))
        {
          continue;
        }
        if (score < best_score)
        {
          best_score = score;
          best_j = static_cast<int>(j);
        }
      }

      if (best_j < 0 && enable_unique_tight)
      {
        unsigned int broad_count = 0;
        int unique_j = -1;
        bool unique_tight = false;

        for (unsigned int j : order)
        {
          if (used[j] || seeds[j].nblobs >= m_robustPieceMinBlobs)
          {
            continue;
          }
          bool broad = false;
          bool tight = false;
          double cheap_score = 0.0;
          if (!cheap_candidate_relation(current, seeds[j], allow_same_sector,
                                        window_scale, broad, tight, cheap_score) ||
              !broad)
          {
            continue;
          }
          ++broad_count;
          if (broad_count == 1)
          {
            unique_j = static_cast<int>(j);
            unique_tight = tight;
          }
        }

        if (broad_count == 1 && unique_j >= 0 && unique_tight)
        {
          // Cheap matching is only a preselection.  Every accepted connection
          // gets the same exact score definition so h_score and reassociation
          // are directly comparable.
          double exact_score = 0.0;
          if (candidates_can_connect(current,
                                     seeds[static_cast<unsigned int>(unique_j)],
                                     exact_score, allow_same_sector, window_scale))
          {
            best_j = unique_j;
            best_score = exact_score;
          }
        }
      }

      // Ambiguous/non-tight cases use the exact sagitta matcher. If the cheap
      // unique path did not select anything this is exactly the old behavior
      // apart from optional same-sector participation and window scaling.
      if (best_j < 0)
      {
        for (unsigned int j : order)
        {
          if (used[j] || seeds[j].nblobs >= m_robustPieceMinBlobs)
          {
            continue;
          }

          double score = 0.0;
          if (!candidates_can_connect(current, seeds[j], score,
                                      allow_same_sector, window_scale))
          {
            continue;
          }

          if (score < best_score)
          {
            best_score = score;
            best_j = static_cast<int>(j);
          }
        }
      }

      if (best_j < 0)
      {
        break;
      }

      std::vector<unsigned int> trial_indices = current_indices;
      const Candidate& accepted_seed = seeds[static_cast<unsigned int>(best_j)];
      trial_indices.insert(trial_indices.end(),
                           accepted_seed.piece_indices.begin(),
                           accepted_seed.piece_indices.end());

      Candidate refit;
      if (!refit_candidate(pieces, trial_indices, refit))
      {
        break;
      }

      const unsigned int accepted_gap =
          minimum_layer_mask_gap(current.layer_mask, accepted_seed.layer_mask);
      const unsigned int qa_region = region_index(
          current.first_layer <= accepted_seed.first_layer
              ? accepted_seed.first_region
              : current.first_region);
      if (m_doDebugHistograms)
      {
        std::lock_guard<std::mutex> lock(m_debugMutex);
        if (m_h_score[qa_region]) m_h_score[qa_region]->Fill(best_score);
        if (m_h_layer_gap[qa_region]) m_h_layer_gap[qa_region]->Fill(static_cast<double>(accepted_gap));
        if (m_h_matched_sector_delta[qa_region])
        {
          m_h_matched_sector_delta[qa_region]->Fill(
              static_cast<double>(wrapped_sector_delta(current.last_sector,
                                                       accepted_seed.first_sector)));
        }
      }

      current = std::move(refit);
      current_indices.swap(trial_indices);
      used[static_cast<unsigned int>(best_j)] = 1;
    }

    output.push_back(std::move(current));
  }
}

void Tpc_AssembledTrackReco::rescue_side_candidates(const std::vector<Piece>& pieces,
                                                     const std::vector<Candidate>& seeds,
                                                     int side,
                                                     double window_scale,
                                                     std::vector<Candidate>& output) const
{
  std::vector<unsigned int> order;
  for (unsigned int i = 0; i < seeds.size(); ++i)
  {
    if (seeds[i].side == side)
    {
      order.push_back(i);
    }
  }
  if (order.empty())
  {
    return;
  }

  std::sort(order.begin(), order.end(), CandidateStartSort(&seeds));
  std::vector<uint8_t> used(seeds.size(), 0);

  for (unsigned int iseed : order)
  {
    if (used[iseed])
    {
      continue;
    }

    std::vector<unsigned int> current_indices = seeds[iseed].piece_indices;
    used[iseed] = 1;
    Candidate current;
    if (!refit_candidate(pieces, current_indices, current))
    {
      continue;
    }

    while (true)
    {
      int best_j = -1;
      double best_score = std::numeric_limits<double>::max();

      for (unsigned int j : order)
      {
        if (used[j])
        {
          continue;
        }

        const Candidate& candidate = seeds[j];
        // Rescue is intentionally a small-fragment stage.  Robust+robust
        // connections must have been decided by the normal exact pass; do not
        // join two established tracks merely because the rescue window is wider.
        if (current.nblobs >= m_robustPieceMinBlobs &&
            candidate.nblobs >= m_robustPieceMinBlobs)
        {
          continue;
        }
        if (current.nsegments > m_rescueMaxSegments &&
            candidate.nsegments > m_rescueMaxSegments)
        {
          continue;
        }

        double score = 0.0;
        if (!candidates_can_connect(current, candidate, score,
                                    true, window_scale))
        {
          continue;
        }
        if (score < best_score)
        {
          best_score = score;
          best_j = static_cast<int>(j);
        }
      }

      if (best_j < 0)
      {
        break;
      }

      const Candidate& accepted = seeds[static_cast<unsigned int>(best_j)];
      std::vector<unsigned int> trial_indices = current_indices;
      trial_indices.insert(trial_indices.end(),
                           accepted.piece_indices.begin(),
                           accepted.piece_indices.end());
      Candidate refit;
      if (!refit_candidate(pieces, trial_indices, refit))
      {
        break;
      }
      current = std::move(refit);
      current_indices.swap(trial_indices);
      used[static_cast<unsigned int>(best_j)] = 1;
    }

    output.push_back(std::move(current));
  }
}

void Tpc_AssembledTrackReco::final_reassociate(const std::vector<Piece>& pieces,
                                                  int side,
                                                  std::vector<Candidate>& tracks) const
{
  if (!m_enableFinalReassociation || m_finalReassociationIterations == 0 ||
      tracks.size() < 2)
  {
    return;
  }

  for (unsigned int iter = 0; iter < m_finalReassociationIterations; ++iter)
  {
    std::vector<int> owner(pieces.size(), -1);
    for (unsigned int it = 0; it < tracks.size(); ++it)
    {
      if (tracks[it].side != side) continue;
      for (unsigned int ip : tracks[it].piece_indices)
      {
        if (ip < owner.size()) owner[ip] = static_cast<int>(it);
      }
    }

    bool found = false;
    unsigned int best_piece = 0;
    unsigned int best_donor = 0;
    unsigned int best_target = 0;
    double best_ratio = std::numeric_limits<double>::max();
    double best_target_score = std::numeric_limits<double>::max();

    for (unsigned int ip = 0; ip < pieces.size(); ++ip)
    {
      const int donor_i = owner[ip];
      if (donor_i < 0 || pieces[ip].side != side)
      {
        continue;
      }
      const unsigned int donor = static_cast<unsigned int>(donor_i);
      if (tracks[donor].piece_indices.size() <= 1U)
      {
        continue;
      }

      Candidate single;
      if (!refit_candidate(pieces, std::vector<unsigned int>{ip}, single))
      {
        continue;
      }

      std::vector<unsigned int> donor_indices;
      donor_indices.reserve(tracks[donor].piece_indices.size() - 1U);
      for (unsigned int old_ip : tracks[donor].piece_indices)
      {
        if (old_ip != ip) donor_indices.push_back(old_ip);
      }
      Candidate donor_without;
      if (!refit_candidate(pieces, donor_indices, donor_without))
      {
        continue;
      }

      double donor_score = std::numeric_limits<double>::max();
      const bool donor_relation_ok =
          candidates_can_connect(donor_without, single, donor_score, true, 1.0);

      for (unsigned int target = 0; target < tracks.size(); ++target)
      {
        if (target == donor || tracks[target].side != side)
        {
          continue;
        }

        double target_score = 0.0;
        if (!candidates_can_connect(tracks[target], single,
                                    target_score, true, 1.0))
        {
          continue;
        }

        if (donor_relation_ok)
        {
          if (!(target_score < m_finalReassociationImprovementRatio * donor_score))
          {
            continue;
          }
        }

        std::vector<unsigned int> target_indices = tracks[target].piece_indices;
        target_indices.push_back(ip);
        Candidate target_with;
        if (!refit_candidate(pieces, target_indices, target_with))
        {
          continue;
        }

        const double ratio = donor_relation_ok
                                 ? target_score / std::max(donor_score, 1.0e-12)
                                 : 0.0;
        if (!found || ratio < best_ratio ||
            (std::fabs(ratio - best_ratio) < 1.0e-12 &&
             target_score < best_target_score))
        {
          found = true;
          best_piece = ip;
          best_donor = donor;
          best_target = target;
          best_ratio = ratio;
          best_target_score = target_score;
        }
      }
    }

    if (!found)
    {
      break;
    }

    std::vector<unsigned int> donor_indices;
    for (unsigned int ip : tracks[best_donor].piece_indices)
    {
      if (ip != best_piece) donor_indices.push_back(ip);
    }
    std::vector<unsigned int> target_indices = tracks[best_target].piece_indices;
    target_indices.push_back(best_piece);

    Candidate donor_refit;
    Candidate target_refit;
    if (!refit_candidate(pieces, donor_indices, donor_refit) ||
        !refit_candidate(pieces, target_indices, target_refit))
    {
      break;
    }

    tracks[best_donor] = std::move(donor_refit);
    tracks[best_target] = std::move(target_refit);
  }
}

Tpc_AssembledTrackReco::SeedParameters
Tpc_AssembledTrackReco::make_seed_parameters(const Candidate& c) const
{
  SeedParameters seed;
  if (!m_idealPadMap)
  {
    return seed;
  }

  const double r_first = m_idealPadMap->get_radius(c.first_layer);
  const double r_last = m_idealPadMap->get_radius(c.last_layer);
  if (!std::isfinite(r_first) || !std::isfinite(r_last))
  {
    return seed;
  }

  const double radius = 0.5 * (r_first + r_last);
  const double phi = wrap_to_pi(predict_phi(c, radius));
  const double dphi_dr = predict_phi_slope(c, radius);
  const double dtbin_dr = c.tbin_slope_r;
  const double tbin = c.tbin_slope_r * radius + c.tbin_intercept_r;
  if (!std::isfinite(radius) || !std::isfinite(phi) ||
      !std::isfinite(dphi_dr) || !std::isfinite(dtbin_dr) || !std::isfinite(tbin))
  {
    return seed;
  }

  seed.x = radius * std::cos(phi);
  seed.y = radius * std::sin(phi);
  seed.z = tbin;

  const double dx_dr = std::cos(phi) - radius * std::sin(phi) * dphi_dr;
  const double dy_dr = std::sin(phi) + radius * std::cos(phi) * dphi_dr;
  const double dz_dr = dtbin_dr;
  const double norm = std::sqrt(dx_dr * dx_dr + dy_dr * dy_dr + dz_dr * dz_dr);
  if (!std::isfinite(norm) || norm <= 0.0)
  {
    return seed;
  }

  seed.px = dx_dr / norm;
  seed.py = dy_dr / norm;
  seed.pz = dz_dr / norm;

  for (auto& i : seed.cov)
  {
    for (double& j : i)
    {
      j = 0.0;
    }
  }
  seed.cov[0][0] = m_seedSigmaX * m_seedSigmaX;
  seed.cov[1][1] = m_seedSigmaY * m_seedSigmaY;
  seed.cov[2][2] = m_seedSigmaZ * m_seedSigmaZ;
  seed.cov[3][3] = m_seedSigmaPx * m_seedSigmaPx;
  seed.cov[4][4] = m_seedSigmaPy * m_seedSigmaPy;
  seed.cov[5][5] = m_seedSigmaPz * m_seedSigmaPz;
  seed.ok = std::isfinite(seed.x) && std::isfinite(seed.y) && std::isfinite(seed.z) &&
            std::isfinite(seed.px) && std::isfinite(seed.py) && std::isfinite(seed.pz);
  return seed;
}

int Tpc_AssembledTrackReco::process_event(PHCompositeNode* /*unused*/)
{
  reset_tree_vars();
  if (m_assembledTrackContainer)
  {
    m_assembledTrackContainer->Reset();
  }

  std::vector<Piece> pieces;
  if (m_tpcModuleTrackContainer)
  {
    const unsigned int n = m_tpcModuleTrackContainer->size();
    pieces.reserve(n);
    for (unsigned int i = 0; i < n; ++i)
    {
      Piece p;
      if (make_piece(i, p))
      {
        pieces.push_back(std::move(p));
      }
    }
  }

  std::vector<std::vector<Candidate> > sector_outputs(24);
  std::vector<std::thread> workers;
  workers.reserve(24);
  for (int side = 0; side < 2; ++side)
  {
    for (unsigned int sector = 0; sector < 12; ++sector)
    {
      const unsigned int index = static_cast<unsigned int>(side) * 12 + sector;
      if (m_enableFastSameSectorPass)
      {
        workers.emplace_back(&Tpc_AssembledTrackReco::connect_sector_pieces_fast,
                             this, std::cref(pieces), side, sector,
                             std::ref(sector_outputs[index]));
      }
      else
      {
        workers.emplace_back(&Tpc_AssembledTrackReco::connect_sector_pieces,
                             this, std::cref(pieces), side, sector,
                             std::ref(sector_outputs[index]));
      }
    }
  }
  for (std::thread& worker : workers)
  {
    worker.join();
  }

  std::vector<Candidate> sector_tracks;
  for (auto& sector_output : sector_outputs)
  {
    sector_tracks.insert(sector_tracks.end(),
                         std::make_move_iterator(sector_output.begin()),
                         std::make_move_iterator(sector_output.end()));
  }

  std::vector<Candidate> assembled_tracks;
  // Short fragments are intentionally left out of the sector pre-pass, so
  // same-sector connections must remain available globally.
  const bool allow_same_sector_global = true;
  connect_side_candidates(pieces, sector_tracks, 0, assembled_tracks,
                          allow_same_sector_global, 1.0,
                          m_enableNormalUniqueTight);
  connect_side_candidates(pieces, sector_tracks, 1, assembled_tracks,
                          allow_same_sector_global, 1.0,
                          m_enableNormalUniqueTight);

  if (m_enableRescuePass && m_rescueIterations > 0 && !assembled_tracks.empty())
  {
    for (unsigned int iter = 0; iter < m_rescueIterations; ++iter)
    {
      const double fraction = static_cast<double>(iter + 1) /
                              static_cast<double>(m_rescueIterations);
      const double scale = 1.0 + fraction * (m_rescueMaxWindowScale - 1.0);

      std::vector<Candidate> rescued;
      rescue_side_candidates(pieces, assembled_tracks, 0, scale, rescued);
      rescue_side_candidates(pieces, assembled_tracks, 1, scale, rescued);
      assembled_tracks.swap(rescued);
    }
  }

  if (m_enableFinalReassociation && !assembled_tracks.empty())
  {
    final_reassociate(pieces, 0, assembled_tracks);
    final_reassociate(pieces, 1, assembled_tracks);
  }

  for (unsigned int it = 0; it < assembled_tracks.size(); ++it)
  {
    const Candidate& c = assembled_tracks[it];

    // Final-track cleanup only:
    // allow every module-track Piece to participate in assembly, then
    // suppress short/trash assembled tracks before they are written to
    // TPC_ASSEMBLEDTRACKS and therefore before CrossingFinder sees them.
    //
    // Candidate::nblobs is the sum of the reconstructed blob/cluster
    // counts of all Pieces belonging to the final assembled track.
    if (m_minClustersPerTrack > 0 &&
        c.nblobs < m_minClustersPerTrack)
    {
      continue;
    }

    const unsigned int assembled_id =
        m_assembledTrackContainer ? m_assembledTrackContainer->size() : it;
    if (m_h_nsegments)
    {
      m_h_nsegments->Fill(static_cast<double>(c.nsegments));
    }
    if (c.nsegments == 3 && m_idealPadMap)
    {
      const double rfirst = m_idealPadMap->get_radius(c.first_layer);
      const double rlast = m_idealPadMap->get_radius(c.last_layer);
      if (std::isfinite(rfirst) && std::isfinite(rlast))
      {
        const double first_tbin = c.tbin_slope_r * rfirst + c.tbin_intercept_r;
        const double last_tbin = c.tbin_slope_r * rlast + c.tbin_intercept_r;
        const double dtbin_track = last_tbin - first_tbin;
        if (m_h_track_tbin_slope_vs_tbin_span_3modules)
        {
          m_h_track_tbin_slope_vs_tbin_span_3modules->Fill(dtbin_track, c.tbin_slope_r);
        }
        if (m_h_track_tbin_slope_vs_first_tbin_3modules)
        {
          m_h_track_tbin_slope_vs_first_tbin_3modules->Fill(first_tbin, c.tbin_slope_r);
        }
        if (m_h_track_tbin_slope_vs_last_tbin_3modules)
        {
          m_h_track_tbin_slope_vs_last_tbin_3modules->Fill(last_tbin, c.tbin_slope_r);
        }
      }
    }

    Tpc_AssembledTrackv1* out = new Tpc_AssembledTrackv1();
    out->set_event(static_cast<unsigned int>(m_event));
    out->set_track_id(assembled_id);
    out->set_side(c.side);
    out->set_nsegments(c.nsegments);
    out->set_nblobs(c.nblobs);
    out->set_nrawhits(c.nrawhits);
    out->set_first_layer(c.first_layer);
    out->set_last_layer(c.last_layer);
    out->set_first_sector(c.first_sector);
    out->set_last_sector(c.last_sector);
    out->set_first_region(c.first_region);
    out->set_last_region(c.last_region);

    const SeedParameters seed = make_seed_parameters(c);
    if (seed.ok)
    {
      out->set_seed_valid(1);
      out->set_seed_x(seed.x);
      out->set_seed_y(seed.y);
      out->set_seed_z(seed.z);
      out->set_seed_px(seed.px);
      out->set_seed_py(seed.py);
      out->set_seed_pz(seed.pz);
      for (unsigned int iseed = 0; iseed < 6; ++iseed)
      {
        for (unsigned int jseed = 0; jseed < 6; ++jseed)
        {
          out->set_seed_cov(iseed, jseed, seed.cov[iseed][jseed]);
        }
      }
    }

    m_tree_track_id.push_back(assembled_id);
    m_tree_side.push_back(c.side);
    m_tree_nsegments.push_back(c.nsegments);
    m_tree_nblobs.push_back(c.nblobs);
    m_tree_nrawhits.push_back(c.nrawhits);
    m_tree_first_layer.push_back(c.first_layer);
    m_tree_last_layer.push_back(c.last_layer);
    m_tree_first_sector.push_back(c.first_sector);
    m_tree_last_sector.push_back(c.last_sector);
    m_tree_first_region.push_back(c.first_region);
    m_tree_last_region.push_back(c.last_region);

    for (unsigned int piece_indice : c.piece_indices)
    {
      const Piece& p = pieces[piece_indice];
      out->add_source_track(p.source_track_id, p.region, p.sector);
      m_tree_source_assembled_track_id.push_back(assembled_id);
      m_tree_source_inmodule_track_id.push_back(p.source_track_id);
      m_tree_source_region.push_back(p.region);
      m_tree_source_sector.push_back(p.sector);
      m_tree_source_side.push_back(p.side);
    }

    // Candidate no longer duplicates all hit keys during every trial refit.
    // Recover them once from the accepted pieces when writing the final track.
    for (unsigned int piece_index : c.piece_indices)
    {
      const Piece& p = pieces[piece_index];

      for (const auto& hit_index : p.hit_indices)
      {
        out->add_hit_index(hit_index.first, hit_index.second);
        m_tree_hit_assembled_track_id.push_back(assembled_id);
        m_tree_hit_hitsetkey.push_back(
            static_cast<unsigned long long>(hit_index.first));
        m_tree_hit_hitkey.push_back(
            static_cast<unsigned long long>(hit_index.second));
      }
    }

    if (m_assembledTrackContainer)
    {
      m_assembledTrackContainer->add_track(out);
    }
    else
    {
      delete out;
    }
  }

  if (m_tree)
  {
    m_tree->Fill();
  }

  if (Verbosity() > 0)
  {
    std::cout << Name() << "::process_event - event " << m_event
              << " input pieces=" << pieces.size()
              << " assembled_tracks=" << assembled_tracks.size() << std::endl;
  }

  ++m_event;
  return Fun4AllReturnCodes::EVENT_OK;
}