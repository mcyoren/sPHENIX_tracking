#ifndef TPCPOLYCLUSTERRESIDUALS_H
#define TPCPOLYCLUSTERRESIDUALS_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <string>
#include <vector>

class FinalTrackContainer;
class FinalTrackVertexContainer;
class PHCompositeNode;
class TFile;
class TTree;
class TpcPolyClusterTrackContainer;

class TpcPolyClusterResiduals : public SubsysReco
{
 public:
  explicit TpcPolyClusterResiduals(const std::string& name = "TpcPolyClusterResiduals",
                                   const std::string& outfilename = "tpcpolycluster_residuals.root");
  ~TpcPolyClusterResiduals() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setClusterNodeName(const std::string& n) { m_clusterNodeName = n; }
  void setFinalTrackNodeName(const std::string& n) { m_finalTrackNodeName = n; }
  void setFinalTrackVertexNodeName(const std::string& n) { m_finalTrackVertexNodeName = n; }
  void setOutputFileName(const std::string& n) { m_outfilename = n; }
  void setMagneticFieldTesla(double b) { m_magneticFieldTesla = b; }
  void setMinPt(double v) { m_minPt = v; }
  void setMaxPt(double v) { m_maxPt = v; }
  void setMinTpcClusters(unsigned int v) { m_minTpcClusters = v; }
  void setMaxTpcClusters(unsigned int v) { m_maxTpcClusters = v; }
  void setDedxThicknessPerRegion(double inner_odd, double inner_even, double mid, double outer)
  {
    m_dedxThicknessPerRegion = {{inner_odd, inner_even, mid, outer}};
  }

 private:
  bool get_nodes(PHCompositeNode* topNode);
  void reset_tree_values();

  std::string m_outfilename;
  std::string m_clusterNodeName;
  std::string m_finalTrackNodeName;
  std::string m_finalTrackVertexNodeName;

  double m_magneticFieldTesla {1.4};
  double m_minPt {0.0};
  double m_maxPt {1.0e30};
  unsigned int m_minTpcClusters {0};
  unsigned int m_maxTpcClusters {0xffffffffu};
  std::array<double, 4> m_dedxThicknessPerRegion {{1.0, 1.0, 1.0, 1.0}};

  unsigned int m_evt {0};
  TFile* m_outfile {nullptr};
  TTree* m_tree {nullptr};
  TpcPolyClusterTrackContainer* m_clusterTracks {nullptr};
  FinalTrackContainer* m_finalTracks {nullptr};
  FinalTrackVertexContainer* m_finalTrackVertices {nullptr};

  unsigned int m_event {0};
  unsigned int m_finalTrackId {0};
  unsigned int m_clusterTrackId {0};
  unsigned int m_sourceFullTrackId {0};
  int m_side {0};
  unsigned int m_ntpcClusters {0};
  int m_fitStatus {0};
  double m_pt {0.0};
  double m_px {0.0};
  double m_py {0.0};
  double m_pz {0.0};
  double m_eta {0.0};
  double m_theta {0.0};
  double m_charge {0.0};
  double m_chi2 {0.0};
  double m_ndf {0.0};
  double m_quality {0.0};
  double m_dedx {0.0};
  double m_vertexX {0.0};
  double m_vertexY {0.0};
  double m_vertexZ {0.0};
  double m_vertexR {0.0};
  double m_pcaX {0.0};
  double m_pcaY {0.0};
  double m_pcaZ {0.0};
  double m_zDCA {0.0};
  double m_rDCA {0.0};
  double m_rDCAZero {0.0};
  double m_R {0.0};
  double m_rzSlope {0.0};
  std::vector<unsigned int> m_clusterIndex;
  std::vector<unsigned int> m_sector;
  std::vector<unsigned int> m_layer;
  std::vector<double> m_clusterX;
  std::vector<double> m_clusterY;
  std::vector<double> m_clusterZ;
  std::vector<double> m_clusterR;
  std::vector<double> m_clusterPhi;
  std::vector<double> m_clusterAdc;
  std::vector<unsigned int> m_clusterPadSize;
  std::vector<double> m_stateX;
  std::vector<double> m_stateY;
  std::vector<double> m_stateZ;
  std::vector<double> m_stateZDca;
  std::vector<double> m_stateR;
  std::vector<double> m_statePhi;
  std::vector<double> m_deltaPhi;
  std::vector<double> m_residualRPhi;
  std::vector<double> m_residualZ;
};

#endif
