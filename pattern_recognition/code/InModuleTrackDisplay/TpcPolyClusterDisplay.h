#ifndef TPCPOLYCLUSTERDISPLAY_H
#define TPCPOLYCLUSTERDISPLAY_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class FinalTrackContainer;
class FinalTrackVertexContainer;
class TFile;
class TpcPolyClusterTrackContainer;

class TpcPolyClusterDisplay : public SubsysReco
{
 public:
  TpcPolyClusterDisplay(const std::string& name = "TpcPolyClusterDisplay",
                        const std::string& outfilename = "tpcpolycluster_display.root",
                        const std::string& clusterNodeName = "TPCPOLYCLUSTERTRACKS",
                        unsigned int maxEventDisplays = 5);
  ~TpcPolyClusterDisplay() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setClusterNodeName(const std::string& n) { m_clusterNodeName = n; }
  void setFinalTrackNodeName(const std::string& n) { m_finalTrackNodeName = n; }
  void setFinalTrackVertexNodeName(const std::string& n) { m_finalTrackVertexNodeName = n; }
  void setMagneticFieldTesla(double b) { m_magneticFieldTesla = b; }
  void setZRange(double zmin, double zmax) { m_zmin = zmin; m_zmax = zmax; }
  void setTrackVertexZRange(double zmin, double zmax) { m_trackVertexZMin = zmin; m_trackVertexZMax = zmax; }
  void setXYRange(double xymax) { m_xymax = xymax; }

 private:
  bool get_nodes(PHCompositeNode* topNode);

  std::string m_outfilename;
  std::string m_clusterNodeName;
  std::string m_finalTrackNodeName;
  std::string m_finalTrackVertexNodeName;
  unsigned int m_maxEventDisplays;
  unsigned int m_evt;
  unsigned int m_eventsSaved;
  double m_zmin;
  double m_zmax;
  double m_trackVertexZMin;
  double m_trackVertexZMax;
  double m_xymax;
  double m_magneticFieldTesla;

  TFile* m_outfile;
  TpcPolyClusterTrackContainer* m_clusterTracks;
  FinalTrackContainer* m_finalTracks;
  FinalTrackVertexContainer* m_finalTrackVertices;
};

#endif
