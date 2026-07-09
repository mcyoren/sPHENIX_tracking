#ifndef TPCPOLYTRACKDISPLAY_H
#define TPCPOLYTRACKDISPLAY_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;
class TFile;
class TpcPolyTrackContainer;

class TpcPolyTrackDisplay : public SubsysReco
{
 public:
  TpcPolyTrackDisplay(const std::string& name = "TpcPolyTrackDisplay",
                      const std::string& outfilename = "tpcpolytrack_display.root",
                      const std::string& trackNodeName = "TPCPOLYTRACKS",
                      unsigned int maxEventDisplays = 5);
  ~TpcPolyTrackDisplay() override;

  int Init(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;
  int End(PHCompositeNode*) override;

  void setTrackNodeName(const std::string& n) { m_trackNodeName = n; }
  void setSide(int) {}
  void setZRange(double zmin, double zmax) { m_zmin = zmin; m_zmax = zmax; }
  void setXYRange(double xymax) { m_xymax = xymax; }
  void setLowPtMax(double ptmax) { m_lowPtMax = ptmax; }
  void setHighPtMin(double ptmin) { m_highPtMin = ptmin; }
  void setHighPtMaxAbsD0(double d0max) { m_highPtMaxAbsD0 = d0max; }
  void setMagneticFieldTesla(double b) { m_magneticFieldTesla = b; }

 private:
  bool get_nodes(PHCompositeNode* topNode);

  std::string m_outfilename;
  std::string m_trackNodeName;
  unsigned int m_maxEventDisplays;
  unsigned int m_evt;
  unsigned int m_eventsSaved;
  double m_zmin;
  double m_zmax;
  double m_xymax;
  double m_lowPtMax;
  double m_highPtMin;
  double m_highPtMaxAbsD0;
  double m_magneticFieldTesla;

  TFile* m_outfile;
  TpcPolyTrackContainer* m_tracks;
};

#endif
