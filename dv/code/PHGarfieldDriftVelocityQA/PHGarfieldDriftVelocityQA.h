// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef PHGARFIELDCALIBRATION_PHGARFIELDDRIFTVELOCITYQA_H
#define PHGARFIELDCALIBRATION_PHGARFIELDDRIFTVELOCITYQA_H

#include <fun4all/SubsysReco.h>

#include <array>
#include <memory>
#include <unordered_map>
#include <cstdint>
#include <string>

class PHCompositeNode;
class PHGarfield;
class TH1;
class TH2;
class TrkrHitSetContainer;

class PHGarfieldDriftVelocityQA : public SubsysReco
{
 public:
  explicit PHGarfieldDriftVelocityQA(const std::string& name = "PHGarfieldDriftVelocityQA");
  ~PHGarfieldDriftVelocityQA() override;

  int InitRun(PHCompositeNode*) override;
  int process_event(PHCompositeNode*) override;

  void setHitSetNodeName(const std::string& name) { m_hitSetNodeName = name; }
  void setUsePTCorrectedPHGarfield(bool value) { m_usePTCorrectedPHGarfield = value; }
  void setReferenceCMFieldVcm(double value) { m_referenceCMFieldVcm = value; }
  void setReferencePressure(double value) { m_referencePressure = value; }
  void setReferenceTemperatureK(double value) { m_referenceTemperatureK = value; }
  void setTpcTimeBinNs(double value) { m_tpcTimeBinNs = value; }
  void setPHGarfieldT0Bins(double value) { m_phgT0Bins = value; }
  void setPHGarfieldStepNs(double value) { m_phgStepNs = value; }

 private:
  static constexpr int kNFees = 312; // 12 sectors x 26 physical FEEs per side
  bool initializeDriftCoordinateCalibration(PHCompositeNode*);
  bool buildFeeMap(PHCompositeNode*);
  std::string getHistoName(const std::string& base) const;
  double timeBinToZ(int side, int fee, double tbin) const;
  struct PadInfo { int fee; };
  std::unordered_map<std::uint64_t, PadInfo> m_padToFee;
  std::array<std::array<double, kNFees>, 2> m_feeX{}, m_feeY{}, m_feeCount{};
  std::array<std::array<double, kNFees>, 2> m_phgDriftTimeNs{}, m_phgDriftVelocityCmNs{};
  std::array<TH2*, 2> m_hAdcTbinVsFee{}, m_hAdcZVsFee{};
  TH1* m_hNEvents{nullptr};
  TH1* m_hMissingPads{nullptr};
  std::string m_hitSetNodeName{"TRKR_HITSET"};
  bool m_usePTCorrectedPHGarfield{true};

  // Run-79513 reference used in the drift-velocity scan.
  double m_referenceCMFieldVcm{377.0};
  double m_referencePressure{30.0968906};
  double m_referenceTemperatureK{299.98646566};

  // Raw TPC timing convention used in this study.
  double m_tpcTimeBinNs{56.8};
  double m_phgT0Bins{8.0};
  double m_phgStepNs{5.0};
  double m_zPadPlaneCm{102.605};

  std::unique_ptr<PHGarfield> m_driftPHGarfield;
};

#endif
