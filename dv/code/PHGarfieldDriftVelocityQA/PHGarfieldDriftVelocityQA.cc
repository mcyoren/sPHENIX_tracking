#include "PHGarfieldDriftVelocityQA.h"

#include <qautils/QAHistManagerDef.h>

#include <cdbobjects/CDBTTree.h>
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllHistoManager.h>
#include <fun4all/Fun4AllReturnCodes.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <phgarfield/PHGarfield.h>



#include <trackbase/TpcDefs.h>
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrHit.h>
#include <trackbase/TrkrHitSet.h>
#include <trackbase/TrkrHitSetContainer.h>

#include <tpcconditions/TpcConditions.h>

#include <TH1F.h>
#include <TH2F.h>

#include <TPolyLine3D.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace
{

  // ============================================================
  // Raw-channel waveform cleaning
  // ============================================================
  constexpr unsigned int kMaxTimeBins = 2000;
  constexpr unsigned int kNominalMinFrameBins = 700;

  // Persistent channel noise: channel fires through a sizable fraction
  // of the frame (including the beginning of the frame).
  constexpr double kPersistentOccupancy = 0.35;
  constexpr double kPersistentEarlyOccupancy = 0.30;
  constexpr unsigned int kPersistentEarlyBins = 100;
  constexpr unsigned int kMinPersistentHits = 20;

  // A real pulse on top of a noisy/saturated baseline must exceed both
  // an absolute and a relative threshold.
  constexpr double kSignalAboveBaselineADC = 30.0;
  constexpr double kSignalAboveBaselineFraction = 0.50;

  // Saturated SAMPA: after the onset keep only the first N time bins,
  // then suppress the plateau until a new pulse rises above it.
  constexpr unsigned int kKeepBinsAfterFire = 5;
  constexpr double kSaturatedTailOccupancy = 0.35;
  constexpr unsigned int kMinSaturatedTailBins = 20;

  struct TimeAdc
  {
    unsigned int tbin = 0;
    double adc = 0.0;
  };

  double median(std::vector<double> values)
  {
    if (values.empty())
    {
      return 0.0;
    }
    const auto middle = values.begin() + values.size() / 2;
    std::nth_element(values.begin(), middle, values.end());
    double value = *middle;
    if (values.size() % 2 == 0)
    {
      const auto lower = std::max_element(values.begin(), middle);
      value = 0.5 * (value + *lower);
    }
    return value;
  }

  double signal_threshold(const double baseline)
  {
    return std::max(kSignalAboveBaselineADC, kSignalAboveBaselineFraction * baseline);
  }

  bool is_persistent_noise(const std::vector<TimeAdc>& waveform, const unsigned int frame_nbins)
  {
    if (waveform.size() < kMinPersistentHits)
    {
      return false;
    }

    unsigned int early_hits = 0;
    for (const auto& sample : waveform)
    {
      if (sample.tbin < kPersistentEarlyBins)
      {
        ++early_hits;
      }
    }

    const double occupancy = static_cast<double>(waveform.size()) / static_cast<double>(std::max(1U, frame_nbins));
    const double early_occupancy = static_cast<double>(early_hits) / static_cast<double>(kPersistentEarlyBins);
    return occupancy >= kPersistentOccupancy && early_occupancy >= kPersistentEarlyOccupancy;
  }

  double persistent_baseline(const std::vector<TimeAdc>& waveform)
  {
    // Prefer the beginning of the frame, before a physics pulse is likely
    // to dominate. The median is robust against occasional large ADC pulses.
    std::vector<double> values;
    values.reserve(waveform.size());
    for (const auto& sample : waveform)
    {
      if (sample.tbin < kPersistentEarlyBins)
      {
        values.push_back(sample.adc);
      }
    }
    if (values.size() < 5)
    {
      values.clear();
      for (const auto& sample : waveform)
      {
        values.push_back(sample.adc);
      }
    }
    return median(std::move(values));
  }

  int find_saturation_onset(const std::vector<TimeAdc>& waveform, const unsigned int frame_nbins)
  {
    // Find the earliest hit for which the remainder of the frame has a
    // persistent occupancy. A normal short pulse therefore does not qualify.
    if (waveform.size() < kMinSaturatedTailBins)
    {
      return -1;
    }

    for (std::size_t i = 0; i < waveform.size(); ++i)
    {
      const unsigned int onset = waveform[i].tbin;
      const unsigned int tail_start = onset + kKeepBinsAfterFire;
      if (tail_start >= frame_nbins)
      {
        break;
      }

      const auto first_tail = std::lower_bound(
          waveform.begin(), waveform.end(), tail_start,
          [](const TimeAdc& sample, const unsigned int tbin)
          { return sample.tbin < tbin; });

      const unsigned int tail_hits = static_cast<unsigned int>(waveform.end() - first_tail);
      const unsigned int tail_bins = frame_nbins - tail_start;
      if (tail_hits >= kMinSaturatedTailBins &&
          static_cast<double>(tail_hits) / static_cast<double>(tail_bins) >= kSaturatedTailOccupancy)
      {
        // Do not call a channel saturated if it was already persistently noisy
        // from the beginning of the frame; that is handled separately.
        const unsigned int before_bins = std::max(1U, onset);
        const unsigned int before_hits = static_cast<unsigned int>(i);
        const double before_occupancy =
            static_cast<double>(before_hits) / static_cast<double>(before_bins);
        if (before_occupancy < kPersistentEarlyOccupancy)
        {
          return static_cast<int>(onset);
        }
      }
    }
    return -1;
  }

  std::vector<TimeAdc> clean_waveform(std::vector<TimeAdc> waveform, const unsigned int frame_nbins)
  {
    if (waveform.empty())
    {
      return {};
    }

    std::sort(
        waveform.begin(), waveform.end(),
        [](const TimeAdc& lhs, const TimeAdc& rhs)
        { return lhs.tbin < rhs.tbin; });

    std::vector<TimeAdc> cleaned;
    cleaned.reserve(waveform.size());

    // ----------------------------------------------------------
    // Case 1: channel fires throughout the frame.
    // Remove its robust per-channel baseline and keep only a pulse
    // that rises significantly above that baseline.
    // ----------------------------------------------------------
    if (is_persistent_noise(waveform, frame_nbins))
    {
      const double baseline = persistent_baseline(waveform);
      const double threshold = signal_threshold(baseline);

      for (const auto& sample : waveform)
      {
        const double excess = sample.adc - baseline;
        if (excess > threshold)
        {
          cleaned.push_back({sample.tbin, excess});
        }
      }
      return cleaned;
    }

    // ----------------------------------------------------------
    // Case 2: quiet channel fires and then stays on (SAMPA saturation).
    // ----------------------------------------------------------
    const int saturation_onset = find_saturation_onset(waveform, frame_nbins);
    if (saturation_onset < 0)
    {
      return waveform;
    }

    const unsigned int onset = static_cast<unsigned int>(saturation_onset);
    const unsigned int first_keep_end = std::min(frame_nbins, onset + kKeepBinsAfterFire);

    // Keep the original first five bins after the onset.
    for (const auto& sample : waveform)
    {
      if (sample.tbin >= onset && sample.tbin < first_keep_end)
      {
        cleaned.push_back(sample);
      }
    }

    // Estimate the saturated plateau after the first kept bins.
    std::vector<double> plateau_values;
    for (const auto& sample : waveform)
    {
      if (sample.tbin >= first_keep_end)
      {
        plateau_values.push_back(sample.adc);
      }
    }
    double plateau = median(std::move(plateau_values));
    double threshold = signal_threshold(plateau);

    // While suppressed, a new pulse above the plateau re-opens a 5-bin
    // acceptance window. Its ADC is baseline-subtracted.
    unsigned int accept_until = first_keep_end;
    bool in_retrigger_window = false;

    for (const auto& sample : waveform)
    {
      if (sample.tbin < first_keep_end)
      {
        continue;
      }

      if (sample.tbin < accept_until && in_retrigger_window)
      {
        const double excess = sample.adc - plateau;
        if (excess > 0.0)
        {
          cleaned.push_back({sample.tbin, excess});
        }
        continue;
      }

      const double excess = sample.adc - plateau;
      if (excess > threshold)
      {
        in_retrigger_window = true;
        accept_until = std::min(frame_nbins, sample.tbin + kKeepBinsAfterFire);
        cleaned.push_back({sample.tbin, excess});
      }
      else
      {
        in_retrigger_window = false;
      }
    }

    return cleaned;
  }

  double temperature_kelvin(double t) { return t < 200.0 ? t + 273.15 : t; }
  bool is_good_number(double v) { return std::isfinite(v) && std::fabs(v) < 1.e30; }
}

PHGarfieldDriftVelocityQA::PHGarfieldDriftVelocityQA(const std::string& name): SubsysReco(name) {}
PHGarfieldDriftVelocityQA::~PHGarfieldDriftVelocityQA() = default;
std::string PHGarfieldDriftVelocityQA::getHistoName(const std::string& base) const { return "h_" + Name() + "_" + base; }

bool PHGarfieldDriftVelocityQA::buildFeeMap(PHCompositeNode* topNode)
{
  auto* geom = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  const std::string url = CDBInterface::instance()->getUrl("TPC_FEE_CHANNEL_MAP");
  if (!geom || url.empty()) return false;
  CDBTTree map(url); map.LoadCalibrations();
  // Same physical FEE -> CDB mapping and sector convention as the TPC unpacker.
  const int feeMap[26]{4,5,0,2,1,11,9,10,8,7,6,0,1,3,7,6,5,4,3,2,0,2,1,3,5,4};
  const int feeR[26]{2,2,1,1,1,3,3,3,3,3,3,2,2,1,2,2,1,1,2,2,3,3,3,3,3,3};
  const int sectors[12]{5,4,3,2,1,0,11,10,9,8,7,6};
  m_padToFee.clear(); m_feeX = {}; m_feeY = {}; m_feeCount = {};
  for (int side=0; side<2; ++side) for (int sec=0; sec<12; ++sec)
    for (int fee=0; fee<26; ++fee) for (int ch=0; ch<256; ++ch)
    {
      const int mapped = feeMap[fee] + (feeR[fee]==2 ? 6 : feeR[fee]==3 ? 14 : 0);
      const int key = 256*mapped + ch;
      const int layer = map.GetIntValue(key,"layer");
      if (layer<7 || layer>54) continue; // antenna pads are not in TRKR_HITSET
      auto* lg = geom->GetLayerCellGeom(layer); if (!lg) return false;
      const double phi = (side==1 ? 1. : -1.)*(map.GetDoubleValue(key,"phi")-M_PI/2.) + sec*M_PI/6.;
      if (!is_good_number(phi)) return false;
      const unsigned int pad = lg->get_phibin(phi,side);
      const auto hs = TpcDefs::genHitSetKey(layer,sectors[sec],side);
      const std::uint64_t lookup = (std::uint64_t(hs)<<16) | pad;
      const int id=26*sec+fee;
      auto result=m_padToFee.emplace(lookup,PadInfo{id});
      if (!result.second && result.first->second.fee != id) return false;
      if (!result.second) continue;
      m_feeX[side][id] += lg->get_radius()*std::cos(phi);
      m_feeY[side][id] += lg->get_radius()*std::sin(phi);
      m_feeCount[side][id] += 1.;
    }
  for (int s=0;s<2;++s) for(int f=0;f<kNFees;++f)
  {
    if (m_feeCount[s][f]<=0.) return false;
    m_feeX[s][f]/=m_feeCount[s][f]; m_feeY[s][f]/=m_feeCount[s][f];
  }
  return true;
}
bool PHGarfieldDriftVelocityQA::initializeDriftCoordinateCalibration(PHCompositeNode* topNode)
{
  for (auto& side_values : m_phgDriftTimeNs)
  {
    side_values.fill(0.0);
  }
  for (auto& side_values : m_phgDriftVelocityCmNs)
  {
    side_values.fill(0.0);
  }

  auto* conditions = findNode::getClass<TpcConditions>(topNode, "TpcConditions");
  if (!conditions || !conditions->get_ConditionsAvailable())
  {
    std::cout << PHWHERE << " " << Name()
              << " TpcConditions unavailable for PHGarfield P/T calibration" << std::endl;
    return false;
  }

  auto* geom_container = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  if (!geom_container)
  {
    std::cout << PHWHERE << " " << Name() << " missing TPCGEOMCONTAINER" << std::endl;
    return false;
  }

  auto* layer_geom = geom_container->GetLayerCellGeom(20);
  if (!layer_geom)
  {
    std::cout << PHWHERE << " " << Name() << " missing TPC geometry for layer 20" << std::endl;
    return false;
  }

  m_zPadPlaneCm = layer_geom->get_max_driftlength() + layer_geom->get_CM_halfwidth();
  if (!is_good_number(m_zPadPlaneCm) || m_zPadPlaneCm <= 0.0)
  {
    std::cout << PHWHERE << " " << Name() << " invalid pad-plane z = " << m_zPadPlaneCm << std::endl;
    return false;
  }

  const double pressure = conditions->get_Pressure();
  const double temperatureK = temperature_kelvin(conditions->get_Temperature());

  double effective_field_vcm = m_referenceCMFieldVcm;
  if (m_usePTCorrectedPHGarfield)
  {
    if (!(pressure > 0.0 && temperatureK > 0.0 &&
          m_referencePressure > 0.0 && m_referenceTemperatureK > 0.0))
    {
      std::cout << PHWHERE << " " << Name()
                << " invalid P/T for PHGarfield scaling:"
                << " P=" << pressure << " T=" << temperatureK << std::endl;
      return false;
    }

    // Match the run-dependent reduced field E/N with the reference gas table:
    // E_eff = E_ref * (P_ref/P_run) * (T_run/T_ref).
    effective_field_vcm = m_referenceCMFieldVcm *
                          (m_referencePressure / pressure) *
                          (temperatureK / m_referenceTemperatureK);
  }

  std::cout << Name() << "::initializeDriftCoordinateCalibration"
            << " P=" << pressure
            << " T=" << temperatureK << " K"
            << " Eref=" << m_referenceCMFieldVcm << " V/cm"
            << " Eeff=" << effective_field_vcm << " V/cm"
            << " zPad=" << m_zPadPlaneCm << " cm"
            << " t0=" << m_phgT0Bins << " bins"
            << " dt=" << m_tpcTimeBinNs << " ns/bin"
            << std::endl;

  m_driftPHGarfield = std::make_unique<PHGarfield>(Name() + "_DriftPHGarfield");
  m_driftPHGarfield->SetCMVoltageDefault(effective_field_vcm);

  const int garfield_ret = m_driftPHGarfield->InitRun(topNode);
  if (garfield_ret != Fun4AllReturnCodes::EVENT_OK)
  {
    std::cout << PHWHERE << " " << Name()
              << " PHGarfield::InitRun failed with code " << garfield_ret << std::endl;
    m_driftPHGarfield.reset();
    return false;
  }

  for (int side=0;side<2;++side) for(int fee=0;fee<kNFees;++fee)
  {
    PHGarfield::ReverseDriftStatus status = PHGarfield::ReverseDriftStatus::Running;
    std::unique_ptr<TPolyLine3D> path(m_driftPHGarfield->ReverseDrift(
      m_feeX[side][fee],m_feeY[side][fee],side==0 ? -m_zPadPlaneCm : m_zPadPlaneCm,m_phgStepNs,&status));
    if (!path || path->GetN()<1 || status!=PHGarfield::ReverseDriftStatus::CentralMembrane) return false;
    const double time = path->GetN()*m_phgStepNs;
    m_phgDriftTimeNs[side][fee]=time;
    m_phgDriftVelocityCmNs[side][fee]=m_zPadPlaneCm/time;
  }
  m_driftPHGarfield.reset(); // integration only at InitRun, never per hit
  return true;
}
double PHGarfieldDriftVelocityQA::timeBinToZ(const int side, const int moduleIndex, const double tbin) const
{
  if (side < 0 || side > 1 || moduleIndex < 0 || moduleIndex >= kNFees)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double drift_velocity_cm_ns = m_phgDriftVelocityCmNs[side][moduleIndex];
  if (!(drift_velocity_cm_ns > 0.0) || !std::isfinite(drift_velocity_cm_ns))
  {
    return std::numeric_limits<double>::quiet_NaN();
  }

  // tbin = 8 maps to the pad plane.  Increasing raw time moves toward the CM.
  // Do not clip at z=0: an observed CM edge displaced from zero is exactly the
  // diagnostic we want when testing the PHGarfield drift velocity.
  const double drift_time_ns = (tbin - m_phgT0Bins) * m_tpcTimeBinNs;
  const double drift_distance_cm = drift_time_ns * drift_velocity_cm_ns;
  const double abs_z = m_zPadPlaneCm - drift_distance_cm;
  return side == 0 ? -abs_z : abs_z;
}


int PHGarfieldDriftVelocityQA::InitRun(PHCompositeNode* topNode)
{
  if (!(m_tpcTimeBinNs>0.) || !(m_phgStepNs>0.) || !std::isfinite(m_phgT0Bins) ||
      !buildFeeMap(topNode) || !initializeDriftCoordinateCalibration(topNode))
  { std::cerr << Name() << ": FEE mapping/drift calibration failed" << std::endl; return Fun4AllReturnCodes::ABORTRUN; }
  auto* hm=QAHistManagerDef::getHistoManager();
  m_hNEvents=new TH1F(getHistoName("h_nEvents").c_str(),"Processed events;;Events",1,.5,1.5);
  m_hMissingPads=new TH1F(getHistoName("h_unmapped_samples").c_str(),"Unmapped samples;;Samples",1,.5,1.5);
  hm->registerHisto(m_hNEvents); hm->registerHisto(m_hMissingPads);
  for(int s=0;s<2;++s)
  {
    const std::string suffix="_side"+std::to_string(s);
    m_hAdcTbinVsFee[s]=new TH2F(getHistoName("h_adc_tbin_vs_fee"+suffix).c_str(),
       "ADC weighted raw time;FEE ID (26*hardware sector + physical FEE);Raw time bin;Sum ADC",312,-.5,311.5,1000,0.,1000.);
    m_hAdcZVsFee[s]=new TH2F(getHistoName("h_adc_z_vs_fee"+suffix).c_str(),
       "ADC weighted PHGarfield z;FEE ID (26*hardware sector + physical FEE);PHGarfield z [cm];Sum ADC",312,-.5,311.5,1000,-130.,130.);
    hm->registerHisto(m_hAdcTbinVsFee[s]); hm->registerHisto(m_hAdcZVsFee[s]);
    auto* dv=new TH1F(getHistoName("h_phg_vdrift"+suffix).c_str(),"PHGarfield average drift velocity;FEE ID;cm/us",312,-.5,311.5);
    auto* dt=new TH1F(getHistoName("h_phg_drift_time"+suffix).c_str(),"PHGarfield drift time;FEE ID;ns",312,-.5,311.5);
    for(int f=0;f<kNFees;++f) { dv->SetBinContent(f+1,1000.*m_phgDriftVelocityCmNs[s][f]); dt->SetBinContent(f+1,m_phgDriftTimeNs[s][f]); }
    hm->registerHisto(dv); hm->registerHisto(dt);
  }
  auto* count=new TH1F(getHistoName("h_calibration_count").c_str(),"Number of segment calibrations;;Segments",1,.5,1.5);
  count->SetBinContent(1,1.); hm->registerHisto(count);
  return Fun4AllReturnCodes::EVENT_OK;
}
int PHGarfieldDriftVelocityQA::process_event(PHCompositeNode* topNode)
{
  auto* hits=findNode::getClass<TrkrHitSetContainer>(topNode,m_hitSetNodeName);
  if (!hits) return Fun4AllReturnCodes::ABORTRUN;
  m_hNEvents->Fill(1.);
  unsigned int frame=kNominalMinFrameBins;
  const auto range=hits->getHitSets(TrkrDefs::TrkrId::tpcId);
  for(auto it=range.first;it!=range.second;++it)
  {
    if(!it->second) continue;
    const auto hr=it->second->getHits();
    for(auto j=hr.first;j!=hr.second;++j) { const auto t=TpcDefs::getTBin(j->first); if(t<kMaxTimeBins) frame=std::max(frame,t+1U); }
  }
  for(auto it=range.first;it!=range.second;++it)
  {
    if(!it->second) continue;
    const int side=TpcDefs::getSide(it->first); if(side<0 || side>1) continue;
    std::map<unsigned int,std::vector<TimeAdc>> waves;
    const auto hr=it->second->getHits();
    for(auto j=hr.first;j!=hr.second;++j)
    { const auto t=TpcDefs::getTBin(j->first); if(j->second && t<kMaxTimeBins) waves[TpcDefs::getPad(j->first)].push_back({t,double(j->second->getAdc())}); }
    for(auto& channel:waves)
    {
      const auto found=m_padToFee.find((std::uint64_t(it->first)<<16)|channel.first);
      if(found==m_padToFee.end()) { m_hMissingPads->Fill(1.,channel.second.size()); continue; }
      const int fee=found->second.fee;
      for(const auto& sample:clean_waveform(std::move(channel.second),frame))
      {
        if (!(sample.adc>0.)) continue;
        // Retain out-of-range samples in ROOT under/overflow bins.
        m_hAdcTbinVsFee[side]->Fill(fee,sample.tbin,sample.adc);
        m_hAdcZVsFee[side]->Fill(fee,timeBinToZ(side,fee,sample.tbin),sample.adc);
      }
    }
  }
  return Fun4AllReturnCodes::EVENT_OK;
}
