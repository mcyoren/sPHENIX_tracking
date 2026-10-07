# Standalone PHGarfieldDriftVelocityQA

This package builds only PHGarfieldDriftVelocityQA.cc into libPHGarfieldDriftVelocityQA.so. It does not build or link the old PHGarfieldCalibrationQA library or its tracking QA modules. Standard installed sPHENIX dependencies, ROOT, PHGarfield, TpcConditions and CDB access are still required.

After sourcing your usual sPHENIX environment, unpack this folder and build out of source:

```bash
cd PHGarfieldDriftVelocityQA
mkdir -p build
cd build
../autogen.sh --prefix="$MYINSTALL"
make -j4
make install
cd ..
```

Use your usual install prefix for MYINSTALL and ensure its include/lib paths are active in your login setup. Header installs to include/phgarfielddriftvelocityqa/PHGarfieldDriftVelocityQA.h; library installs to lib/libPHGarfieldDriftVelocityQA.so. The macro already uses those names. No ROOT dictionary is needed for this non-persistent SubsysReco class.

Run:

```cpp
.x Fun4All_raw_hit_TPC_DriftVelocityQA.C(10,79513,0,".",0,"run3pp","ana532_nocdbtag_v001","DriftVelocityQA")
```

The wrapper keeps your original environment setup and argument order. Change Initialdir in the included Condor file to your work directory. Only the 48 TPC streams are loaded. No silicon/TPOT unpacking or tracking is scheduled. Geometry is read from CDB.

The four large TH2F histograms have 312 x 1000 bins: two per side, about 5 MB total bin contents before any optional ROOT error arrays. Side 0 is South, side 1 North. Names:

- h_PHGarfieldDriftVelocityQA_h_adc_tbin_vs_fee_side0 (and side1)
- h_PHGarfieldDriftVelocityQA_h_adc_z_vs_fee_side0 (and side1)

X is `26 * hardware_sector + physical_FEE`, 0..311 per side. Hardware sector is 0..11 within each endcap, before the unpacker's mc_sectors permutation. Recover sector with `id/26`, physical FEE with `id%26`. The CDB hardware map is inverted to identify the FEE of each unpacked pad; the result is cached. No FEE lookup or PHGarfield integration is performed inside the sample loop. Unmapped samples are counted in h_unmapped_samples and should be zero.

Time axis: 1000 bins in [0,1000); z axis: 1000 bins in [-130,130) cm. Both are weighted by cleaned ADC, using exactly the previous persistent-noise and saturation cleaning. The input ADC is already processed by the unpacker. Underflow/overflow is retained. Antenna pads are excluded by the standard unpacker.

PHGarfield is integrated once per FEE at its mean instrumented pad position. As in the old module, z uses the *average* PHGarfield velocity from the pad plane to CM, rather than interpolating the nonlinear trajectory: z = sign(side) * [z_pad - (tbin-8)*56.8 ns*v]. Increasing time moves toward and beyond the CM; z is never clipped at zero. P/T correction and reference values match your original QA. Failed mapping or any failed FEE drift aborts the run.

The small h_phg_vdrift_side* (cm/us), h_phg_drift_time_side* (ns), and h_calibration_count histograms record calibration. After hadd, divide calibration bin contents by h_calibration_count bin 1 to obtain the segment average. Main ADC histograms sum normally.

Output: QA/tpc_drift_velocity_DriftVelocityQA_<N>evt_<skip>skip_<run>_<segment>_qa.root. Update your merge glob from `tpc_poly_track_RawHitQA_*` to `tpc_drift_velocity_DriftVelocityQA_*`.

To extract side profiles in ROOT:

```cpp
auto* h = dynamic_cast<TH2*>(file->Get("h_PHGarfieldDriftVelocityQA_h_adc_tbin_vs_fee_side0"));
auto* all = h->ProjectionY("time_south",1,312);
int feeID = 42;
auto* fee = h->ProjectionY("time_fee42",feeID+1,feeID+1);
```

The old notebook expects PHGarfieldRawHitsQA names and module bins. It needs these histogram names and the FEE axis substituted before use; it has not been modified in this package.

Validation: source mapping checked against the standard TpcCombinedRawDataUnpacker implementation, all 26 physical FEEs map uniquely onto 26 CDB FEEs, shell syntax checked. ROOT/sPHENIX is unavailable here, so compilation, a short run, event timing, and input memory still require validation in your environment. Start with 10 events and inspect the unmapped sample count and both side projections.
