# TrackingDiagnostics migration to the new TpcTrackReco

## Changes made

- `TpcTrackFit.h`, `TpcTrackHelixFitter.h`, `TpcTrackKalmanFitter.h`, and
  `TpcTrackKalmanFitter.cc` are now owned by `TrackingDiagnostics`.
- `TpcV0CandidateTree` includes the fitter headers locally rather than through
  `tpctrackreco/`.
- The build installs the fitter headers and compiles `TpcTrackKalmanFitter.cc`.
- V0 pattern tracks are matched to `TpcCrossingDecisionContainer` entries using
  `Tpc_PolyTrack::get_source_assembled_track_id()` and
  `TpcCrossingDecision::get_assembled_track_id()`.
- The default crossing node is `TPC_CROSSING_DECISIONS`, matching
  `TpcCrossingFinder` in the new `TpcTrackReco`.
- Pair and track trees now store crossing validity, selected crossing, and the
  source assembled-track ID.
- By default, pairs with two valid but different crossings are rejected. This
  can be disabled with `set_require_same_crossing(false)`.

## Missing uploaded source

The old TpcTrackReco `Makefile.am` references `TpcTrackHelixFitter.cc`, but that
file was not present in the uploaded old-code ZIP. The copied header declares
non-inline functions used by both the V0 module and Kalman fitter, so the
library cannot link until `TpcTrackHelixFitter.cc` is added to this directory
and to `libTrackingDiagnostics_la_SOURCES`.
