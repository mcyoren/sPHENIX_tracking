Hi Ben,

Thanks again for the macro and the PR! I have it running on 79516 and it works nicely.

We found one small thing that should be an easy fix. With `set_useGlobal(true)`, the cluster positions are moved into global coordinates, which include the TPC survey shift and tilt. But the truth CM pattern stays in the ideal TPC frame. So the survey shift of the TPC (about 1.7 mm in x) shows up as a fake distortion: mostly a sin/cos pattern in phi, and bigger on the South side (around 3 mm with Garfield on). The laser can't really see a shift of the whole TPC, because the pads and the membrane move together.

Two simple ways to fix it:
1. Run the matching with `set_useGlobal(false)`, so the clusters and the truth are both in TPC coordinates. Or:
2. Move the truth positions to global with `transformTpcEnvelopeToWorld` before matching.

With either one, the survey and ideal geometry give the same maps. I checked this with the no-Garfield maps.

Two tiny things I noticed on the way:
- In `Fun4All_LaserAggregatedClustering_simple.C`, the clusterizer is registered as `clus2`, but the variable is `clus`.
- In `LaserClusterHelper`, the hardcoded PHGarfield TPC rotation uses an x-tilt of +0.000298, while the G4 survey uses -0.00029817. Calling `SetUseSurveyGeometry(true)` instead of the hardcoded `MoveTpc`/`RotateTpc` would keep PHGarfield and Acts on the same numbers. This one has almost no effect on the maps.

I also have a small patch on top of your branch that lets the macro choose the PHGarfield field map, keff and TPC geometry from a config file. Happy to share it if it's useful.

Thanks!
Yuri
