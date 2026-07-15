# TpcDstV0Finder

Standalone Fun4All analysis module that builds K0S, Lambda, and anti-Lambda candidates directly from the fitted states stored in `FINALTRACKS`.

## Important interpretation

`FinalTrack::{x,y,z,px,py,pz,charge}` is treated as one complete fitted track state at a common reference point on the trajectory. No cluster fit is repeated. In a uniform solenoidal field the module analytically propagates this state as a helix, finds the 3D closest approach of two helices, and evaluates each daughter momentum at its own PCA point.

The candidate vertex is the midpoint of the two PCA points. Invariant masses, Armenteros variables, DIRA, and V0 momentum are all calculated with the propagated daughter momenta at that candidate vertex.

## Inputs

Required:
- `FINALTRACKS` (`FinalTrackContainer`)

Optional:
- `FINALTRACKVERTICES` (`FinalTrackVertexContainer`), collision vertex index 0

There is no truth-node dependency and no dependency on `TpcPolyClusterTrackContainer`.

## Build

```bash
source /opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL=/path/to/your/install
export OFFLINE_MAIN=/path/to/the/install/that/contains/inmoduletracks

./autogen.sh
./configure --prefix="$MYINSTALL"
make -j8
make install
```

Make sure both `$MYINSTALL/lib` and the reconstruction install containing `libinmoduletracks.so` are in `LD_LIBRARY_PATH`.

## Run

```bash
root -l -b -q 'macro/Fun4All_TpcDstV0.C("input.dst","TpcDstV0.root",0)'
```

## Trees

- `pairTree`: pair PCA, daughter momenta at PCA, K0S/Lambda/anti-Lambda masses, Armenteros variables, DIRA, decay radius, track quality, stored reference states.
- `trackTree`: optional QA tree containing exactly the stored final-track state plus derived pT, eta, and DCA to the chosen primary vertex.

`refit_used` is currently always zero. `enable_vertex_constrained_refit(true)` only prints a warning and is reserved for a later controlled extension.

## First QA checks

```cpp
pairTree->Draw("mass_Kshort>>h(200,0.3,0.7)",
               "charge1*charge2<0&&pairDCA<3&&decay_radius>4&&abs(pca_z)<10&&cosThetaReco>0.7");

pairTree->Draw("qT:alpha>>hap(160,-1,1,160,0,0.4)",
               "charge1*charge2<0&&pairDCA<3", "COLZ");

trackTree->Draw("dedx:charge*pt>>hdedx(200,-5,5,200,0,2000)",
                "nclusters>30", "COLZ");
```

## Sign convention check

The propagation uses the same signed-radius convention as `FinalTrackVertexer`: `R_signed = pt/(0.003*q*B)`, with circle center `xc=x+R_signed*py/pt`, `yc=y-R_signed*px/pt`. Run the single-track QA and compare propagated points against your display for one positive and one negative track. If your reconstruction defines the stored charge/B convention oppositely, change only the sign of `omega` in `make_state`; do not refit the tracks.
