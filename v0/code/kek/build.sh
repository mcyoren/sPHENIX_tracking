#!/usr/bin/env bash

cd ../tpctrackreco/build/

make -j8 install 

cd ../

cd ../offline/packages/TrackingDiagnostics/build/

make -j8 install 

cd ../../../../work
