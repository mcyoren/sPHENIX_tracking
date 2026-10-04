#!/bin/bash

#runs=(79507 79508 79509 79510 79511 79512 76905 82626 75405 75391)
#runs=(76905  75405 75391)
#runs=(82626)
#runs=(79507 79508 79509 79510 79511 79512)
#runs=(79508 79509 79510 79511 79512 79513 79514 79515 79516 79522 79523 79524 79525 79526 79527 79528 79529 79530 81558 76905 67785)
 runs=(80892 80909 81566)

for run in "${runs[@]}"; do
    hadd -f -k -O \
        /sphenix/user/mitrankov/garf/input/QA/QA2_${run}.root \
        /sphenix/tg/tg01/hf/mitrankov/ibf/output/QA/tpc_poly_track_RawHitQA_10evt_0skip_${run}* &
done
