#!/bin/bash

python ../transform_traces.py -a mito10_backchecked.amod -s G3-SA_modeled_cropped_1464_alignem_project.json -r rev_transforms.json -y 29 -o out_mito10_backchecked.amod
EC=$?
if [[ "$EC" != "0" ]]; then
    echo "Test failed: transform_traces.py returned non-zero exit code."
    exit 1
fi

diff out_mito10_backchecked.amod ref_m28_mito10_backchecked.amod
EC=$?
if [[ "$EC" != "0" ]]; then
    echo "Test failed: output file differs from reference."
    exit 1
else
    echo "Ok, test passed" 
fi
