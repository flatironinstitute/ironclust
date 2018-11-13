#!/bin/bash

cd `dirname "$0"`

octave -q --eval "matlab/ml_ironclust_spec"
#matlab -nodisplay -nosplash -r "matlab/ml_ironclust_spec; quit();"