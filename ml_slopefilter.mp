#!/bin/sh

cd `dirname "$0"`
octave -q --eval "matlab/ml_slopefilter_spec"
#matlab -nodisplay -nosplash -r "matlab/ml_slopefilter_spec; quit();"
