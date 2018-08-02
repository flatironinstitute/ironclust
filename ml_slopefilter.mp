#!/bin/sh

cd `dirname "$0"`
octave -q --eval "ml_slopefilter_spec"
#matlab -nodisplay -nosplash -r "ml_slopefilter_spec; quit();"
