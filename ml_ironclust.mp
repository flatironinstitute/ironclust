#!/bin/bash

cd `dirname "$0"`

octave -q --eval "ml_ironclust_spec"
#matlab -nodisplay -nosplash -r "ml_ironclust_spec; quit();"