#!/bin/sh

cd `dirname "$0"`
octave --eval "jrclust_spec; quit;"
