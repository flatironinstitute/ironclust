## run `. setup_python.sh (env_name)`
## (env_name) is optional, defaulted to `spikeforest`

CONDA_ENV=${1:-irc}
conda deactivate
conda env remove -n $CONDA_ENV -y || echo "$CONDA_ENV conda environment not removed. Try closing other terminals using $CONDA_ENV"

conda create -ny $CONDA_ENV python=3.6 jupyterlab nodejs ipywidgets
conda activate $CONDA_ENV


## Install mountainsort
git clone https://github.com/flatironinstitute/mountainlab-js sorters/mountainsort
cd sorters/mountainsort
npm install -g .
conda install -c flatiron mountainsort mountainlab -y


## install other sorters from github 
git clone https://github.com/MouseLand/Kilosort2.git sorters/kilosort2
git clone https://github.com/cortex-lab/KiloSort.git sorters/kilosort


## install spyking-circus
#pip install ml_ms4alg
pip install spyking-circus


## install other sorters from github
pip install git+git://github.com/paninski-lab/yass@master
