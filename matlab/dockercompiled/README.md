# Ironclust Docker

This folder contains a Dockerfile along with instructions on how to build a Docker Image that runs matlab implemented ironclust sorter

1. Compile `p_ironclust.m` script
2. Pack the compiled project into a base image
3. Extend from the base to run ironclust

## Pre steps

Creating a base image to run ironclust


### Matlab Requirements
- JSONLab toolbox
- MATLAB Compiler
- Statistics and Machine Learning Toolbox
- Signal Processing Toolbox
- Image Processing Toolbox

### Steps to compile p\_ironclust and create the base image
- Open Matlab 
- Set Matlab's workspace folder to: `ironclust/matlab`
- Open Matlab's `Application Compiler`
- Click `Open Project`
- Select `p_ironclust.prj` and click `Open`
- Click on `Package` button and wait the packaging process to finish
- Close `Package` and `Application Compiler` windows
- In matlab console run: ```compiler.package.docker('p_ironclust_compiler/for_testing/p_ironclust', 'p_ironclust_compiler/for_testing/requiredMCRProducts.txt', 'ImageName', 'ironclust-base')```

## Building the image

The `Dockerfile` in this folder extends the generated image in previous section (ironclust-base). This is needed in order to properly run ironclust sorter.

To build the image run:

```
docker build -t spikeinterface/ironclust .
```

## Pull the image

Currently the image is pushed in `chyumin` dockerhub profile, it'll be changed to `spikeinterface` in the future

```
docker pull chyumin/ironclust
```

## Running the container

```
docker run -v <host-data-folder>:<docker-data-folder> -it chyumin/ironclust p_ironclust [ARGS]
```

Sample run:
```
docker run -v /opt/data:/opt/data -it chyumin/ironclust p_ironclust /opt/data/raw.mda /opt/data/geom.csv '' '' /opt/data/tmp/firing.mda /opt/data/argfile.txt
``` 

