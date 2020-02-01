# IronClust

Terabyte-scale, drift-resistant spike sorter for multi-day recordings from [high-channel-count probes](https://www.nature.com/articles/nature24636)

## Getting Started

## Probe drift handling
IronClust tracks the probe drift by computing the anatomical similarity between time chunks (~20 sec) where each chunk contains approximately equal number of spikes. For each chunk, the anatomical snapshot is computed from the joint distribution bwteen spike amplitudes and positions. Based on the anatomical similarity, each chunk is linked to 15 nearest chunks (self included) forming ~300 sec duration. The linkage is constrained within +/-64 steps (~1280 sec) to handle rigid drift occuring in faster time-scale while rejecting slower changes. The KNN-graph between the spikes is constrained to the linked chunks, such that the neighborhood of spikes from a given chunk is restricted to the spikes from the linked chunks. Thus, drift-resistance is achieved by including and excluding the neighbors based on the activity-inferred anatomical landscape surrounding the probe.

### Prerequisites

- Matlab 
- Matlab signal and image processing toolboxes
- (Optional) CUDA Toolkit (for GPU processing for significant speed-up)
- For terabyte-scale recording: At least 128GB RAM

### Installing
- Clone from Github
```
git clone https://github.com/flatironinstitute/ironclust
```
- (optional) Compile GPU codes (.cu extension)
```
irc2 compile
```

## Quick tutorial

This command creates `irc2` folder in the recording directory and writes output files there.
```
irc2 `path_to_recording_file`
```
Examples 
```
irc2 mydir/myrecoding.bin   # for SpikeGLX format
irc2 mydir/myrecoding.mda   # for .mda format
```

This command writes output files to `output_directory`
```
irc2 `path_to_recording_file` `output_directory`
```

This command shows the parameter file (`.prm` extension) used for sorting
```
irc2 edit `path_to_recording_file`
```

IronClust caches `path_to_recording_file` for subsequent commands. To display the currently selected parameter file, run
```
irc2 which
```

To select a new parameter file, run
```
irc2 select `path_to_prm_file`
```

You can re-run the sorting after updating the sorting parameter by running 
```
irc2 `path_to_recording_file`
```
IronClust only runs the part of sorting pipeline affected by the updated parameters. 

You can initialize the sorting output by running either of the following commands:
```
irc2 clear `path_to_recording_file`
irc2 clear `output_directory`
irc2 clear `path_to_prm_file`
```

## Deployment

- IronClust can run through SpikeForest2 or spikeinterface pipeline
- IronClust output can be exported to Phy, JRClust, Klusters formats for manual clustering

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

To display the current version, run
```
irc2 version
```

## Authors

- James Jun, Center for Computational Mathematics, Flatiron Institute
- Jeremy Magland, Center for Computational Mathematics, Flatiron Institute

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
