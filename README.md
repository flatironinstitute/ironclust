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
irc2 [path_to_my_recording.imec#.bin]   # for SpikeGLX Neuropixels recordings (.bin and .meta files)
irc2 [path_to_my_recording.mda]   # for .mda format
irc2 [path_to_my_recording] (output_dir)  # specify output directory (default location is `irc2` under recording directory)
```

You can import [SpikeGLX](https://github.com/billkarsh/SpikeGLX) format to [MDA format](https://users.flatironinstitute.org/~magland/docs/mountainsort_dataset_format/) by supplying a `.bin` file and [`.prb` (probe) file](https://github.com/JaneliaSciComp/JRCLUST/wiki/Probe-file). Make sure that `.meta` file exists in the same directory. Multiple `.bin` files can be joined if you provide a wild card for `[path_to_my_recording.bin]` or supply a `.txt` (text) file containing   a list of files to be merged.
```
irc2 import-spikeglx [path_to_my_recording.bin] [path_to_probe_file.prb] (path_to_output_dir)
```

Export to [Phy](https://github.com/kwikteam/phy-contrib/blob/master/docs/template-gui.md) for manual curation
```
irc2 export-phy [path_to_prm_file] (output_dir)   # default output location is `phy` under the output folder
```

and run the following python command to open Phy
```
phy template-gui path_to_param.py
```

This command shows the parameter file (`.prm` extension) used for sorting
```
irc2 edit `path_to_recording_file`
```

IronClust caches the `path_to_prm_file` for subsequent commands. To display the currently selected parameter file, run
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

* We thank our collaborators and contributors of the ground-truth datasets to validate our spike sorting accuracy through spikeforest.flatironinstitute.org website.
* We thank Loren Frank's lab for contributing the terabyte-scale 10-day continuous recording data.
* We thank Dan English's lab for contributing four-day uLED probe recordings.