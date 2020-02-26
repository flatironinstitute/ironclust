![IronClust logo](img/ironclust_logo.png)

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

### Installation
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
irc2 [path_to_my_recording.mda] (output_dir)  # for .mda format
irc2 [path_to_my_recording.imec#.bin] (output_dir)  # for SpikeGLX Neuropixels recordings (requires `.meta` files)
irc2 [path_to_my_recording.bin] [myprobe.prb] (output_dir) # specify the probe file, output to `myprobe` the recording directory
irc2 [path_to_my_recording.dat] [myprobe.prb] (output_dir)  # for Intan (requires `info.rhd`) and Neuroscope (requires `.xml`) format
```
* `output_dir` (optional): default output location is `irc2` under the recording directory or `myprobe` if the probe file is specified
* `myprobe.prb`: required for Intan and Neuroscope formats. SpikeGLX does not require it if [Neuropixels probe](https://www.neuropixels.org/) is used.

IronClust caches the `path_to_prm_file` for subsequent commands. To display the currently selected parameter file, run
```
irc2 which
```

To select a parameter file (or a recording file):
```
irc2 select [path_to_my_param.prm]
irc2 select [path_to_my_recording]
```

Rerun using new parameters (up to four parameters can be specified, no spaces between name=value pairs):
```
irc2 rerun [path_to_my_param.prm] [name1=val1] [name2=val2] [name3=val3]
irc2 rerun [name1=val1] [name2=val2] [name3=val3] [name4=val4]  # uses a cached parameter file
```

To visualize the raw or filtered traces and see clustered spikes on the traces, run (press 'h' in the UI for further help)
```
irc2 traces [path_to_my_recording] 
irc2 traces [path_to_my_param.prm]
```

Manual clustering user interface
```
irc2 manual [path_to_my_recording] 
irc2 manual [path_to_my_param.prm]
```

This command shows the parameter file (`.prm` extension) used for sorting
```
irc2 edit `path_to_recording_file`
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

## Importing multiple `.bin` files from [SpikeGLX](https://github.com/billkarsh/SpikeGLX)
```
irc2 import-spikeglx [path_to_my_recording.bin] [path_to_probe_file.prb] (path_to_output_dir)
```
- `path_to_output_dir` (optional): defalt location is 'probe_name' under the recording dorectory.
- Output format is [.mda format](https://users.flatironinstitute.org/~magland/docs/mountainsort_dataset_format/) 
- Probe file (`.prb`) is required unless Neuropixels probe is used. [`.prb` file format](https://github.com/JaneliaSciComp/JRCLUST/wiki/Probe-file)
- `path_to_my_recording.bin`: you may use a '\*' character to join multiple files, or provide a text (`.txt`) file containing a list of files to be merged in a specified order (a text file containing the list is created when you use '\*' character). 

## Importing multiple `.dat` files from [Intan RHD format](http://intantech.com/downloads.html?tabSelect=Software&yPos=0)
```
irc2 import-intan [path_to_my_recording.bin] [path_to_probe_file.prb] (path_to_output_dir)
```
- This step is not necessary if all channels are saved to a single file.
- `path_to_my_recording.bin`: Use '\*' character to join all channels that are saved to separate files.

## Deployment

- IronClust can run through SpikeForest2 or spikeinterface pipeline
- IronClust output can be exported to Phy, Klusters, and JRClust formats for manual clustering

## Export to Phy
Export to [Phy](https://github.com/kwikteam/phy-contrib/blob/master/docs/template-gui.md) for manual curation. You need to clone Phy and set the path `path_phy_x` where x={'pc,'mac','lin'} to open the output automatically.
```
irc2 export-phy [path_to_prm_file] (output_dir)   # default output location is `phy` under the output folder
```

If Phy doesn't open automatically, run the following python command to open Phy
```
phy template-gui path_to_param.py
```

## Export to Klusters
Export to [Klusters](http://neurosuite.sourceforge.net/) for manual curation. You can set the path `path_klusters_x` in `user.cfg` where x = {'pc', 'mac', 'lin'} to open the output automatically.
```
irc2 export-klusters [path_to_prm_file] (output_dir)
```
* output_dir (optional): default output location is `klusters` under the same directory.

If Klusters doesn't open automatically, open Klusters GUI and open `.par.#` file (#: shank number). 

## Export to JRCLUST
Export to [JRCLUST](https://github.com/JaneliaSciComp/JRCLUST) for manual curation. You need to clone JRCLUST and set the path `path_jrclust` in `user.cfg` (you need to create this file if it doesn't exist).
```
irc2 export-jrclust [path_to_prm_file]
```
* output_dir: it creates a new JRCLUST parameter file by appending `_jrclust.prm` at the same directory.

If JRCLUST doesn't open automatically, run `jrc manual [my_jrclust.prm]`

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
* We thank [Loren Frank's lab](https://www.cin.ucsf.edu/HTML/Loren_Frank.html) for contributing the terabyte-scale 10-day continuous recording data.
* We thank [Dan English's lab](https://www.englishneurolab.com/) for contributing four-day uLED probe recordings.