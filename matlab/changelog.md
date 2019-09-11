# Changelog
IronClust, written by J. James Jun, Flatiron Institute, Simons Foundation
The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [4.9.11] - 2019-9-1
### Added
- Peak memory usage is reported in `irc describe` command
  -Unit is in gigabytes (GiB)

### Changed
- `run_irc` will use `argfile.txt` instead of `params.json` if present


## [4.9.10] - 2019-9-10
### Changed
- Comment changed in `default.prm` regarding `qqFactor` parameter
- Default value changaed for `post_merge_mode0 = [12, 15, 17]`
  - Previously was [15, 19, 17]


## [4.9.9] - 2019-9-10
### Added
- New pre-merging mode based on KNN and drift (`post_merge_mode0 = [15, 19, 17]`)
  - step 1 (mode=15), merge based on KNN and time including the temporal neighbors
  - step 2 (mode=19), merge based on KNN and time
  - step 3 (mode=17), merge based on KNN only


## [4.9.8] - 2019-8-27
### Added
- Memory usage is reduced by ~10x
  - Set `fSave_spkwav` to skip waveform saving for GUI (turned on by default)
  - Set `fSave_spkwav = true` if you want to do waveform similarity-based 
    merging by setting `maxavCor` = (0..1)
  - Turn it off in order to automatically merge using graph-connectivity only (kNN)


## [4.9.7] - 2019-8-26
### Added
- Improved automated merging
  - pre-merging operation uses KNNs for each cluster in each drift snapshot 
    to determine optimal merging partners.
  - For each drift snapshot in each cluster, the KNNs of the subset of spikes are found, 
    and the mode of their cluster IDs are set as the merging pair.


## [4.9.6] - 2019-8-21
### Changed
- Waveform-based merging uses averaged waveform obtained after SVD denoising
  - nComponents=3


## [4.9.5] - 2019-8-19
### Added
- `phy` exporter is added
  - run `irc export-phy myparam.prm`


## [4.9.4] - 2019-8-16
### Fixed
- `irc export-klusters` memory error is fixed when exporting a large recording

### Changed
- Faster automated merging based on knn
  - New knn-graph merging algorithm is used
 

## [4.9.3] - 2019-7-21
### Added
- Automatically remove refractory-period violating spikes
  - Useful when spike waveforms have multiple peaks (complex waveforms)
  - Introduced a new parameter (`spkRefrac_merge_ms=2`)


## [4.9.2] - 2019-7-9
### Changed
- Default setting changed to `maxWavCor=.985` (previously .975)
  - Based on the parameter optimization

### Fixed
- Backward compatibility is added to the output files pre-dating v4.8.0.


## [4.9.1] - 2019-7-8
### Changed
- `post_merge_mode=1` behavior is changed (currently the default option).
  - It only includes neighboring spikes having greater density values (16 nearest neighbors).
  - It previously included all neighboring spikes (8 nearest neighbors)
  - The new method computes average waveforms using higher density events (akin to the core set concept).
- Default setting changed to `maxWavCor=.975` (previously .985).

### Fixed
- Backward compatibility is added
  - Introduced in v4.8.0 (Jun 19, 2019) when `nSites_fet` parameter is introduced.
  - Parameters are now updated when an older file format is loaded.


## [4.8.10] - 2019-7-3
### Changed
- Default setting changed to `post_merge_mode=8`.
  - It performs waveform-based merging first and graph-based merging second.


## [4.8.9] - 2019-7-2
### Added
- New automated merging is added
  - Enable by setting `post_merge_mode=9`
  - This step performs merging based on outgroup/ingroup divergence ratio using the KNN graph.

### Fixed
- NeuroSuite Klusters exporter GUI error is fixed.
  - This happened when the `Klusters` directory didn't exist.


## [4.8.8] - 2019-7-1
### Added
- `nBatch_max_drift` parameter is introduced.
  - This sets the maximum number of batches for drift correction.
  - Default value: 8 (set to `inf` to disable).
- NeuroSuite importer is added for a given .xml parameter file.
  - `irc makeprm myparam.xml chanMap.mat`
    - `chanMap.mat`: Kilosort probe file format.
  - `irc makeprm myparam.xml myprobe.prb`
    - `myprobe.prb`: IronClust probe file format.
- NeuroSuite Klusters exporter is added to the manual GUI.
  - A new menu item is created in `File\Export to NeuroSuite Klusters`.

### Changed
- `Probe map` view (upper left window) restricts the sites from the same shanks.


## [4.8.7] - 2019-6-28
### Added
- `irc export-klusters myparam.prm` command is added
  - It exports IronClust output to klusters GUI compatible format (part of NeuroSuite).
  - Klusters is available from `http://neurosuite.sourceforge.net/`


## [4.8.6] - 2019-6-27
### Fixed
- `Traces` GUI '[F]ilter' command fixed
  - bsxfun mixed data format error is fixed
  - bsxfun_() function is implemented which checkes and converts the data format
  - gpuArray format is correctly handled


## [4.8.5] - 2019-6-24
### Added
- `Trials\Add uLED events` added
  - It accepts .mat file containing cmr_uled cell that contains start/stop time
  - The unit must be in seconds


## [4.8.4] - 2019-6-24
### Changed
- Default changed `fParfor=0`
  - It's not any faster to run parallel processing for merging spikes. Turned off by default now.

### Added
- `fGpu_rho=1`
  - rho computation is carried out in GPU by default
  - This is an approximation of KNN and less accurate than CPU-based KNN computation, 
    however, the performance is comparable.

### Fixed
- `irc manual` GUI: fixed auto-correlogram when one unit is selected.


## [4.8.3] - 2019-6-20
### Fixed
- `View\Show average waveforms on all channel` menu works for merged files

### Added
- `View\Show drift view` and `View\Show average waveforms on all channel` 
  menus toggle the views on and off.
- Ability to add multiple .nev files for loading stimulus timing
  - This feature is enabled when multiple raw recording files are specified in `csFile_merge`


## [4.8.2] - 2019-6-19
### Changed
- `mergeprm` command now supports an optional third argument for a new file name
  - `irc mergeprm myprm.prm files_to_append* newprm.prm` will create newprm.prm 
    in the folder containing myprm.prm

## [4.8.1] - 2019-6-19
### Added
- `irc mergeprm myparam.prm files_to_append*.*` will append other binary files to the prm files.
- `view\show global drift` is added to the manual GUI menu.


## [4.8.0] - 2019-6-19
### Changed
- 'maxSite' and 'nSites_ref' was removed
  - Specify `maxDist_site_um` and `maxDist_site_spk_um` instead to specify inner and outer radius.
  - Inner radius controls the feature dimensions and outer radius controls the waveform extraction.

### Fixed
- `drift view` in manual GUI excludes background spikes from other shanks
- `find_site_spk23_` error is resolved

### Added
- `ledprobe.prb` is added (Euisik Yoon's group)


## [4.7.8] - 2019-6-17
### Added
- x and y position views can be switched by selecting a tab
  - Drift view can be enabled by selecting "View\Show drift view" menu
- Checkbox appears next to the menu when `Drift view` or 
  `show average waveforms on ...` menu is selected
  - Checkbox disappears when the figure is closed.
- `split` and `merge` is enabled in the drift view
  - Select one cluster and press 's' to split by drawing a polygon.
  - Select two clusters and press 'm' to merge


## [4.7.7] - 2019-6-17
### Added
- "view\Show drift view" menu is added 
  - currently view only tool but will soon be able to merge and split in 
    the position view.


## [4.7.6] - 2019-6-14
### Added
- `post_merge_mode` 6 is added (default: 1)
  - experimental feature. this will compute average waveforms per unit per drift-batch
- `sort_mode` is reintroduced (default: 1)


## [4.7.5] - 2019-6-14
### Added
- FFT-denoising: Dual threshold scheme is added
  - `fft_thresh_low=4` is added. Noise peaks that crosses this is removed
    if it is adjacent to higher threshold crossing (`fft_thresh`)


## [4.7.4] - 2019-6-13
### Added
- `viSiteExcl` field is added. It excludes the sites from being analyzed. 
  - `viSiteZero` will set the site values to zero but it won't remove the sites. 
  - Re-indexing will occur in the site numbers after removal

## Fixed
- `irc makeprm myrecording.dat` now works without having to provide a probe file
  - Probe file can be specified later in the dialog box.


## [4.7.3] - 2019-6-13
### Changed
- `irc describe` command displays more info
  - Also fixed `fft_thresh` display error
- `fft_thresh=8` (previously 10)

## [4.7.2] - 2019-6-12
### Changed
- The definition of rho is reverted to defalut (`rho` :=1/d_knn)
  - (default) set `sort_mode=2` to use `rho:=1/dknn` and correct for density variation across time
  - set `sort_mode=3` to use rho:=`network in-degree`
  - set `sort_mode=1` to use `1/dknn` and not correct for the density variation across time 


## [4.7.1] - 2019-6-11
### Changed
- `delta_cut=1.1` (previously 1 by mistake)


## [4.7.0] - 2019-6-11
### Changed
- k-nearest-neighbor computation for computing `rho` is done exactly using CPU by default.
  - CPU-based computation is ~2x slower than GPU but exact.
- `rho` is now computed using the number of neighbors pointing to each event
  - `kNN` matrix is used to compute the number of neighbors pointing to each event
  - Previously `rho` := 1/d_knn where `d_knn` is the k-th neighbor distance
- Introduced a delta cut-off parameter called `delta_cut=1.1`.
  - Previuosly this value was fixed at 1.
  - The new value reduced cluster oversplitting.


## [4.6.9] - 2019-6-11
### Changed
- Density peak merging behavior is changed
  - peaks are merged only when both peaks are connected
  - KNN (k-nearest-neighbors) relations are not symmetric thus both relations need to be checked.
- `fft_thresh=10` 
  - Previously it was 6, but it was too agressively removing noise and affecting signal spectrum.


## [4.6.8] - 2019-6-10
### Changed
- S_clu_peak_merge_: degree of seperation of up to 2 is allowed between density peaks
  - previously separation of 1 was allowed.
  - 16 nearest neighbors were considered (previously all knn neighbors, knn=30)


## [4.6.7] - 2019-6-10
### Changed
- `fft_thresh=6` (previously 10).
  - Detrending of the power is done using sliding median filter (medfilt1.m)


## [4.6.6] - 2019-6-7
### Changed
- raw data and `_spkraw.jrc` files are stored in its original format
  - raw spike waveforms are now stored in the original format regardless of `fft_thresh`.
  - previously, raw spike waveforms were saved after performing `fft_thresh` operation.
- `fft_thresh=10;` by default (previously 0)
  - `fft_thresh` operation is now performed in the `fft_filter` function
  - Only one fft and ifft operation is performed


## [4.6.5] - 2019-6-5
### Changed
- `vcDataType_filter` sets the internal number representation of 
  the raw and filtered waveforms
  - The new default is `single` (previusly `int16`)
- `single` data type for the raw recording doesn't get scaled down to int16 format by default.
  - If `vcDataType_filter=single` is set.


## [4.6.4] - 2019-6-4
### Fixed
- `single` data type for the raw recording is correctly scaled to int16 format.
  - max and min value of the recording is used to fill +/-2^14 int16 range
  - The actual microvolts scaling may be off using this approach.
    However, most single-type recordings are from simulated studies 
    for a ground-truth validation purpose, where uV scaling isn't necessary.


## [4.6.3] - 2019-5-30
### Changed
- `maxWavCor` set to .985 (previously .975)
- `spkLim_ms` set to [-.25, .75] (previously [-.25, 1])


## [4.6.2] - 2019-5-29
### Changed
- `gain_boost` parameter is removed.
- `scale_filter` paramter is set to 8.


## [4.6.1] - 2019-5-28
### Changed
- `scale_filter` parameter is set to 10 (previously 1)
  - This parameter scales up the int16 representation after filtering
  - Substantially improves lower SNR performance by preventing data clipping
- `maxWavCor = .975` (previously .99)
  - Controls automated merging threshold

### Added
- `batch_sec_drift`: controls time duration per cluster batch size
- `step_sec_drift`: controls time step per anatomical spatial distribution calculation


## [4.5.11] - 2019-5-17
### Fixed
- monotrode (single-channel electrode) sorting error is fixed


## [4.5.10] - 2019-5-16
### Fixed
- `convert-mda` error is fixed for exporting to SpikeForest groundtruth dataset (.mda format)
- Unused parameters are disabled in `default.prm`


## [4.5.9] - 2019-5-14
### Changed
- `default.prm` paramters updated based on parameter search using SpikeForest dataset
  - `maxWavCor=.99`, previously .98
  - `nTime_clu=2`, previously 4
  - `vcCommonRef='mean';`, previously 'none'

### Fixed
- `fft_filter` error fixed when `vcFilter='wiener';` option is used.


## [4.5.8] - 2019-5-7
### Fixed
- `Show averaged waveforms on all channels` error fixed
  - Previously setting `vcFilter='bandpass';` caused error.
- `template file not found` error is resolved.
  - Python 'None' type is now correctly translated to `[]` empty array in Matlab.


## [4.5.7] - 2019-5-1
### Added
- manual GUI: `Show raw waveforms` toggle will influence the all-channel waveform view
  when `view\Show averaged waveforms on all channels` is selected.
- manual GUI: Keyboard shortcut is assigned to `reset window positions`. 
  - Press `1` (number one) to activate. 
  - Focus is returned to the main window after rearranging.
  - This features doesn't eactly align window boundaries in Linux.
- `irc scoreboard`: New GUI tool that works with parallel computing resources
  to find optimal set of sorting parameters
- `irc convert-mda-buzsaki [in] [out]`: convertion utility that exports to mda format
  - [in]: .mat file that follows `*.spikes.cellinfo.mat` format.
  - [out]: output directory
  - Different shanks are stored in different subfolders.


### Fixed
- manual GUI: all figures will close when the main window is closed.


## [4.5.6] - 2019-4-19
### Fixed
- Error is resolved when using CPU (`fGpu = 0`);
- Wiener filter kernel is corrected. now works with bandpass filter if `freqLim` is set
  - set `freqLim=[300,nan]` to combine with a high-pass filter at 300 Hz
  - set `freqLim=[nan,6000]` to combine with a low-pass filter at 6000 Hz
  - set `freqLim=[300,6000]` to combine with a band-pass filter at 300-6000 Hz


## [4.5.5] - 2019-4-16
### Fixed
- vcFilter = `wiener` option is revised
  - The filter is symmetric in the frequency domain
  - It's computed once and reused throughout the recording

### Changed
- Default value changed: `nPcPerChan=3`


## [4.5.4] - 2019-4-15
### Added
- vcFilter = `wiener` is added (fft_filter.m)
  - it's now set as default
  - It applies high-pass filter and wiener filter
  - low-pass limit on the `freqLim` is ignored


## [4.5.3] - 2019-4-13
### Added
- `fftdiff` filter option is added (band-limited differentiator)
  - it was set as a default for `vcFilter`
  - `freqLim` and `freqLim_width` controls the low-pass filter cut-off
  - high-pass filter cutoff is ignored when `fftdiff` is set.


## [4.5.2] - 2019-4-12
### Changed
- Drift correction performance improved by changing parameters and 
  automated merging algorithm
  - nPcPerChan: 1 -> 2
  - freqLim: [300, 4000] -> [300, 5000] 
  - merging algorithm: 
- `run_irc.m` changed to use cached results

### Added
- `run_ksort2.m` added
  - Usage: `run_ksort2 dir_in dir_out (params.txt)`


## [4.5.1] - 2019-4-11
### Changed
- Drift correction speed and performance improved


## [4.5.0] - 2019-4-11
### Added
- Drift correction performance improved
  - Validated using Kilosort2 simulated data, outperformining Kilosort2


## [4.4.9] - 2019-4-10
### Fixed
- `irc convert-mda-english` bugfix
  - Intracellular time range is corrected in `summarize recording` menu
    - It uses `spkLim_ms_gt` constant stored in `default.cfg` for the spike 
      waveform time range in milliseconds
  - `convert` menu correctly exports .mda files containing extracellular 
    channels (`raw.mda`) and stores intracellular channel (`raw_true.mda`)
    in correct dimensions (1 x nSamples).

### Added
- `ms_bandpass` filter is optimized by making efficient use of GPU memory.


## [4.4.8] - 2019-4-9
### Added
- Merged Zach Sperry's contribution on spike triggered visualization
- Added GUI menu (view\Show averaged waveforms on all channels)
  - It displays median waveforms on all channels for a selected cluster
- `fZach` switch is added. It's turned off by default. 
  - To enable, create `user.cfg` file in ironclust folder and set `fZach=1;`


## [4.4.7] - 2019-4-8
### Added
- Added GUI functionalities on `irc convert-mda-english`
  - `summarize recording` menu is added, which plots average waveforms on 
    intracellular and extracellular channels
  - Ability to customize probe files on the GUI table

### Fixed
- fastplot.m error fixed when plotting multiple channels


## [4.4.6] - 2019-4-1
### Fixed
- Singular dimension error fixed when no adjacent channel is found per spiking event
  - `squeeze_` function explicitly takes which dimension to squeeze


## [4.4.5] - 2019-4-1
### Added
- Fast plot of large dataset is implemented in `fastplot.m`
  - It uses x-axis crop subsampling to reduce the number of points to be plotted
  - It supports multiplot and amplitude rescaling
- `irc convert-mda-english` GUI tool
  - It is used to export paired recording groundtruth from Dan English's dataset
  - The settings are stored in `dan_english.cfg`


## [4.4.4] - 2019-3-22
### Fixed
- `fRemove_duplicate` works correctly for channels with different noise baseline
  - the peak channel location is determined based on the waveform amplitude divided by detection threshold'
- Crash is averted if only one cluster is found


## [4.4.3] - 2019-3-13
### Added
- `irc convert-h5-mda recording.h5` works on Matthias Henning lab format

### Changed
- Drift correction now uses centroid vs amplitude histogram
  - Previously used channel vs amplitude histogram (quantized by electrodes)


## [4.4.2] - 2019-2-22
### Fixed
- Filtering and detection operation became more GPU-memory efficient.

### Changed
- Default parameter changed: `vcFilter = 'bandpass';`
  - previous value was `vcFilter = 'ndiff';`


## [4.4.1] - 2019-2-22
### Added
- GUI cluster splitting based on the waveform feature (contributed by Zach Sperry)
  - This feature is enabled when a user selects `manual` splitting and `Waveform` projection.
  - Currently selects the peak channel and a user can draw a line to group spikes.
  - Future version will allow user to select a channel using arrow keys.
- ViSAPy validation importer added (used by YASS paper)
  - `irc convert-mda-yass [input_dir] [output_dir]`
 
### Fixed
- SpikeForest interface works correctly for a float32 (single) type
  - It correctly reads `scale_factor` field from `params.json` file used in 
    spikeforest format


## [4.4.0] - 2019-2-8
### Added
- Splitting in the PSTH window
  - Right-click and select either horizontal or vertical cut on the raster window


## [4.3.9] - 2019-2-7
### Added
- `irc makeprm myfile.ns5 myprobe.prb` command is supported for nsx format
  - previously user had to run a command `irc import-ns5 myfile.ns5 myprobe.prb`
- `irc makeprm myfile.ns5` now takes myprobe.prb from the default_prm file
  - `default_prm` points to the default parameters file (`default.prm by default`)
  - `default.cfg` or `user.cfg` stores `default_prm` setting
  - `user.cfg` settings takes precedence over `default.cfg` and doesn't get 
    overwritten during updates


## [4.3.8] - 2019-2-7
### Fixed
- PSTH count histogram non-zero error is fixed
  - Previously, PSTH histgram count was non-zero where it does't contain any spike.

### Added
- `irc github` opens a ironclust github website


## [4.3.7] - 2019-2-6
### Fixed
- Plotting PSTH doesn't take the focus away from the main waveform window
- PSTH histogram displays the correct unit

## [4.3.6] - 2019-2-4
### Added
- New template added `static_diff_template.prm`
  - static probe with a differentiation filter

## [4.3.5] - 2019-1-31
### Added
- `irc summarize-study path_study` outputs a summary table of ground truth units
  - It works for paired recording studies with one ground truth unit.
  - path_study contains folders containing either firings_true.mda or _gt1.mat files

### Changed
- Adding PSTH channel allows user to specify min_interval_psth

### Fixed
- Floating point data type doesn't saturate when gets converted to int16
  - Now subtracts the median of each channel to center the data
  - This prevents data saturation


## [4.3.4] - 2019-1-25
### Added
- `irc convert-mda-crcns` opens a GUI that lets users to select and export Buzsakilab 
  groundtruth dataset (paired intra/extra recordings) from crcns.org (hc-1)
  - syntax: `irc convert-mda-crcns path_in path_out`
  - Data source: https://crcns.org/data-sets/hc/hc-1 (log-in required to download)

### Fixed
- `irc install` compiles CUDA code using the GPU compute capability 
  installed on the local system
  - Previously it compiled for compute capability 3.5 (Kepler and above) 
    for backward compatibility
  - run `irc install 0` to compile CUDA code for older GPU cards for 
    backward compatibility (compute capability 3.5)

### Changed
- syntax updated for `run_irc path_in path_out template_file`
  - template_file can be either full file name or simplified name
  - If template file is `tetrode` the full name is auto-complted to `tetrode_template.prm`
  - Autocompletion only works when the simplified name does not contain a file extension    


## [4.3.3] - 2019-1-8
### Added
- `Trials\add PSTH channel` menu validates the file format
  - Operation is aborted if no events are found.
  - Number of events are shown after successful validation
- User file selection dialogs were added for `add analog channel` and 
  `add PSTH channel` under `Trials` menu

### Changed
- Defalut changed: `thresh_core_knn = .9` (previously .875)
  - improves real groundtruth performance
- Added a paramter for post-clustering: `fTemplateMatch_post` (default value: 1)
  - controls template matching for the outer part of the cluster cloud.

### Fixed
- `manaul GUI`: PSTH bar plot is now correctly centered.
  - Previous version displayed incorrectly time-shifted bars on the upper plot.


## [4.3.2] - 2018-12-13
### Added
- Settings menu is added for the `PSTH` plot
- `Open prm folder` command is added under `File` menu.
  - This opens a system file browser containing the .prm file
- `Plot PSTH` command is added under `Trials` menu.
  - User needs to add a `PSTH channel` under `Trials` menu and select it

### Fixed
- `irc mcc` and `run_irc.m` scripts were made compatible with the spikeforest singularity container.
  - `geom.prb` file is generated in the output folder
  - `tetrode_template.prm` file is correctly located in the irc path
  - Plot error is caught if plots can't be generated inside of the container
- run_irc.m (used by the python wrapper) and run_irc executable files were made consistent
  - Syntax: `run_irc dir_in dir_out myparam.prm`


## [4.3.1] - 2018-11-23
### Added
- `preview` command under `Trials` menu displays the channel label and unit
- Blackrock Neuroshare reader upgraded to v4.5.3.0
  - Previous version: v4.5.2.0
  - https://github.com/BlackrockMicrosystems/NPMK/releases/tag/4.5.3.0
- `split` menu now displays error bounds (mean+/-SD) of the unit waveforms 
  and channel numbers

### Fixed
- `add analog` command under `Trials` menu
  - Correct channel number is loaded when the channel number doesn't match 
    the file storage order
  - Previously it was selecting a channel using an electrode number

### Changed
- `add analog` command under `Trials` menu asks for `Gain` factor instead of
  the bit conversion factor
  - The timeseries later gets multiplied by `uV_per_bit` and `Gain` factor


## [4.3.0] - 2018-11-21
### Fixed
- `irc manual` analog channel add error
  - makeStruct() renamed to makeStruct_()


## [4.2.9] - 2018-11-19
### Added
- `trials` GUI menu: Added a menu to import event information from .nev file
  - Contributed by Zach Sperry

### Fixed
- `uimenu_()` function is created to make uimenu() version-neutral
  - Since Matlab R2017b `Label` field is renamed to `Text` field
- Parallel processing is disabled during GUI operation (`irc manual`)
- `auto-merge` command in `irc manaul` correctly recomputes quality scores and position


## [4.2.8] - 2018-11-13
### Added
- `trials` menu is added in `irc manual` GUI
  - renamed from `plots`
  - ability to add multiple trial types
    - `event` trial type: user supplies 'on' and 'off' times in seconds
    - `analog` trial type: user supplies file name and channel number 
  - users can select, edit, remove and preview the trial variable
  - `plot` menu correlates the unit firing rate with a selected trial variable

### Fixed
- `memory` function in Mac OS is added, which returns the free memory.


## [4.2.7] - 2018-11-07
### Fixed
- `addpath` correctly finds the current source path when a relative path is given
  and directory is not found.

### Added
- `irc addpath` adds the current source path to the search path


## [4.2.6] - 2018-11-07
### Fixed
- `sbatch-mda` command creates .log files in the destination folder
  - `template_prm` file is found correctly if absolute path is given

### Added
- `sbatch` in `default.cfg` stores cluster run settings
  - see examples: `sbatch_bionet.m` and `sbatch_magland_synth.m`
- `copyfrom-voms` command copies files (`copyfrom_voms`) in (`path_voms`) under default.cfg


## [4.2.5] - 2018-11-02
### Added
- `irc sbatch-mda` can run in non-blocking mode to launch multiple sbatch script
  - see `sbatch_magland_synth.m` as an example


## [4.2.4] - 2018-11-02
### Added
- `irc sbatch-mda` uses Flatiron compute-cluster (sbatch + disBatch)
  - syntax: `irc batch-mda input_dir output_dir (template_file)`
- `irc batch-mda` runs 
  - syntax: `irc batch-mda input_dir output_dir (template_file)`


## [4.2.3] - 2018-11-02
### Added
- Dockerfile added
- `irc mcc` command creates run_irc binary using matlab compiler
  - matlab compiler licence is required
  - usage: `./run_irc 

### Changed
- All matlab scripts are moved to /matlab


## [4.2.2] - 2018-10-29
### Fixed
- mcc output is supported (matlab compiler)
  - addpath is protected 
- `validate` plot and parfor errors captured during compute-cluster mode
- fGpu flag is correctly applied for clustering


## [4.2.1] - 2018-10-23
### Changed
- vcDetect='xcov' uses improved detection method based on time-smoothed xcov

### Added
- `irc mda-cut mylist.txt` command added
  - It copies a folder containing raw.mda files after taking half the time-series
  - It outputs to new folders ending with `_half`

### Fixed
- Validation figure names are set correctly
  - `filename_timestamp_` function adds 'yymmdd-HHMMSS` before the file extension


## [4.2.0] - 2018-10-19
### Changed
- vcFilter='bandpass' now uses fft-based bandpass filter 
  - copied from `ms_bandpass_filter.m` from www.github.com/magland/mountainsort/old
  - this function utilizes GPU if available
  - signal processing toolbox is not required to use the `bandpass` filter.
- Spike waveform is subtracted by the average across all channels 
  - It uses 'meanSubt_spk_' function
  - Previously spike waveforms were subtracted by each channel's time-averages
  - Less false positives since average of each channel is not forced to zero

### Added
- `irc verify` command displays and saves a figure that embedes the code (irc.m), parameters, and scores
- Faster KNN computation 
- `irc gtsort` command uses spike timing from the ground-truth file (vcFile_gt)
- `irc showgt` command plots clusters based on the ground-truth file (vcFile_gt)


## [4.1.9] - 2018-10-12
### Added
- npy-matlab added: https://github.com/kwikteam/npy-matlab

## [4.1.8] - 2018-10-12
### Changed
- Removed 'sample.bin' and 'sample.meta' under ./prb folder

## [4.1.7] - 2018-10-12
### Bugfix
- 'makeStruct not found error' fixed

## [4.1.6] - 2018-10-11
### Changed
- `irc export-mda` command updated
  - `irc export-mda myparam.prm myoutput.mda` option created
  - `myoutput.mda`: output file fullpath

## [4.1.5] - 2018-10-10
### Changed
- Default values changed to optimize out-of-box performance
  - `vcSpkRef='none';` which performs local referencing
  - `maxDist_site_spk_um=100;` (previously 75) which uses more number of channels for waveform features
  - `thresh_core_knn=.9;` (previously .75) which is used for post-hoc template match operation

### Fixed
- Blackrock format (.ns5) compatibility issues resolved

## [4.1.4] - 2018-10-03
### Changed
- Default value changed: `maxWavCor = .97;`
  - Previously `maxWavCor = .98;`
- Waveform merge is now done by comparing `nTemplates_clu` (default:100) 
  representative waveforms per cluster between cluster pairs centered at
  the same site.

### Added
- Post-cluster correction is performed on the outer part of clusters using template matching
  - This operation is performed if `fTemplateMatch_post=1` (default is 1)
  - Templates are chosen using events in the core part of the cluster
  - Core vs. outter part of the cluster is separeated using 
- `makeprm-mda` can now import a list of raw.mda files in a .txt file
  - `irc makeprm-mda mylist.txt mytemplate.prm` creates mylist.batch file
    where .txt file is located
  - irc locates "geom.csv", "params.json" files listed under the same folder 
    as in "raw.mda"
  - irc locates "firings_true.mda" file for the ground-truth comparison if available


## [4.1.3] - 2018-09-26
### Changed
- xcov based detection is now 100x faster

## [4.1.2] - 2018-09-25
### Added
- Cross-covariance (xcov) based detection method is added
  - set `vcDetect='xcov';` to enable this feature
  - This is significantly more sensitive to detecting low SNR spikes
  - This currently runs much slower than the minimum peak detection method (`vcDetect='min'`)
- xcov based feature extraction and clustering method
  - set `vcCluster='xcov';` to use cluster using this method
  - Currently supporting small number of channels to be clustered together

## [4.1.1] - 2018-09-20
### Added
- Tetrode-optimized clustering for small number of channels, smaller amplitudes
  - `vcCluster = 'waveform'`
  - This clustering will create nC*nC*nD features
    - nC: number of channels
    - nD: number of time delays (`nDelays_wavcov` parameter, default 4)
  - Currently using CPU, GPU-version will follow (~40x speed-up expected)
- Memory-optimized KNN using CPU when using `vcCluster='drift-knn'`

### Fixed
- "icl not found error" resolved when capturing the rho-delta plot


## [4.1.0] - 2018-09-19
### Added
- Mountainlab integration: exposed two new parameters
  - `pc_per_chan` controls number of principal components per channel
  - `merge_thresh` controls waveform-based cluster merging (0..1)


## [4.0.9] - 2018-09-19
### Added
- mountainlab pipeline suppport added.
  - `firings_true` input is accepted for ground-truth validation.

### Fixed
- `adjacency_radius=-1` is interpreted as using all the sites. 


## [4.0.8] - 2018-09-19
### Added
- Auto-merging before the waveform merging using knn overlap.

### Changed
- maxWavCor=.98 default value (waveform-based merging threshold
  - Previously .975 was default
- Spikes belonging to the clusters having less than min_count (default 30) are removed
  - These spikes get assigned with "0" cluster value (noise cluster)
- Average waveforms for pair-wise merging is calculated using half the waveforms
  having smaller RMS errors compared to the median waveforms (trimmed-mean)

### Fixed
- Cluster merging error is resolved in manual GUI.


## [4.0.7] - 2018-09-14
### Added
- `frac_equal_knn` parameter is used to identify events whose fraction of 
  their K-nearest neighbors having the same cluster membership.
  - Events whose equal-knn fraction below the threshold (`frac_equal_knn=.5`) 
    are reassigned using the median cluster value of their KNN.
- `frac_shift_merge` parameter added: % of maximum time shift for waveform merge consideration
- Preview GUI asks users to input the filter setting  
  - Differentiator filter `vcFilter='ndiff'` asks for `nDiff_filt` parameter
  - Band-pass filter `vcFilter='bandpass'` and fir1 filter asks for `freqLim` parameter
- `verify` command adds interactive plots where users can click on the cluster to see 
  more on how clustering failed
- Bursting spike and overlapping spike metrics added to `verify` and `batch-verify` command
  - `S_burst` and `S_overlap` struct is included in `_score.mat` output file. 

### Changed
- "fWavRaw_merge" default is now 0 (previously 1)
  - Speed enhancement for auto-merging clusters when fWavRaw_merge=0
  - Average raw waveforms are still computed at the end of merging instead
  of at each merging step.
- Mean cluster waveform is computed at three vertical depths when "fDrift_merge=1"
  - Previously this was only done when fWavRaw_merge=1
- `convert-h5-mda` tool now outputs the maximum channel index
- CUDA compiled for NVIDIA GPU version SM3.5 for backward compatibility
- Defaulg value changed: `maxWavCor = .975` (waveform merge criteria for post-clustering stage)
- Default changed `nDiff_filt=3` (previously 2)
- Deafult changed `maxDist_site_spk_um=100` (previously 75)

### Fixed
- viSiteZero (bad sites) is handled properly such that they are never included as neighboring channels
- Fixed `makeprm` command ".bin file not found" error.


## [4.0.6] - 2018-08-21
### Added
- Convertion tool from .h5 (Boyden lab) to .mda (Mountainlab)
  - `irc convert-h5-mda [mylist.txt] (output_dir)`
  - [mylist.txt] contains a list of .h5 files (shouldn't be _raw.h5 nor _spikes.h5)
  - [mylist.txt] can bs substituted with an individual .h5 file (e.g. myfile.h5)
  - (output_dir) is optional and if empty a new folder is created where .h5 resides
- Import tool from .h5 (Boyden lab)
  - `irc makeprm [mylist.txt] (mytemplate.prm)
  - [mylist.txt] contains a list of .h5 files (shouldn't be _raw.h5 nor _spikes.h5)
  - [mylist.txt] can bs substituted with an individual .h5 file (e.g. myfile.h5)
  = It creates mylist.batch file where mylist.txt resides
  - Output files are created where .h5 files reside
- New command added for batch-creating .batch file and .prm files when
  a list of binary files are provided in "mylist.txt".
  - `irc makeprm mylist.txt mytemplate.prm`
  - creates mylist.batch containing .prm files and individual .prm files 
  each pointing to the mytemplate.prm file
- Importing .raw file (multichannel systems format) using `makeprm` command
  - `irc makeprm myfile.raw myprobe.prb`
  - This extracts a header portion of .raw file, creates `.meta` file 
  (`myfile.meta`) determines the number of header bytes to skip.
- `irc changelog` displays changelog.md on Matlab editor and on the web 
  markdown viewer.
### Changed
- change_log.txt changed to changelog.md and the formatting changed to 
  markdown format (`.md`)

### Fixed
- "fGpu" flag correctly controls whether CPU or GPU is used for clustering 
  when vcCluster = "drift-knn" is set.


## [4.0.5] -2018-08-14
### Added 
- p_ironclust runs ground-truth validation if "firings_true.mda" is present 
  in the directory containing the raw mda file (e.g. "raw.mda").
### Fixed
- tetrode_template.prm default changed to work with the new default 
  clustering algorithm (vcCluster = "drift-knn")


## [4.0.4] - 2018-08-03
### Fixed
- tetrode_template.prm filter option changed to "fir1, [300, 6000] Hz"
  It works even if the Signal Processing Toolbox is absent using 
  "fir1_octave.m" function. Previously used bandpass elliptic filtfilt.  


## [4.0.3] -2018-08-03
### Fixed
1. fParfor flag correctly controls parallel computation toolbox use.
  fParfor=0 if invoked by mountainlab-js


## [4.0.2] -2018-08-03
### Fixed
1. Crash in the cluster prevented by setting fSavePlot_RD=0 when created 
  by the mountainlab-ps called by p_ironclust.m


## [4.0.1] -2019-08-03
- Editor opening disabled in the cluster node.


## [4.0.0] -2018-08-02
### Added
- New clustering algorithm: New density estimate using K-nearest neighbor distance (fast estimation).

---

## [3.3.0] -2018-05-03
### Fixed
- mda integration error (-p prm_template_name).

## [3.2.9] -2018-05-03
### Added
- Clustering using CPU. GPU is optional and Parallel computing toolbox is not necessary.
- Set `fGpu=0` to disable GPU use. If GPU processing fails, CPU wlil be 
  automatically used for clustering.

## [3.2.8] -2018-04-20
### Added
- mda support added (mountainlab)
- Output directory option added to "makeprm" command
  -`irc makeprm [myrecording.bin] [myprobe.prb] (mytemplate.prm) (myoutput_dir)`
  - `myoutput_dir`: if specified, output files to be written there. 
  If not specified, the default output directory is where `.bin` file is located
  - `mytemplate.prm`: if specified, this overwrites parameters in `default.prm`
- `irc makeprm-mda [raw.mda] [myprobe.prb/geom.csv] (params.txt) (myoutput_dir) (mytemplate.prm)`
  - It creates a `.prm` file from the mountainlab format (`.mda`)
   
## [3.2.7] -2018-03-20
### Added
- default.cfg (default configuration file): Added "reset_path" = 1 by default. 
  Set to 0 to disable reset matlab path when running JRCLUST
- `irc makeprm [myfile.rhs] [myprobe.prb]` converts .rhs to .bin and .meta files

### Fixed
- `matlab -nodesktop -nodisplay` safe: edit and figure creation errors are captured.
- `irc compile` compiles GPU codes (.cuda) to .ptx files using the installed GPU
- Support for `.rhs` format (Intan RHS series)  
  - `irc makeprm [myrecording.rhs] [myprobe.prb]`

## [3.2.6] -2018-02-27
### Bugfixes for Manual GUI
- Running a new manual clustering after closing a previous one is handled correctly.
  - Previously it required `irc clear` to clear the residual memory from the previous run.
- Faster plotting after merge, split, delete operations.
- Splitting clusters with less than three doesn't crash anymore.

## [3.2.5] -2018-01-05
### Fixed
- Centered cluster average waveforms are properly saved. 
  - `trWav_spk_clu` for filtered waveforms and `trWav_raw_clu` for raw waveforms
  - `tmrWav_spk_clu` and `tmrWav_raw_clu` have worked properly, which saves waveforms at all sites.
- Selecting 'vpp' (peak-to-peak voltage) calcultes different values
  depending on whether you are displaying filtered or raw waveforms.
  When 'pca' or 'ppca' (private pca) features are selected, PCA values are computed for the filtered waveforms whether or not the raw waveforms are displayed.
- Manual GUI: The x-axis labels on the waveform view are displayed without overlapping for smaller sized monitors.
  To avoid overlapping, the text size is decresased and text angle is rotated (new feature since Matlab 2016b).

### Added
- Manual GUI: 'ppca' (private PCA) option is added under "Projection" menu.
  When one unit is selected, filtered spike waveforms are projected to the first and second principal vectors of the selected unit for each site.
  When two units are selected, filtered spike waveforms are projected to the first principal axes of the two selected units for each site.
- irc export-chan myparam.prm chan-list" exports channel(s) on "chan-list" to the Matlab workspace and to a binary file (_ch#.jrc).
  To export more than one channel, use "," or ":" without space. 
  - `irc export-chan 33,34,35` or `irc export-chan 33:35` to export channels 33-35.
- Export LFP waveforms (ordered by the site numbers) to the Matlab workspace
  - `irc export-lfp`
  - mnLfp" has a nSamples x nSites dimension, and stores the LFP waveforms in the original data type (16-bit integer by default).
  - `mrSiteXY`"` has a nSites x 2 dimension, and stores x and y coordinate for the sites (in micrometers)
  - To extract a column of sites located at x=0, run `mnLfp(:, mrSiteXY(:,1)==0)`

## [3.2.4] -2018-01-04
### Fixed
- Manual GUI: projection display is correctly updated.
  - Previously the second unit selection (shown in red) required two right clicks to update.

## [3.2.3] -2018-01-03
### Fixed
- incorrect "gather" function was called for R2017b when Bigdata toolbox is installed.
- Faster and more robust fft_clean operation (run if `fft_thresh` > 0). 
  - It can handle a limited GPU memory.
- `nLoads_gpu` parameter is introduced, which is the ratio of the system memory (RAM) to GPU memory.
  If your GPU memory size is 12 GB and RAM is 96 GB, set nLoads_gpu = 8.
  Increase this value if GPU error occurs (default: 8), do not change the default value otherwise.
  `nSamples_gpu` is now deprecated, and automatically calculated from `nLoads_gpu` parameter.

---

## [3.2.2] -2017-12-30
Bugfix: PSTH histogram plot error is fixed and two PSTH plots are displayed
  top and bottom when two units are selected.

## [3.2.1] -2017-12-29
Old multishank probe file format is supported (v1, jrclust.m). 
  'cviShank' which is a cell containing a list of sites is converted to 'shank' format
FFT clean operation is faster now (used if "fft_clean" parameter is greater than 0)
Manual GUI: Number of spikes per unit displayed the x-axis label can be toggled by pressing 'n' 
Manual GUI: performance improvement
  Faster update when a unit is selected
  Faster PSTH plot
  Default log length decreased to 5 from 10 previously ('MAX_LOG' parameter).
New command: "jrc quality myparam.prm"
  This command displays the unit quality info and exports to a csv file (ends with "_quality.csv").

## [3.2.0] -2017-12-21
Bugfix (Manual GUI): Unit selection cursors do not disappear anymore after file saving.
  This has been an issue since v3.1.8 with an introduction of saving the 
  cluster selection plot (rho vs. delta for the DPCLUS clustering algorithm).
Bugfix ('maketrial' command): Supports binary files with header information and/or
  multiple combined recordings if multiple binary files are set in "csFile_merge" parameter.
Speed improvement: "delete/merge/split" operation for Manual GUI got significantly faster
  Removed unnecessary visual updates and text box display for the unit counts, which is now
  shown on the x-axis as "unit# (counts)". The text box indicating the unit count used to
  block the waveforms and took a long time to plot.
"jrc clear mybatch.batch" will clear the output of the .prm files listed the batch text file (.batch).
LFP-related commands are migrated from jrclust.m (ver. 1)
  "jrc import-lfp" extracts LFP signal by subsampling the raw traces or by loading 
    the LFP file (.imec.lf.bin) for Neuropixels recordings.
  LFP file name (.lfp.jrc) is automatically generated from .prm file and 
    `vcFile_lfp` parameter is updated in .prm file.
  Subsampling is done by averaging the time bin (P.nSkip_lfp) to prevent aliasing 
    as opposed to simply skipping samples.
  "jrc traces-lfp" plots the LFP traces with a larger time window by a factor of P.nSkip_lfp.
    This works similarly to "jrc traces" command for viewing wide-band traces.

## [3.1.10] -2017-12-15
New parameter: `nSamples_gpu` which sets the number of samples to be loaded to GPU memory.
  Current default is 250,000, and previous default was 1000,000 which could not be changed.
  Reduce this number if you run out of GPU memory.  
  This value is used if P.fft_thresh > 0 or "P.vcCommonRef=='whitening'".
New parameter: 'fSpatialMask_clu' is disabled by default (previously enabled)
  Spatial mask applies spatially decaying weight as a function of distance from the center site
  `maxDist_site_um` is used as the half-distance scale (1/2 weight is multiplied at that distance)
  Spatial masking is disabled if you have less than "min_sites_mask" number of sites per spiking events.
n-trode geometry is now supported (n>=1). You need to set each tetrode group as separate shanks by setting 
  `shank` parameter in the .prb file for each site number. 
  (e.g.) if site 1:4 belongs to tetrode 1 and site 5:8 belongs to tetrode 2, then set "shank = [1 1 1 1 2 2 2 2];";
Bugfix: Spikes near the beginning of the recordings are now properly merged (spkMerge_single_() function).
New options for common-mode rejection: 'vcCommonRef' now supports 'median' and 'whitening'.
  Setting `vcCommonRef` = 'whitening' will perform spatial whitening between channels.
Trace GUI: color of spikes remains consistent between time window switching.
Manual GUI: Bugfix: Combination of selecting 'raw waveform' view in 'pca' unit works properly.
New feature: "maketrial" command
  This function generates trial timestamps to plot PSTH (either rising or falling edges) and generates .mat file
    containing timestamps and sets 'vcFile_trial' field.
  For IMEC probes, use "maketrial-imec" to extract event input on the digital input pins (1-16). 
    Specify the digital input channel by setting "dinput_imec_trial".

## [3.1.9] -2017-12-06
Bugfix: v3.1.8 incorrectly points to v3.1.7 on Github. Please download this version instead.
Bugfix: Opening previously sorted data (v3.1.1) to the latest vresion now works.
  Parameter version upgrade function (upgrad_param_()) now handles new way to compute 
  maxSite and nSites_ref parameters.

## [3.1.8] -2017-12-06
New feature: Rho-Delta plot is captured as .png file when _jrc.mat file is saved
  Setting the parameter 'fSavePlot_RD = 0' will disable the figure capture (enabled by default).
  You may adjust 'rho_cut' and 'delta1_cut' paramters based on the Rho-Delta plot.
  Run 'jrc plot-rd' command to view the Rho-Delta plot.
Manual GUI
  Unit position view displays number of spikes in a selected unit
  Bugfix: Clusters with less than 6 spikes after manual splitting doesn't crash any more.
  Bufix: Second manual splitting erases the first drawn polygon (splitting in the waveform view)
  Bugfix: Switching between filtered and raw waveforms now use the same set of colors
Improved error reporting: "jrc_error_log.mat" records detailed error log.
Raw waveform duration is reverted back to 2 msec
  spkLim_raw_factor = 2 (default) will use 2x of the time range (spkLim_ms) used for the filtered waveform 
  spkLim_factor_merge = 1 (default) will use 1x of the spkLim_ms (1 msec default) for the waveform based auto-merging
New default: maxDist_site_spk_um = 75 (changed from 70) for improved performance for probes with 25 um site spacing

## [3.1.7] -2017-12-01
Bugfix: Feature projection view
  Correct y-axis values (PC2) are plotted when splitting on the diagonal 
  site pairs when PCA feature is selected. 
Speed improvement: fft_clean operation (fft_thresh parameter) is performed faster.
New feature: Auto-splitting now supports manual splitting in three PC dimensions.
  Auto-split now offers three choices: {'yes', 'no', 'manual'}.
  Select "manual" and choose the projection ('PC1 vs PC2', 'PC3 vs PC2', 'PC1 vs PC3').

## [3.1.6] -2017-11-20
Bugfix: Nested S0 struct is removed when older version is loaded.
Drift correction is improved. Three average waveforms are computed by grouping
  spikes from similar center-of-mass. Now using all spikes from the same unit
  previously used spikes that were centered at the peak site only.
New parameter: fSpatialMask_clu (1 by default) applies spatial masking to 
  the features. Half decay length set by the maxDist_site_um.
Added support to import hdf5 dataset from the Ed Boyden lab recordings for 
  paired intracellular recording validation.

## [3.1.5] -2017-11-05
Bugfix: 'jrc describe' failed to report clustering time
Bugfix: Probe file backward compatibility added
  Failed if 'maxSite' exists but 'nSites_ref' doesn't exist in the probe file (.prb)
  Now nSites_ref is computed based on 'maxDist_site_spk_um' (see wiki)
  Probe files are first searched in ./prb/ location and elsewhere if not found
Bugfix: 'jrc update' command now copies sub-folders and backup files correctly
  Previously subfolders didn't get copied and files were not back-up properly.
  **Run 'jrc update' twice to get proper updates, since you need to run a new update script.**

## [3.1.4] -2017-11-03
JRCLUST repository moved to JanalieSciComp
  https://github.com/JaneliaSciComp/JRCLUST

## [3.1.3] -2017-11-01
Bugfix: 'S0' nested saving error fixed. 
Speed improvement: Manual GUI mouse click doesn't slow down over time.

## [3.1.2] -2017-11-01
New parameter: maxDist_site_merge_um = 35 (default).
  Max. distance separation of the peak sites to consider merging based on the unit waveform similarity
Bugfix: auto-merge now computes the self-correlation of the waveform correlation matrix

## [3.1.1] -2017-10-27
Bugfix: PSTH time offset error is removed (0.01 s previously)
Manual GUI: Speed improvement when selecting a unit.
Bugfix: raw waveform display fixed for the manual GUI. 
  Waveforms did not display if you used an older .prm file
New paramter: fInterp_fet = 1 (default), set to 0 to disable this new feature since v3.0.8
  This interpolates the spike waveforms 2x and finds the optimal time-lag to project the principal vector.
New paramter: fDrift_merge = 1 (default). Set to 0 to disable this new feature since v3.0.8
  This compares mean unit waveforms at three drift locations based on the spike position (CoM).
New paramter: nInterp_merge = 2 (default), set to 1 to disable this new feature since v3.0.8
  This compares mean unit waveforms at three drift locations based on the spike position (CoM).
Bugfix: Annotations previously got deleted when deleting clusteres using "delete-auto" command.
Manual GUI: "auto-merge" to automatically merge units based on the waveform correlations.

## [3.1.0] -2017-10-26
Bugfix: raw waveform display fixed for the manual GUI. 
  Waveforms did not display if you used an older .prm file
'fir1' filter is supported for "vcFilter" parameter.
  Finite Impulse Response filter
  The frequency range is sett by "freqLim=[f_low, f_high]" parameter.
New default: Detection threshold set to "qqFactor = 5".
  Increasing the detection threshold can incorrectly merge large amplitude units with background noise.
New default: Auto-calculation of maxSite and nSites_ref based on the physical distance.
  maxSite=[] and nSites_ref=[] will automatically set it based on the site spacing and layout.
  The parameter "maxDist_site_um" is used to determine the number of sites used for features (default 50 um)
  The parameter "maxDist_site_spk_um" is used to determine the number of reference sites (default 75 um)

## [3.0.9] -2017-10-25
Faster auto-splitting. Now displaying three principal component axis.
  Individual spike waveforms are no longer displayed but mean waveforms are displayed
  nSites_ref (reference sites) are not included in the auto-splitting projection
Faster updates for manual GUI. Refresh speed increased when user selects a unit.
Bugfix: Waveform correlation matrix error fixed.

## [3.0.8] -2017-10-24
Faster manual GUI loading speed. Optimized caching
Improved drift resistance. Three average unit waveforms are computed based on the vertical spike positions,
  and waveform similarities are compared for spatial and temporal shifts.
Faster post-hoc automated merging by using parfor (parallel computing toolbox)
Bugfix: manual GUI waveform similiarity matrix calculation.
Bugfix: Removed a horizontal curve from the rho-delta plot.
Bugfix: Manual GUI will clear the unsaved results when user selects "do not save".
  Previously cached information was incorrectly shown when the manual GUI was restarted.
New default: Raw spike waveform duration is now the same as the filtered spike waveform length.
  spkLim_ms = [-.25, .75]; sets the spike waveform range to 0.25 ms before the negative peak and 0.75 ms after the peak. 
  Setting spkLim_raw_ms = []; will set spkLim_raw_ms = spkLim_ms.

## [3.0.7] -2017-10-08
Default changed: nPcPerChan=1 (previously 3), since other PC dimensions explain less than
  2% of the total variance when differentiator filter is used. Using only the first PC is 
  also more robust to the presence of overlapping spikes.
New default filter: vcFilter = 'ndiff' (convolves with [1,2,-2,-1] if nDiff_filt=2 chosen). 
  The new differentiator filter is less attenuating than 'sgdiff' for small spikes,
    which convolves with [2,1,0,-1,-2] if nDiff_filt=2 is chosen.
Default changed: fRepeat_clu=0 (previously 1)
  Repeating clustering is not necessary if lower feature dimension is used (nPcPerChan=1).
Default changed: delta1_cut=0.6 (previously 0.5) if lower feature dimension is used.
Improved auto-merging based on mean unit waveforms using interpolation
  Unit waveform similarity matrix is computed by finding the maximum correlation score 
    by relatively time-shifting waveform pairs after interpolating 4x in time.
Bugfix: Diagonal waveform similarity is calculated correctly after merge and subtract
Reduced the peak memory use by 50%. Spike waveforms are directly saved to the disk.
  fRamCache = 1; is set to default, which stores both the raw and filtered spike waveforms in RAM.
  Setting "fRamCache = 0" will only load either raw or filtered spike waveforms in RAM to save space.
Improved automated clustering for small amplitude units
  For feature extraction, vcSpkRef='nmean' now subtracts the average spike waveform across sites 
  excluding the nearest half of the sites from the center site. Previously the Local average was 
  calculated by excluding only the center site.

## [3.0.6] -2017-10-08
Improved automated clustering
  Setting the flag fRepeat_clu=1 (default) will re-runs clustering for the lower half of the cluster amplitudes.
  For Neuropixels probes, the default number of sites for feature extraction is increased:
    nSites_spk = 1+2*maxSite-nSites_ref = 1+6.5*2-4 = 10
CUDA codes are recomplied to accept up to 45 feature dimensions (previously 30).
  The max number of dimensions is stored in "nC_max" parameter (read-only, do not change).

---

[2017/9/29]
Bugfix: negative rho is fixed. This was caused by dc==0 (radius of the neighbor-finding spheres)
  when the median distance-cutoff was zero. This is fixed by discarding zero-distance values before
  calculating the quantile values.
Display fix: curves in rho-delta plot (bottom-right plot in the manual view) is suppressed.
  These spikes come from delta==0 values whose nearest pair has exactly same feature scores.
  These spikes are now hidden from the display.
"jrc gui" menu items are disabled if can't run given the current state, just like the buttons do.
Current version number of JRCLUST is overwritten in P.version after running.

[2017/9/28]
New user interface: jrc gui
  It opens a GUI that allows users to access majority of the JRCLUST commands via menu or buttons.
  This is recommended for new users who are not familiar with JRCLUST commands.
"jrc git-pull version#" will download a specific version # tag.
  "jrc git-pull" will download the latest release.
Bugfix: "version" field in .prm file is correctly updated.

[2017/9/27]
Added a parameter validater. Incorrect parameter values will be rejected with warning messages.
Bugfix: "vcFilter" issue is resolved for older format.

[2017/9/26]
jrc can run from any path after adding path to the jrclust location
  Run "addpath('jrclust_location')". Note that this doesn't persist after Matlab restart.
"jrc git-pull" updates JRCLUST from the Github repository
  https://github.com/jamesjun/JRCLUST
  git must be installed in Windows (https://git-scm.com/download/win).
Probe file moved to "./prb" folder. 
  The prm files is backward compatible. 
  .prb files will be searched in './prb/' folders if not found in the main JRCLUST folder.
Added "jrc wiki" command. 
  It opens a JRCLUST wiki webpage, a new home for the JRCLUST documentation.
Unit test is expanded.
  "jrc preview" GUI is added.
"jrc version" command is added.
  It displays the current version and the date updated.

[2017/9/17]
Bugfix: Program crashed when the number of clusters were less than 7.
  The color grouping size is 7 (# colors for cluster mean waveforms).
  Now it's possible to have any number of clusters.
Speed improvement (~40%) when fft cleanup is disabled (set fft_thresh=0).
  Previously unnecessary data conversion was performed even when fft cleanup is not used.
jrc import-nsx myrecording.ns5 myprobe.prb
  JRC imports Neuroshare format and export the analog channels to .bin file.
  JRC generates .prm file by combining the .bin and .prb file names (e.g. binfile_prbfile.prm) 

[2017/9/15]
jrc import-intan mychannels-A*.dat myprobe.prb
  JRC imports intan recordings, which saved each channel as -A###.dat files.
  JRC combines .dat files to a single .bin file and saves to the directory above.
  JRC generates .prm file by combining the .bin and .prb file names (e.g. binfile_prbfile.prm) 

[2017/9/14]
Bugfix: Waveform correlation based auto-merging
  Some cluster pairs were not merged although their correlation is above the maxCorWav threshold
jrc makeprm myfile.bin mytemplate.prm
  This create a prm file based on a template file for a given .bin file
jrc batch filelist.batch mytemplate.prm
  This creates .prm files for the list of .bin files using mytemplate.prm
  If mytemplate.prm is not specified, default.prm is used
jrc batch fileliist.batch command
  This command works for batch file containing .prm files. 
  Default command is spikesort if omited
jrc batch fileliist.batch command
  This command works for batch file containing .prm files. 
  Default command is spikesort if omited
jrc3 export-prm myparams.prm (myparam_full.prm)
  Export complete list of parameters to a new file.
  If the second argument is omitted, the current parameter file is updated.
  If some parameter values are missing, they will be copied from the default template file (default.prm).
New feature: when myparam_jrc.mat file is written, myparam_full.prm is created that contains the full parameter settings
  If the parameter file name is myparam.prm, jrc creates myparam_full.prm whenever myparam_jrc.mat is updated

[2017/9/12]
Improved rho-delta detrending
  lin-log fit using quadatic function, delta1_cut=0.6, rho_cut=-3 recommended
  Previously used log-log fit using linear function
Waiting hour-glass icon error is resolved.

[2017/8/18]
"jrc preview" UI is backward compatible with Matlab R2015b.
User interface: waiting dialog while plotting to prevent keyboard press while busy.

[2017/8/17]
Display bugfix: jrc traces spike display fixed
New user interface: jrc preview
  Change and preview threshold and sorting paramters
  Find optimal thresholds for motion rejection (blank_thresh), adaptive notch filter (fft_thresh), spike threshold (qqFactor).
  Save spike threshold per channel to a file ("vcFile_thresh = myparam_thresh.mat") and apply them to the entire recordings.  
  Auto-commit changes to the .prm file  
  Detect bad sites using max. pairwise correlation of cleaned raw traces (after notch filter and before AP-band filter).
Spike detection improved for smaller and narrower spikes
  Negative turning point detection uses different algorithm by checking left and right neighbors.
    Previously it was possible to randomly miss some of a negative turning point that has less than two neighbors below threshold.

[2017/8/7]
Bugfix: when multiple files are merged, the file ordering will follow the correct order
  by using the file creation time is stored in the .meta file (SpikeGLX generated).
    If .meta file is not availble, file creation time is used to sort the files.
"jrc dir file*.bin mylist.txt" will generate a list of binary files and 
  order them according to the file creation time.
You can specify a file sorting option at the end "jrc dir file*.bin mylist.txt [0-4]".
  0: alphabetical order
  1: meta file stored file creation time
  2: file created time
  3: file modified time
  4: numerical order in the file names   
"jrc makeprm mylist.txt myprobe.prb" will create a prm file called "mylist_myprobe.prm".
  'makeprm' now accepts .txt file containing list of binary files
  If .meta file (SpkeGLX) is available, it will be used to auto-fill the .prm file.
  If .meta file is unavailable, user will be asked to provide required info (sampling rate, # channels/file, uV/bit, probe file).

[2017/8/4]
`blank_period_ms` is now using millisecond unit. Previous it was using second unit.

[2017/8/2]
Bugfix: `makeprm` sometimes couldn't create .prm file. Now fixed.

[2017/8/1]
Bugfix: `makeprm` sometimes couldn't locate .meta file. Now fixed.

[2017/7/31]
'makeprm' command now works with a wild-card file path (csFile_merge = 'f*.bin')
  e.g. "jrc makeprm mybin*.bin myprobe.prb" creates a parameter file "mybin_myprobe.prm" by removing '*' characters and adding the probe file name.
  makeprm assumes the file format is the same for all binary files (same number of channels, sampling rate, uV/bit).
Multiple file joining (csFile_merge) now currectly orders files with numbers ('f1.bin', 'f2.bin', 'f10.bin'}
  Previously the files are ordered according to the modified timestamp. 
    The listing order can fail if the files get copied to another location, by alphabetically sorting the file order {'f1.bin', 'f10.bin', 'f2.bin'}
csFile_merge now accepts cell arrrays of wild-card
  csFile_merge = {'f1*.bin', 'f2*.bin'} is now valid.
  csFile_merge = {'f1*.bin', 'f2.bin'} is now valid.
  csFile_merge = {'f1*.bin', 'f2.bin', 'f3*.bin'} is now valid.
  csFile_merge = {'f1.bin', 'f2.bin'} is still valid.
  csFile_merge = 'f1*.bin' is still valid.

[2017/7/26]
vcFile_trial supports csv file format. 

[2017/7/25]
makeprm command can automatically determine the Neuropixels probe option from the meta file.
  If the probe is not a Neuropixels probe, user is asked to provide the probe file name.
  Call "jrc makeprm myrecording.bin myprobe.prb" to explicitly state the probe file.

[2017/7/23]
"jrc traces" UI improvement. 
  Filter on/off is indicated on the figure title bar.
  Mouse pointer changes while loading the data.
  'vcFilter_show' can override the filter setting for display purpose.
    Choices: {'', 'bandpass', 'sgdiff'}. '' uses vcFilter setting.

[2017/7/21]
Default parameter file template (default.prm) is cleaned up.
  Deprecated parameters and developer's parameters are organized separately.
SpikeGLX-generated .meta file is no longer required.
  If the meta file is missing, user is asked to provide sampling rate, uV/bit scaling factor and number of channels stored in the file.
Quality score now includes the microvolt amplitude unit.
  Minimum (S_clu.vrVmin_uv_clu) and peak-to-peak (S_clu.vrVpp_uv_clu) amplitudes from the center site is calculated.
UI: Automated deletion is implemnted. Clusters are automatically deleted based on the criteria
  Based on Unit SNR (Filtered peak amplitude divided by the RMS (Quiroga) for the primary channel
  Based on the # spike/unit.
  Undo is possible after automatic deletion.
UI: PCA feature display uses global principal vectors if vcFet='gpca' feature is used.
  Previously PCA feature display used principal vectors computed for each site.
  Site-specific principal vectors are used if vcFet='gpca' feature is not used.
Bugfix: quality metrics are now updated after manual cluster editing (delete,merge,split).
Bugfix: # of clusters displayed in "Cluster rho-delta" plot is corrected.
Cluster integrity is now automatically checked after each manual operation. 
  All cluster-related fields must have the same number of clusters.
  Manual operation is updated if and only if the automated cluster integrity check succeeds.

[2017/7/19]
New clustering feature is implemented (vcFet='gpca'). This feature is now a default feature.
  The performance of 'gpca' feature resulted in 50% lower median false positives and lower false negatives compared to the 'cov' feature.
  'nPcPerChan=3' sets number of features per channel.
There is now a filter option called (vcFilter='sgdiff'). Set this to (vcFilter='bandpass') if you want to apply a bandpass IIR filter. 
  The option below are for the Savitzky-Golay differentiator filter (vcFilter='sgdiff'), which is set by default.
    Set 'nDiff_filt=2' to set the filter tap order (+/-2 samples). Increasing this number will smooth the signal and decreasing will sharpen the signal.
  The options below are for the bandpass IIR filter (vcFilter='bandpass').
    Set 'freqLim=[300 6000]' to apply a bandpass filter from 300-6000 Hz. 
    Set 'freqLim=[300 inf]' to apply a high-pass filter with a 300 Hz cut-off.
    Set 'freqLim=[0 300]' to apply a low-pass filter with a 300 Hz cut-off.
    Set 'fEllip=1' to apply an Elliptic filter. Set it to 0 if you want to use a Butterworth filter.
    Set 'filtOrder=3' to apply a third-order bandpass filter.
'vcDetrend_postclu=global' option results in a more reliable detrending for the cluster center finding (rho-delta plot). 
  Default parameter values are set to delta1_cut=0.6, rho_cut=-3.
  delta1_cut indicate log10(z-score). For example, delta1_cut=0.6 corresponds to z_thresh = 10^0.6 = 4.0 after detrending.

[2017/7/14]
"jrc batch my_prm_files.batch" command now lists list of errors found during the batch process
dc_percent behavior changed. Now dc_percent is applied to local neighbors sharing the same
  primary site only. Previously neighbors sharing primary or secondary sites were included,
  leading to overestimation of the distance-cutoff.
  dc_percent default value is set back to 2.
"vcFilter" paramter is added to clarify whether differentiation filter or bandpass filter is used.
  vcFilter='sgdiff' uses nDiff_filt value to run the Savitzky-Golay differentiator.
  vcFilter='bandpass' uses freqLim value to run a band-pass filter.
    Butterworth filter is used if "fEllip=0", and otherwise elliptic filter is used.
  'filtOrder=3' sets the filter order to 3.
Drift view shows spikes from a specific shank. Modify iShank_show=1 in the prm file (shows the first shank by default).


[2017/7/13]
Bugfix: L12677 index range error is fixed
 viKeep_spk = find(viTime_spk >= ilim_spk(1) & ilim_spk <= ilim_spk(2));

[2017/7/12]
Matlab path is temporarily reset to the system default when JRCLUST starts 
  to prevent namespace collision error. The path is reset to user-specified 
  value after restarting the Matlab session.

[2017/7/7]
Waveform merging now checks other clusters whose primary or secondary sites 
  overlap with the cluster's secondary site (second center).
default value of dc_percent changed (2->1). This results in lower false positives without affecting false negatives.
fft_thresh GPU error is corrected. 

[2017/7/6]
Quality metrics added: ISI ratio, L-ratio, Isolation Distance
  Variable names: "S_clu.{vrIsiRatio_clu, vrLRatio_clu, vrIsoDist_clu}"
  Local neighborhood is selected by searching for all spikes sharing the same 
    primary or secondary sites with the cluster center.
Bugfix: average waveform is recalculated after removing outlier spikes.
Bugfix: Multi-shank neighboring sites are calculated correctly, now always taken from the same shank.
ISI violation is automatically fixed by removing spikes that are double-detected.
"jrc manual": deleted spikes are now hidden.

[2017/7/5]
Neuropixels probe layout corrected. Channel 1 is in the middle column but previously it was on the edge
  (imec3.prb, imec3_opt1, imec3_opt2 imec3_opt3 imec3_opt4 corrected)

[2017/7/4]
Default parameter changed to improve spike detection and clustering.
  qqThresh = 5; nDiff_filt = 2;
Improved post-cluster merging based on the location of the primary and secondary peaks.
  Previously the merging decision was considered for cluster pairs within 50 um.
  The new scheme yields more specific cluster merging, thur reducing false cluster merges.
Removal of outlier spikes based on the Mahalnobis distance.
  thresh_mad_clu = 7.5; (default). Set to 0 to disable this.
  This will significantly reduce the false positive spikes with minimal increase of missed spikes.

[2017/6/29]
Matched filter detection is supported (set vcFilter_detect = 'matched';)
  This convolves the spike wave shape to the filtered traces. 

[2017/6/27]
Bugfix: JRC3 now supports energy and PCA features.
Significantly more reliable spike detection at the edges of memory loading blocks. 
  Extra spike detected due to the edge effect is ~0.1 ppm per channel per memory load.
Hard threshold setting is fixed (spkThresh_uV).

[2017/6/23]
Bugfix: "traces" command multiple time segment viewing now adds time separators
  jrc2 and jrc3 updated and time offset values are shown.

[2017/6/22]
"traces" command now supports multiple time segments (nTime_traces parameter).
  nTime_traces = 1 shows a single time segment as previously.
GPU is disabled for manual curation (jrc manual)

[2017/6/19]
Code dependencies to external functions removed in jrc2 and jrc3.
Run "jrc3 dependencies" to show dependencies to toolboxes and other functions.

[2017/6/15]
Paged load uses overlapping time to reduce missed spikes at the edges due to filter edge effect.
  "nPad_filt" controls number of samples for overlap.
  This feature is added to both jrc2 and jrc3.

[2017/6/14]
JRCLUST version 3 is released. Run "jrc3" or "jrc" to use this version.
  Substantial improvement with dataset with probe drift. 
    Set "nTime_clu" to 1 if no drift occurs, 2 or higher if drift occurs fast
  "jrc3 detect" runs spike detection and feature extraction.
  "jrc3 sort" runs clustering and posthoc merging.
  "jrc3 auto" runs posthoc merging.
  "jrc3 manual, spikesort-manual, sort-manual, or auto-manual" runs manual UI.
Better global noise rejection by blanking noisy period.
  "blank_thresh" sets the MAD threshold for the channel mean-based activity rejection
  "blank_period_ms" sets the duration of the blanking. Useful to remove licking or motion artifacts  


[2017/6/6]
"jrc2 auto" command is added. This performs re-clustering after updating the
  post-clustering parameters such as "delta1_cut" (cluster confidence).
  This must be done after running "jrc2 sort" or "jrc2 spikesort".
  This command uses the features previously computed.
Local CAR now uses "nSites_ref" channels with the least activities (lowest stdev).
  Previously reference sites were selected by choosing four furthest sites from the center.
If vcSpkRef="nmean" is used, the features from the local reference sites are set to zero.
  This excludes features from the sites with low SNR.

[2017/6/4]
Robust GPU memory processing. If GPU memory error occurs during spike detection the 
  processing is performed in CPU after transfered the data to the main memory 

[2017/5/29]
Default setting changed to produce lower false negative under drift.
  Now checks 2x greater spatial range. This doesn't hurt performance for non-drifting dataset
Post-cluster waveform-based merging became 2x faster.
  Now checks 2x greater spatial range for merging. 
  A cluster pair must have at least 4 overlapping sites

[2017/5/24]
Document updated. Output files and variables described. 
  Type "jrc2 doc" to view the PDF manual.
Improved history logging behavior
  All manual operations are stored in "cS_log" variable.
  Manual operations are reset after automated sorting.
  Ten last performed operations are stored in RAM. 
  The number of logged entries can be changed in "MAX_LOG" parameter.
  Lastest program state is automatically saved to the disk for each manual operation.

[2017/5/12]
Bugfix: Spike detection error fixed when refrac_factor >1 was used
  The previous version assigned an incorrect site # for the spike being detected.
  This is a serious bug and update is recomended.
Spike merging speed is improved.

[2017/5/6]
Bugfix: Multishank post-hoc merging error fixed
Added jrc2 batch-mat myfile_batch.mat command
 myfile_batch.mat contains csFile_bin (binary files) and csFile_template (template files)
 csFile_template can be a string or cell strings
 Optionally provide csFile_prb (probe files), which can be a string or cell strings (for each .bin file)

[2017/5/5]
History menu is added to GUI. This feature allows you to undo/redo manual actions,
  and restore to previous state if JRCLUST crashes. Restart JRCLUST and select 
  the last action saved in the "History" menu (left of Help menu).

[2017/5/3]
Bugfix: New update (4/26) can display older format saved.
Faster post-hoc merging
Faster merging and splitting
Improved visualization for the cluster correlation matrix (FigWavCor)

[2017/4/26]
'nneigh_min_detect' parameter is introduced. This sets the minimum number of neighbors around spikes
  below the negative detection threshold. Set between [0,1,2]. 1 is recommended.
  Previous version used hard-coded value of 2 which lead to missing spikes.

[2017/4/25]
'nDiff_filt' parameter is intorduced to control the smoothing parameter of the differentiator filter.   
  Set nDiff_filt between 1-4. 3 is recommended. Higher number results in more agressive smoothing. 
  Previous version used 4 but it smoothed too aggressively leading to loss of small and narrow spike waveforms
  To skip differentiator sgfilter, set nDiff_filt = []; 
'tlim_load' parameter is now enabled, which sets the time range of data to load (in sec unit).
  For multiple file merging, the time base refers to each individual file and it applies for all files being loaded

[2017/4/21]
Bugfix: Split selection refinement fixed (polygon adjustment didn't take effect)
Bugfix: Raw waveform display works with the dataset sorted by the previous version 
Spike position view (bottom-left) shows 4x more spikes
Color can be changed for projection view by changing 'mrColor_proj' in the .prm file
  mrColor_proj = [r0 g0 b0; r1 g1 b1; r2 g2 b2]; % 0: background; 1: Cluster selected; 2: Cluster compared

[2017/4/20]
Multiple feature viewing enabled. Go to menu and select 'Projection' to change features to display
  available features: {'vpp', 'cov', 'pca'}. This affects the projection and time figures.
  To change a default feature to display, set 'vcFet_show' to one of {'vpp', 'cov', 'pca'} in the .prm file.

[2017/4/19]
Raw waveform is enabled for manual clustering. set 'fWav_raw_show=1' to enable raw waveform view.

[2017/4/5]
Instruction manual updated (install.txt). Updated installation from the website (jrc2.zip).
Added feature export command. 
  Run "jrc2 export-fet" and this outputs mrFet, miFet_sites to the workspace.

[2017/4/4]
Trace view: For multi-file sorting, you can select files from the command line 
    instead of from the dialog box, which sometimes run out of the monitor space.
    Alternatively, you can directly select a file to display by running
    "jrc2 traces myparam.prm File#"
Bugfix: Projection view polygon cutting resulted in empty view when diagonal 
    channel views are selected. Now fixed.

[2017/4/2]
Bug-fix: cluster restore error resolved in "jrc2 manual".
Speed-up: post-cluster step is 2x faster thanks to parallelization.

[2017/03/31]
Added feature: exporting spike amplitudes in uV unit, organized by clusters
    "jrc2 export-spkamp myparam.prm (clu#)"
    This outputs Vpp, Vmin, Vmax to the workspace (cmrVpp_clu, cmrVmin_clu, cmrVmax_clu)

[2017/03/29]
Compiled for GPU with >4GB memory (nvcc -m 64 -arch sm_35)
    Supporting NVIDIA GPU with compute capatiliby 3.5 or later
    Kepler, Maxwell and Pascal architecture (circa 2013, GeForce 700 series or later)
UI time view debugged. Channel switching behavior fixed.

[2017/03/28]
UI error fixed: Waveform error after splitting/merging
Kilosort compilation script
"jrc2 export-spkwav" command added. This exports spike waveforms from all clusters.
    run "jrc2 export-spkwav unit#" to export and plot spikes from a specified unit#.

[2017/03/26]
Significant improvement in UI startup and cluster selection speed.
Waveform display fixed. Previously displayed differentiated waveform.

[2017/03/24]
Fixed a slow startup for jrc2 manual

[2017/03/23]
Faster plot in "jrc2 manual". Number of plot objects in FigWav reduced to 10.
Export to version 1 format supported via "jrc2 export-jrc1 myparam.prm" command.
Menu commands are included in the unit test.
CUDA compile script supported (jrc2 compile).

[2017/03/22]
Waveform fix after splitting in the FigPos window (left-bottom).
"jrc2 probe" remembers the current prm file
Split in the projection view, polygon adjustment now editable
Kilosort integrated. run "jrc2 kilosort myparam.prm", and run "jrc2 manual" for manual merging/splitting


[2017/03/17]
Export-csv feature fixed
Clusters are indicated as randomized but persistent color and line width for each spike

[2017/03/14]
Code clean-up. jrc2.m dependency was reduced to ~5 external functions.

[2017/03/12]
UI bug fix (jrc2 manual)
    UI automated testing implemented
    Merge and split errors fixed
UI feature improvements
    Amplitude zoom change decoupled beteween waveform, time and projection views
Multiple file processing supported (set csFiles_merge = {'file1', 'file2, ...} or csFiles_merge = myfiles*.bin).
Multiple shanks supported (edit "shank" field in .prb file, see sample.prb file for example).
    "jrc2 probe" command displays different shanks in unique colors.
    

[2017/02/24]
JRCLUST ver 2 update (run jrc2)
    Unit test added
    Memory-efficient sorting (jrc2)
    Improved trace view (jrc2 traces)

[2017/01/31]
sample.prb fixed
Image processing toolbox requirement specified.

[2017/01/23]
File merging error fixed (LFP files are now merged for Neuropixels probe).
Selective site loading by using a custom-prb file

[2016/12/13]
Added IMEC sync export function. It outputs vnSync variable that contains 16-bit uint16.

[2016/10/30]
Fixed GUI command "File>export waveform". Now correct waveform is saved.

[2016/10/25]
Default parameters changed (delta1_cut, dc_percent) to improve the sorting outcome, 
    now splits clusters less.

[2016/10/20]
jrc makeprm supports merging multiple files for multi-shank probe. 
    You need to create a probe file for each shank to use this feature.
    Copy the probe file(myprobe_shank1.prb) to the jrclust folder and run
    "jrclust makeprm 'c:\directory\data_1*.bin' myprobe_shank1.prb
    This will merge the .bin files to data_1all.bin file and create 
    'data_1all_CNT_4x16probe_shank1.prm' in where the .bin file resides.
    you can then run "jrclust spikesort" to sort the data
You don't need to specify .prm file after you specify once. for example,
    After running "jrclust makeprm myparam.prm", you can simply run
    "jrclust spikesort" and it will automatically use "myparam.prm".
You can use "jrc" instead of "jrclust".
Template file is now supported. Try
    "jrclust makeprm rawdata.bin myprobe.prb my_template.prm"

[2016/10/11]
jrc manual debugged. The waveforms are now scaled and displayed correctly. Previously the software showed differentiated waveform when the differentiation filter was enabled.

[2016/10/9]
Added FFT-cleanup procedure which removes narrow-band noise. set 'fft_thresh=0' to disable. fft_thresh=10 recommended (z-score threshold).
More efficient handle on the raw waveform (S0.mrWav moved to global variable).

[2016/10/5]
More efficient local common average referencing
Improved feature used by default (spacetime), which uses physical position of spikes to cluster.

[2016/9/29]
Default settings updated
Spacetime feature uses minimum amplitude

[2016/9/27]
Waveform error fixed after split and merge in UI.
Faster waveform calculation.

[2016/9/22]
Fixed waveform shape. Now it shows non-differentiated waveform when 'spacetime' feature is used.
Faster load for raw traces.
Improved time alignment

[2016/9/20]
GPU errror during filtering fixed.
Default settings changed.

[2016/8/30]
Projection split error fixed.
Added "auto-correlation" feature. (vcFet = 'acor')

[2016/8/28]
Background noise cancellation (set fMeanSite=2 for local median subtraction)
Added 'moment' feature. This projects spikes to mean, sqrt(mean), mean.^2

[2016/8/21]
2x faster file loading using GPU for combining loading/processing operation.
2x faster and more accurate spike detection.
Currently process 5x realtime speed for 120 channels, 1.2x realtime for 374 channels

[2016/8/18]
2x faster filtering and spike detection using multiple CPU cores.

[2016/8/17]
"jrclust maketrial" fixed for a version older thatn R2015a.

[2016/8/15]
"jrclust auto" algorithmic update
    Graded rho/delta calculation
    Determine clustered distance distribution
    Use gradual drop-off spatial windowing instead of shart drop-off and increaes the range.
Improved automated clustering (Waveforms are now centered before projecting principal components).
Faster feature computation for manual curation.

[2016/8/14]
Documentation updated. Show the latest documentation by running ("JRCLUST doc").
Fixed "jrclust install" and "jrclust update" commands.
"jrclust makeprm" now supports multi-file merging using a wild-card. 
    Run "jrclust makeprm myfiles_*.bin myprobe.prb" and JRCLUST will create 
    "myfiles_all.bin", "myfiles_all.meta" and "myfiles_all.prm"

[2016/8/12]
Added auto split/merging when 'spca' feature is used ("postCluster.m").
Added an option to use amplitude-scaled dc (distance-cutoff) parameter for R-L clustering.
    Try dc_frac = 1/3; and use (vcCluDist = 'neucldist';). This sets a radius parameter to 
    compute "rho" parameter in R-L clustering for each individual spikes. Small spikes will
    use smaller radius and large spikes will use larger radius scaled by the norm of vector.
    Previously I used the same dc_frac parameter for all spike amplitudes which resulted in
    over-merging for smaller spikes and over-splitting for larger spikes.
Added "jrclust compile" for recompiling the CUDA GPU codes.
Features are not saved but computed on-fly. This saves lots of memory and time.
default.prm updated. Try the new default configuration which allow much 
    1. more sensitive spike detection
    2. much faster and memory-efficient feature handling
    3. more accurate automated sorting.

[2016/8/5]
sPCA feature added (single-channel feature and clustering).
prm loading error fixed.
"jrclust install" script debugged

[2016/8/3]
Debuged spike display of raw, unfiltered waveform.
Debugged traces view, transition from 'a' to 'f'.

[2016/8/2]
Added "Save figures" under File menu ("jrclust manual").
More sensitive spike detection. Set vcSpatialFilter = "mean" and fPcaDetect = 1.
Time view can display second and third principal components ("jrclust manual").
sPCA feature added. Cluster by channel.

[2016/8/1]
Fixed "jrclust makeprm" command
Fixed "jrclust download" command

[2016/7/28]
Reference can be a channel (vcCommonRef = "1" will subtract channel 1 from the rest of channels)

[2016/7/27]
Cut in PCA view (Press 'F' in amplitude projection view or Time view);

[2016/7/15]
Dipole-based unit shift (in time view, press "e")
Unit track and shift feature (in time view, press "c")

[2016/7/22]
Debugged jrclust loadparam behavior (prb file did not load)
Faster Vpp computation and mean waveform bcomputation
Clusters are initially sorted by their centroid locations
A cluster or a pair of clusters can be selected by clicking on the traces or correlation matrix.
"jrclust manual" layout is rearranged

[2016/7/21]
Export waveform in the time window (jrclust manual)
plot traces debugged.
Changed behavior of "jrclust makeprm myrecording.bin myprobe.prb".
    The prb file can contain settings to be copied to myparam.prm.
"jrclust trackset mycollection.set" plots collection of prm files. 
    Refer to "ph3o3.set" for an example.

[2016/7/19]
Debugged trace view zoom and startup ("jrclust traces").
Large time view by skipping spikes every n samples (nSkip_show).
Faster spike display ("jrclust traces").

[2016/7/18]
Feature plot enabled. PCA or features other than Vpp can be plotted.
Speed improvement in correlation matrix computation.
Speed imporvement in projection view.
X-position view added in the time plot (press 'x', 'y', 't' to switch between x,y position and time views).
Common zoom behavior between the waveform and time view (press up/down arrows).
Click on he Wavform correlation view to select clusters.
Waveform view now permits clicking on the traces. Previously this could not select waveforms.
Local common reference debugged (set vcCommonRef='mean' and nSite_ref>0.
Faster disk loading and filtering (~x2 speed-up).
"jrclust cabletest myrecording.bin"  added for Neuropixels phase 3 probe (BIST #7?).

[2016/7/17]
Added track depth menu in "Time vs. Amplitude" view. Press "T" to switch to the "track depth" view.

[2016/7/16]
Debugged jrclust ui, performance improvement for projection window

[2016/7/15]
Cluster exporting menu added in "jrclust manual".
"jrclust download sample" command added. This downloads sample.bin and sample.meta files from Neuropixels Phase 2 probe (128 sites).

[2016/7/13]
vpp feature debug fixed.
Now supports "jrclust makeprm my_recording.bin my_probe.prb"
Debug fix for "jrclust traces" command. Now supports aux. channel view.

[2016/7/12]
Faster trace view interface. Example: "jrclust traces my_param.prm"
UI supports dynamic parameter file updating. No longer need to restart the UI to update prm file.
Edit command added. Example: "jrclust edit myparam.prm".

[2016/7/10]
Much more responsive UI. Significantly faster to update graphics
Cluster Position calculation
Cluster annotation ability

[2016/7/7]
Faster local common average referencing using GPU memory caching (x10 faster).

[2016/06/28]
LFP-based depth tracking. Parameters have "_track" postfix in .prm file.

[2016/06/21]
Local mean subtraction syntax changed (vcCommonRef='mean', nChans_
Goto time added in "traces" view. Press 'G' to go to the start time to plot.
Bad channels being excluded from the computation all together

[2016/06/20]
Memory and speed optimization for

[2016/06/19]
Memory optimization for common mean reference
Local mean subtraction is supported (vcCommonRef = 'localmean')

[2016/06/17]
Can work with files that were stopped abruptly during recording (multiple of nChans requirement relaxed).

[2016/06/10]
+GUI: Press 's' to autosplit in the waveform view
+ added a feature to reject spikes during movement. Movement is detected FROM the global mean across channel (common reference). Motion is defined when the absolute value of the mean across sites exceed 5 SD. 
    Set "thersh_mean_rejectSpk=5" under "detection" section to change the threshold. Set this to 0 to skip motion rejection.
+File saving speed-up. Feature tensor (trFet) is saved in a binary format (.fet).
+PCA speed-up and reduction in memory usage.

[2016/06/09]
Decreased memory use by half. Loading a fraction of a bin file (~1.2 GB) at a time and transposes.

[2016/06/07]
LFP and Aux channel handling updated
Phase 3 probe updated

[2016/05/21]
run "jrclust update" to update from dropbox

[2016/05/18]
Otsu's graythresh method for automatically determining dc_cut. This is recommend instead of distribution based sorting (set vcDc_clu = 'Otsu').
Channel query when viewing traces (jrclust traces). Press 'c' and draw a rectangle to find out a channel number. (useful in LFP view mode).
New clustering feature (vcFet = 'diff248' see http://www.synopsys.com/community/universityprogram/pages/neuralspikesort.aspx). 
    This works very well when combined with Otsu auto-thresholding method.

[2016/05/08]
tlim_clu fixed when fDetrend_postclu is enabled


[2016/04/28]
memory use optimization. load one channel at a time.


--------------------------
[ToDo]
Drift compensation
Choose spike center according to Vpp
projection view using features used for spike-sorting
auto-split keyboard shortcut and make it faster
memory optimization by loading bit of waveform at a time
annotate clusters
warning message for overwrite
LFP phase analysis integration for jrclust_ui
 place cell integration
cut in rho-delta view
click on a matrix view to show tw pairs
quality score display
store metafile in P structure. save as raw text
realtime spike sorting
manual command startup speedup

[ToDo, Tahl Holtzman]
Waveform width measurement
Toggle between positive and negative clusters (give access to deleted clusters)
relate spike color in the trace view and cluster number

[Todo, Joao, 2016 06 17]
resample in xy projection 'a' 
self similarity in time by dividing into different epocs. 
spike sorting mean waveform output
display the split in the features space

[Anupan, 2016 07 15]
jrclust manual: debug when split
export text field to csv file

---------------------------------------------

@TODO: Before clustering determine viSite_spk based on x,y centroid (nearest site)
@TODO: Depth shift
@TODO: cut in the position view
@TODO: Trace view spike in lfp view and hiding spikes
@TODO: Fix sitestat
@TODO: Auto-split and merging after running clusering
@TODO: R-L cluster: add max. time neighborhood, fix the neighbor size. 
@TODO: extreme-drift correction. Track global drift and shift.
