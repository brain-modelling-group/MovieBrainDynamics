# Scripts for analysis of rest and movie MRI data with Hidden Markov Modeling


## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [User Guide](#user-guide)
- [License](./LICENSE)
- [Citation](#citation)


# Overview

This repository contains the files and Matlab scripts needed to reproduce the Hidden Markov Model results in the paper "Movie viewing elicits rich and reliable brain state dynamics".

This work examines brain state dynamics while watching a movie and compares them to pre-movie resting-state. This movie itself can be found on youtube: https://www.youtube.com/watch?v=p98KAEif3bI, and on its Wikipedia page: https://en.wikipedia.org/wiki/The_Butterfly_Circus.

Functional imaging was performed with a 3 Tesla Siemens scanner to obtain one whole-brain BOLD-weighted image every 0.82 seconds, for a total of 535 functional images (20 min) throughout the entire movie. Prior to the movie a resting-state scan of 225 functional images (8 min) was performed. This was repeated (within-subject) after about 3 months. Sessions are code-named resta, restb, mov1a and mov1b. This repository provides extracted timeseries data and matlab scripts to re-run the analysis to reproduce all the figures in the manuscript.

<a name="repo-contents"></a>
# Repo Contents

- [data](./data): extracted timeseries; pupil diameter; heart rate; story annotations
- [scripts](./scripts): scripts to run the HMM analysis and produce figures on brain states; state paths; state dynamics; state associations
- [results_10](./results_10): empty directory that will be filled with (intermittent) analysis results and figures.
- [figures](./figures): empty directory that will be filled with the figures that can be also found in the main manuscript. Some of these figures were later assembled or modified in Inkscape for presentation in the main manuscript; however the raw data is not changed.
- [neurosynth](./neurosynth): folder containing all python requirements for the production of Figure 3. This is a separate type of analysis from the Matlab-based HMM-MAR analyses.
- [source-data-files](./source-data-files) source data file underlying each figure (saved in excel format)


<a name="system-requirements"></a>
# System Requirements

## Hardware Requirements

This will run on laptop/desktop with:

RAM: 16+ GB  
CPU: 4+ cores

The runtimes below are generated on a desktop computer (16 GB RAM, 4 cores@3.3 GHz)

## Software Requirements

The Matlab scripts have been tested on on Matlab 2018a. It has not been tested for any other Matlab versions.

The "Toolbox for segmentation and characterisation of transient connectivity" (**HMM-MAR**) toolbox is already incorporated (see [all-on-one-commented-script.m]()) in the [scripts](./scripts) folder. The commit number of this HMM-MAR version is: [7a5915c](https://github.com/OHBA-analysis/HMM-MAR/commit/7a5915c8efb899ac3860a21a1e29ffd9b995d6f6). Newer, up-to-date versions can be downloaded from its Github page: [https://github.com/OHBA-analysis/HMM-MAR](https://github.com/OHBA-analysis/HMM-MAR). Check especially its [wiki](https://github.com/OHBA-analysis/HMM-MAR/wiki) page for more information. The folder is [hmm_new](./hmm_new), which contain some additional scripts to save HMM results.

The scripts of this repository are basic Matlab scripts and should run on any OS on any computer. They have been tested on a Linux (Ubuntu 18.04) operating system. 
- You also need the "Statistics and Machine Learning Toolbox". Matlab is not able to do standard PCA analyses (which the HMM-MAR requres) out-of-the-box. Approximate cost: 600 AUD.

The Neurosynth analysis runs as a python jupyter lab notebook. You need to install conda/python on your computer to run these analyses. See below for installation instructions.

<a name="installation-guide"></a>
# Installation Guide

## Matlab code
- git-clone this repository
- download and install [Matlab version 2018a](https://www.mathworks.com/products/new_products/release2018a.html)
- download the SPM scripts: [SPM version 12](https://www.fil.ion.ucl.ac.uk/spm/)
- unzip extracted_timeseries.zip in [data](./scripts/extracted_timeseries.zip) folder
- unzip rois.zip in [data](./data/rois.zip) folder
- Start Matlab and navigate to the github repository
- add SPM scripts folder to the path
- open, view & run the [all-on-one-script.m]()


## Python Neurosynth code
- download and install conda version 4.5.4 from the [miniconda archives](https://repo.anaconda.com/archive/)
- go into the [neurosynth](./neurosynth) folder with a terminal.
- create the 'neurosynth' environment with conda from the provided .yml file: `conda env create -f environment.yml`
- the neurosynth module is already downloaded from the [Neurosynth github repository](https://github.com/neurosynth/neurosynth/commits/master), the commit number was: [b795589](https://github.com/neurosynth/neurosynth/commit/b795589ef42839124dd55a59eca45d510d2487a9)
- start up jupyter lab: `jupyter lab`
- open up the file "neurosynth_analysis.ipynb"
- Follow instructions in this notebook; it will make Figure 3 from the main manuscript.

Installation time is mainly determined by copying and unzipping files and downloading the required components; Given the number of steps it should take approximately 2-3 hours.

<a name="user-guide"></a>
# User Guide

## Functions

1. unzip in the data folder the extracted timeseries
2. unzip in the data folder the rois
3. start Matlab
4. go to the 'scripts' folder
5. Run in sequence the following scripts:

|Matlab Script|Function|
|---|---|
|`Step1_create_dirs_and_run_hmm.m`| Makes all the directories and run the HMM; outputs HMMrun_rep_X and SummaryStatistics_X .mat files in `results_10` folder. The number of inferences is set to 15 repetitions.|
|`Step1b_Check_for_expressions.m`| Checks for each inference if all the brain states are expressed at least *once* during the entire data|
|`Make_Fig1_create_brain_state_figure.m`| Makes pictures of how the states are observed - inside the `results_10` folder|
|`Make_Fig2_state_path_figure.m` | Makes pictures of the viterbi state path for all data, and how these states are similar across subjects|
|`Step4a_make_subspec_matrices_from_statepath.m` | Pre-analysis for dynamic assessment - makes `outr.mat` file containing intermediate data|
|`Make_Fig4_OccDwell.m` | Makes Fractional Occupancy and Dwell Time figure|
|`Step6a_calculate_NBS_mov_rest_sesa.m` | Pre-analysis that runs NBS (Zalesky method) on the state transition matrices - intermediate data are saved in .mat files with the same name as the script - move > rest analysis|
|`Step6c_calculate_NBS_mov_rest_sesb.m` | Pre-analysis which is the same as before, but for the rest > movie contrast|
|`Make_Fig5_Dynamics.m` | Makes the state transition/dynamics figure|
|`Make_Fig7_FO.m` | Runs analysis of IS-RSA comparing dynamics and questionnaire answers - Dynamics are assessed with F0|
|`Make_Fig7_transitions.m` | Runs analysis of IS-RSA comparing dynamics and questionnaire answers - Dynamics are assessed with State Transition matrices (nonzero elements) |
|`f11_radio_annots_electrophys_now_the_annots_II.m` | Matches HMM brain state paths to the movie Annotations; both of which are 0-1 vectorized. See Manuscript Figure 8|
|`f12_new_analysis.m` | Checks if the HR during brain state visit X is higher/lower than the mean HR (across time) for all subjects|
|`f12_new_analysis_PD.m` | Same script as above, but modified to find if the PD brain state visit X is higher/lower|
    
6. For Neurosynth analysis of the brain states - see Neurosynth folder.

*comment: The naming of the scripts and how the functions have been written is rather organic; in future work we would use Python Notebooks to make for cleaner coding practices*

<a href="citation"></a>
# Citation

For usage of the scripts and the associated manuscript, please use the following:

Meer, J.N. van der, Breakspear, M., Chang, L.J., Sonkusare, S., Cocchi, L., 2020. Movie viewing elicits rich and reliable brain state dynamics. Nat. Commun. 11, 5004. [https://doi.org/10.1038/s41467-020-18717-w](https://doi.org/10.1038/s41467-020-18717-w)
