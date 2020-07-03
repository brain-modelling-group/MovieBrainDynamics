# Scripts for analysis of rest and movie MRI data with Hidden Markov Modeling


## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
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

<a name="demo"></a>
# Demo

## Functions

See [all-on-one-script.m]() for the complete analysis that will produce analyze the data and produce the figures.

The main fuctions in the scripts directory are:
- Step1
- Step2
- Step3
- 


<a href="citation"></a>
# Citation

For usage of the scripts and the associated manuscript, please use the following:

Johan N. van der Meer, Michael Breakspear, Luke J. Chang, Luca Cocchi (2020). Movie viewing elicits rich and reliable brain state dynamics. Nature Comms (under review).

