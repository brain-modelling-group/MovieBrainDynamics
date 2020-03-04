# Brain Dynamics during Movie Viewing assesed with HMM-MAR Generative Model

This repository contains the files needed to reproduce the Hidden Markov Model results in the paper **Movie viewing elicits rich and reliable brain state dynamics**

This work examines brain state dynamics while watching a movie. This movie can be found on youtube: https://www.youtube.com/watch?v=p98KAEif3bI, and on its Wikipedia page: https://en.wikipedia.org/wiki/The_Butterfly_Circus. 

Functional imaging was performed with a 3 Tesla Siemens scanner to obtain one whole-brain BOLD-weighted image every 0.82 seconds, for a total of 535 functional images (~20 min) throughout the entire movie. Prior to the movie a resting-state scan of 225 functional images (~8 min) was performed. This repository provides extracted timeseries data and matlab scripts to re-run the analysis to reproduce Figure 1 (Brain States) and Figure 2 (Brain State Dynamics).

This repository consist of:
 - data folder containing movie annotations and 2 zip files. One zip file contains the extracted timeseries. These series have names formatted as: "aroma-ts-<COND>-s<PP>-r<ROI>.txt", where COND is resta, restb, mov1a, mob1b; PP is participant number (1-20); ROI is the network. The other zip file contains the ROI network masks.
 - Matlab scripts for preprocessing. The three main scripts are: 
   - Step1_create_dirs_and_run_hmm.m
   - Step2_create_brain_state_figure.m
   - Step3_create_state_path_figure.m


The directory structure is as follows:

```console
├── data
│   ├── annotations
│   ├── extracted_timeseries.zip --> extracted timeseries
|   └── rois.zip --> rois
|
├── results_10 --> output folder of Step1 script
│   ├── aroma
│   └── aroma_paper --> the aroma results used in the original paper
+
└── scripts
    └── hmm_new
        └── HMM-MAR-master --> a clone of the HMM-MAR repository version 
```

 
 
Compared to resting state, with a more bi-stable configuration, movie viewing is characterized with richer dynamics where:
 - 5-6 brain states are visited
 - brain states are more differentiated in their BOLD activation across networks
 - Participants are more consistent w.r.t. each other in terms of the visisted state

