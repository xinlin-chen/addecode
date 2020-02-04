# addecode
Bass Connections project on Exercise Therapy and Alzheimer's Disease

This repository contains the code for reading in MRI k-space data from the Bruker machine and reconstructing images from undersampled k-space data. Additional toolboxes are also needed and can be found in BassData/MRI Toolboxes on Box. Previously extracted k-space data is also available in BassData, saved as brukerdata.mat.

The main script in this repository is run_bart.m, which walks you through reading in data, generating a sampling pattern, and reconstructing images from undersampled data, all using the Berkeley Advanced Reconstruction Toolbox. You may elect to use our data or use data from the espirit workshop, in which case the BART folder at https://github.com/mikgroup/espirit-matlab-examples must be downloaded.
