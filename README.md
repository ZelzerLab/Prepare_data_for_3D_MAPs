# Prepare_data_for_3D_MAPs
These scripts extract morphometric measurements from 3D binary images of cells and nuclei and prepares them as MATLAB files for [3D MAPs](https://github.com/ZelzerLab/3D-MAPs/blob/main/README.md/). 

These codes were written by Tomer Stern and Ankit Agrawal. To learn about 3D MAPs, read the [preprint](https://www.biorxiv.org/content/10.1101/2020.07.28.225409v1/).

## System Requirements

Linux, Windows and Mac OS are supported for running the code. These scripts were written and compiled with MATLAB version R2016b and depend on MATLAB toolbox "Statistics and Machine Learning Toolbox". They require at least 4GB RAM, but was tested primarily on a workstation with Windows 10 and 512 GB RAM. 

## Installation

This process should take less than 10 minutes.

1. Install MATLAB version 2016b or later.
2. Download the files from this repository. You can do this by cloning the repository or by downloading it as a ZIP folder.
3. Unzip the folder and run.



## Creating (Surfaces).mat files

This script converts 3D binary tif images into a (Surfaces).mat file containing a triangulated mesh of each object. 
