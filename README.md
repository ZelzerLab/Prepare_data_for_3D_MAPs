# Prepare_data_for_3D_MAPs
![image](https://user-images.githubusercontent.com/58815992/112968847-5e94f480-9155-11eb-98d1-dcd0a43269b1.png)


These scripts extract morphometric measurements from 3D binary images of cells and nuclei and prepares them as MATLAB files for [3D MAPs](https://github.com/ZelzerLab/3D-MAPs/blob/main/README.md/). 

These codes were written by Tomer Stern and Ankit Agrawal. To learn about 3D MAPs, read the [preprint](https://www.biorxiv.org/content/10.1101/2020.07.28.225409v1/).

## System Requirements

Linux, Windows and Mac OS are supported for running the code. These scripts were written and compiled with MATLAB version R2016b and depend on MATLAB toolbox "Statistics and Machine Learning Toolbox". They require at least 16GB RAM (depending on number of objects in each image). For the data size analyzed in [3D MAPs](https://github.com/ZelzerLab/3D-MAPs/blob/main/README.md/), 256 GB RAM was needed. It was tested primarily on a workstation with Windows 10 and 512 GB RAM.

## Installation

This process should take less than 10 minutes.

1. Install MATLAB version 2016b or later.
2. Download the files from this repository. You can do this by cloning the repository or by downloading it as a ZIP folder.
3. Unzip the folder and run. An example dataset can be found in the supplementary section of the manuscript.



## Creating (Surfaces).mat files

This script converts 3D binary tif images into a (Surfaces).mat file containing a triangulated mesh of each object. 

1. Open the script "prepare_segmentation_file_for_analysis" by double clicking it.

2. Set the number of parallel workers to use (using more workers will speed up runtime).

![image](https://user-images.githubusercontent.com/58815992/112824484-4066c000-9093-11eb-9998-77b9b31bf9e5.png)

3. Set the minimum and maximum volume range, to eliminate noisy objects of abnormal size.

![image](https://user-images.githubusercontent.com/58815992/112824756-96d3fe80-9093-11eb-8ddc-2ba7aaa047ff.png)

4. Set the voxel size of the image in microns [x,y,z].

![image](https://user-images.githubusercontent.com/58815992/112825156-18c42780-9094-11eb-9506-9789a5886726.png)

5. Create a variable with the data path containing all of the images to process by typing the following in the command window



 ```
path = 'F:\data'
```

6. Run the script from the command window by typing the following
 ```
prepare_segmentation_file_for_analysis(path)
```

## Merging multiple (Surfaces).mat files
This script merges multiple (Surfaces).mat files into one. This should only be performed if there are multiple surfaces files for the same imaging coordinates. i.e. two segmentation thresholds were used on the same image. Make sure to copy the images to be merged into a separate folder, as all files in the folder will be merged. 


1. Open the script "mergingSurfaceFile" by double clicking it.

2. Set the data path.

![image](https://user-images.githubusercontent.com/58815992/112827256-c0daf000-9096-11eb-8e57-147f8635f496.png)

3. Run the script from the command window by typing the following
 ```
mergingSurfaceFile
```

4. Replace the multiple (Surfaces).mat files with the single merged file.

## Creating (Characteristics).mat files
This script reads (Surfaces).mat files and extracts morphometrics on each object. If cells and nuclei are given as an input it will store the nuclei data in G.nuc, cell data in G.cell, and the cells which have a matching nucleus in G.inter. Only cells with a matching nucleus less than the volume of the cell will be used for cell morphometric analysis in 3D MAPs, i.e. G.inter. 

1. Set the voxel overlap threshold for matching nuclei and cells by opening the script "edu_add_calc_morphological_characteristics_remove_overlap" and entering the desired threshold value. 

![image](https://user-images.githubusercontent.com/58815992/112834331-870ee700-90a0-11eb-8f2d-e2f9d90053de.png)

2. Open the script "main_remove_overlap" and set the data path. Nuclei and cell (Surfaces).mat files should be in the same directory.


 ```
data_dir_name = 'F:\data';
```
3. If you are running only nuclei then set 

```
ignore_cells = true;
ignore_edu = true;
```

If you are running both cells and nuclei then set


```
ignore_cells = false;
ignore_edu = true;
```

*Notes: The EDU field is a third binary field that can be activated. For example, whether the cells or nuclei have an EDU signal. This was not implemented in the manuscript.
 
4. If there is a known shift between the two imaging channels in Z from the microscope than set the Z shift as a function of the nuclei slices. Otherwise, leave as 0.

![image](https://user-images.githubusercontent.com/58815992/112835803-6f386280-90a2-11eb-925b-1f58224a8cdd.png)

5.  Run the script from the command window by typing the following and follow the instructions in the command window
 ```
main_remove_overlap
```

*Notes: ellipsoid_fit_new.mat is modified from Yury (2021). Ellipsoid fit (https://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit), MATLAB Central File Exchange. Retrieved April 1, 2021.

## License
[Apache-2.0](https://opensource.org/licenses/Apache-2.0)
