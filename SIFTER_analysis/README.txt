# SIFTER_analysis
SIFTER_analysis is MATLAB software intended for the analysis of single-cell protein complex electrophoretic separations collected with the SIFTER assay. This software is modified from summit (https://github.com/herrlabucb/summit, under the GPL-3.0 License).
# System Requirements
Software dependencies: requires MATLAB version 9.1 (2016b) or higher and the Curve Fitting Toolbox (version 3.5.4 or higher) and Image Processing Toolbox (version 9.5 and higher).
Version tested on: MATLAB 2019b
No required non-standard hardware

# Installation Guide
Instructions: 
1.	Open MATLAB.
2.	In MATLAB use the “Current Folder” navigator to go to the folder in which you want to download SIFTER_analysis.
3.	Right click within the “Current Folder” white space and click on **Source Control > Manage Files**.
4.	In the window titled “Manage Files Using Source Control”, insert https://github.com/herrlabucb/SIFTER_analysis.git into “Repository Path”.
5.	Click “Retrieve” to download the SIFTER_analysis source code. The “SIFTER_analysis” and “demo” folders should now appear in your “Current Folder”.
6.	Finally, add the software to your path by right clicking on the “SIFTER_analysis” folder within “Current Folder” and clicking **Add to Path > Selected Folders and Subfolders**.
Typical installation time is < 1 minute.

# Demo
To try out SIFTER_analysis, please see the script, data and example results in /demo/run_demo.m. Expected run time for the demo analysis is < 30 minutes.

# Instructions for Use
The “Demo” section includes an example of how to use SIFTER_analysis on an image of a SIFTER device. To extract protein expression data from the MATLAB data structure when not running F_G_ratio_SNR (as in the demo), the AUCs for lanes that passed quality control and SNR thresholding are accessible with the following:
struct_ind_dev = struct.index_dev_to_analyze;
struct_snr_ind = struct.inds_good_snr;
struct_snr_good_inds = struct_ind_dev(struct_snr_ind);
 
struct_AUC = struct.AUC;
struct_AUC = struct_AUC(:);
struct_AUC_good = struct_AUC(struct_snr_good_inds);

When analyzing multiple images of the same gel, first align the images in Fiji using ‘Landmark Correspondences”. Run the SIFTER_analysis functions to obtain the data structure for one image. For the remaining images of the gel, when calling roiGeneration2, include the structure name as the final function input. This will apply the same ROIs that were generated for the first analyzed image to the other gel images.
