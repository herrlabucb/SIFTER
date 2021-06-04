% Generate the ROIs. When directed to indicate if wells are right or left
% of the bands, choose 'right' for this example image, as the microwells are to
% the right of the F-actin bands in the example. This will orient the
% F-actin peaks to the right of the microwell in the intensity profile.
[test_struct] = roiGeneration2('2019-07-24_H1_MDAMB231_gfp-actin_F_left_G_right_635_smooth.tif',80,200,1);
pause

% Generate 1D intensity profile for each ROI
[test_struct]=intProf(test_struct,5);
close all;

% Fit the intensity profiles to a Gaussian function. Note: When asked to
% click to select peak boundaries, select the boundaries for F-actin first.
[test_struct]=fitPeaks_beads(test_struct,2,3);

% Perform quality control, eliminating peaks with punctate noise in the
% vicinity of the peak, fit performed on "empty lane" with punctate noise, etc.
% The area-under-the-curve (AUC) of peaks is calculated here to quantify
% protein expression.
[test_struct]=goodProfiles_beads(test_struct,0.7,2);

% Identify the peaks with SNR > 3
[test_struct] = snr_calc(test_struct,2,5);

% Calculate the F-actin ratio from F- and G-actin AUCs (F/F+G). This function
% assumes that the F-actin peak boundaries were selected first in
% fitPeaks_beads.
[test_struct] = F_G_ratio_SNR(test_struct,1);

save('test_struct.mat','test_struct');
