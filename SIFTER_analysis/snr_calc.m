function [struct] = snr_calc(struct,num_peaks,backgroundwidth)
%This function calculates the signal-to-noise ratio (SNR) of each peak from lanes that passed QC from goodProfiles.
%The signal is taken as the Gaussian fit peak amplitude. The noise is
% the standard deviation of a background region (the left-hand 'gutter' area of the
%lane). The background region is limited to pixels in the gutter within two peak widths of
%the peak center.


% Inputs
%   Struct [structure]: The data structure from goodProfiles containing the
%                       intensity profiles and fit parameters for each lane
%   num_peaks [num]: The number of peaks in the separation lane.
%   backgroundwidth [num]: The width of the gutter region (in pixels) to
%   be used for the background. 

% Outputs
%   Struct [structure]: The data structure with added variables:
%   struct.signal [matrix]: An m x p matrix, where m is the number of peaks
%                           (num_peaks), and p is the number of lanes that passed QC in
%                           goodProfiles (i.e., the length of index_dev_to_analyze). Each value in
%                           the matrix is the peak amplitude from the Gaussian fit in fitPeaks.
%   struct.noise [matrix]: An m x p matrix, where m is the number of peaks
%                          (num_peaks), and p is the number of lanes that passed QC in
%                          goodProfiles (i.e., the length of index_dev_to_analyze). Each value in
%                          the matrix is the standard deviation of the background region.
%   struct.snr [matrix]: An m x p matrix, where m is the number of peaks
%                       (num_peaks), and p is the number of lanes that passed QC in
%                       goodProfiles (i.e., the length of index_dev_to_analyze). Each value in
%                       the the signal divided by the noise (SNR).
%   struct.inds_good_snr {matrix]: An m x k matrix, where m is the number of peaks (num_peaks)
%                                  and k is the number of lanes for which
%                                  each peak has an SNR > 3. Note: if all
%                                  lanes had SNR > 3, this matrix will be
%                                  the same length as index_dev_to_analyze.
% v.01: Function created (3.25.21)
good_inds = struct.index_dev_to_analyze;

array_inds = struct.good_indices;

map_inds = array_inds(good_inds);

coeffs = struct.fit_coefficients;

good_coeffs = coeffs(:,:,good_inds);

mat=struct.rois;

intensity_profiles = struct.int_prof;


for i = 1:length(good_inds)
    ind = map_inds(i);
     % Get the image of the lane
    lane = mat(:,:,ind);
    % Get one background ('gutter') region at the left edge of the lane
    left_backgroundregion = lane(:, (1:backgroundwidth)); 
    background_region = [left_backgroundregion];
    x = intensity_profiles(:,1, ind);
      
    for peak = 1:num_peaks
        lane_coeff = good_coeffs(peak,:,i);
        % Get peak center and width to find indices for peak region 
        center=lane_coeff(2);
        sigma=lane_coeff(3);
        width=sigma/sqrt(2);
        
        %Get peak amplitude from fit
        amp = lane_coeff(1);
        
        %determine location of +/- 2 peak widths from the peak center
        auc_left_bound=center-2*width;
        auc_right_bound=center+2*width;
        
        % Determine index of auc_left_bound and auc_right_bound for selection of background values in
        %the region of the peak
        left_diff_auc=abs(x-auc_left_bound);
        left_data_auc=find(left_diff_auc==min(left_diff_auc));
    
        right_diff_auc=abs(x-auc_right_bound);
        right_data_auc=find(right_diff_auc==min(right_diff_auc));
        
        % Make sure the left bound is within the array
        if (left_data_auc < 1)
           
           left_data_auc = 1; 
            
        end
        
        
        % Check to make sure the AUC bounds are within the bounds of the
        % array
        if (right_data_auc > length(x))
            
            right_data_auc = length(x);
            
            
        end
        
        peak_bg = background_region(left_data_auc:right_data_auc,:);
        %Calculate the noise for the lane as the standard deviation of the
        %background gutter region
        noise(peak,i) = std(peak_bg(:));
        %Calculate the SNR by dividing the peak amplitude by the noise
        snr(peak,i) = amp/noise(peak,i);
        signal(peak,i) = amp;
end 


end
struct.snr = snr;
struct.signal = signal;
struct.noise = noise;

good_snr = snr > 3;
%Identify indices of the lanes for which SNR was calculated that have peaks with SNR > 3
if num_peaks == 1
    inds_good_snr = find(good_snr==1);
else
    inds_good_snr = find(sum(good_snr)==num_peaks);

end

struct.inds_good_snr = inds_good_snr;