function [struct] = F_G_ratio_SNR(struct,SNR_threshold)
%This function calculates the fraction of F-actin binding protein bound from a bidirectional separation.
%The function assumes that the F actin is the first row and G actin is the
%second row of the AUC array.
%v2 (3.19.21): Added optional input argument SNR_threshold which is a
%logical (1 = SNR threshold applied). With this input argument, the indices
%for lanes that passed an SNR threshold of 3 (from the SNR calculation
%function snr_calc.m). Note: these indices correspond with the indices in struct.index_dev_to_analyze. 

%% Check input arguments
switch nargin
    
    % If only the data_structure is provided, set good_inds =
    % struct.index_dev_to_analyze
    case 1
        
        good_inds = struct.index_dev_to_analyze;
        
    % If stuct.inds_good_snr provided, good_inds = struct.inds_good_snr
    case 2
        if SNR_threshold==1;
            inds = struct.index_dev_to_analyze;
            snr_inds =struct.inds_good_snr;
            
            good_inds = inds(snr_inds);
        else
            good_inds = struct.index_dev_to_analyze;
        end
    
    otherwise
        
        error('Invalid number of input arguments');
            
        return
    
    
end

auc_values = struct.AUC;

good_auc = auc_values(:,:,good_inds);

f = good_auc(1,:,:);
g = good_auc(2,:,:);
f_actin = f(:);
g_actin = g(:);


for i = 1:length(f_actin)
    f_ratio(i) = f_actin(i)/(f_actin(i)+g_actin(i));
end 
struct.f=f;
struct.g=g;
struct.f_ratio=f_ratio;

end

