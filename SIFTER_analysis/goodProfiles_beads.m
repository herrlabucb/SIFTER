function [data_struct] = goodProfiles_beads(data_struct,r2_threshold,num_peaks)
% Perform quality control on intensity profiles, removing lanes with SNR<3
% and allowing the user to select lanes to remove upon visual inspection
% Generate intensity profiles from the ROI stacks in the output of roiGeneration
% and perform background subtraction on the profiles
% 
% Outputs
% Struct [structure]: A data structure containing objects (Intensity 
%                     profiles for each ROI, 3D matrix with each ROI 
%                     contained in a different z, and coordinates of 
%                     the ROIs).
%   struct.rois: 3D matrix with each ROI contained in a different z.
%   struct.angle: The angle of rotation to straighten the image (number, in degrees).
%   struct.rotate: The angle of rotation required to display the image with
%   separations running vertically instead of horizontally (number, in
%   degrees).
%   struct.array_bounds: User selected boundaries of the array as a 3x2
%   matrix (rows contain upper left, upper right, and lower left
%   coordinates respectively; first column contains x-coordinates; second column contains y-coordinates). 
%   struct.name: The name of the protein target entered by the user
%   (string).
%   struct.wells_per_row: The number of wells per row based on the user
%   selected array bounds and horizontal well spacing.
%   struct.rows: Number of rows in the array
%   struct.int_prof: a 3D matrix containing an array of the intensity 
%                    profiles [x, intensity value] indexed by the 
%                    third dimension
%   struct.good_devices: boolean vector indicating the good devices
%   struct.good_indices: vector of indices of the intensity profiles and
%                        ROIs that were fit with Gaussian peaks 
%                        (i.e. passed the SNR threshold)
%   struc.fit_coefficients: m x 3 x p matrix containing the gaussian fit
%                           coefficients. m is the number of peaks and p                
%                           is the number of good devices (good_devices)                
%   struct.R2: m x 1 matrix containing the R^2 valuesof the Gaussian fits.
%              m is the numberof good devices (good_devices)
%   struct.dev_to_analyze: a bitmask of the good devices The indices refer 
%                          to the lanes that were fit (struct.good_indices)
%   struct.index_dev_to_analyze: a vector of the indices to the devices  
%                                that passed QC. The indices refer to the 
%                                lanes that were fit (struct.good_indices)
%   
%
% Inputs
%   Struct [structure]: The data structure from fitPeaks containing the
%                       intensity profiles and fit parameters for each lane
%      

% v05 (8.18.17) Gaussian fits are now displayed with the intensity profiles
%               during quality control. Also, a "zoom off" button was added
%               so the user may use zoom features in the GUI, and press the
%               "zoom off" button to re-activate the curve selection/next
%               button.

%% Check input arguments
switch nargin
    
    % If only the data_structure is provided, set r2_threshold = 0.7
    case 1
        
        r2_threshold = 0.7;
        
    % If provided, ensure the r2 value is valid
    case 2
        % Exit function if an invalid r2 value is input
        if ((r2_threshold<0) || (r2_threshold > 1))
            
            error('Invalid R^2 value');
            
            return
            
        end
end       

% Get the intensity profiles
int_prof_all = data_struct.int_prof;

% Determine the number and size of the intensity profiles
[x_dim, y_dim, z_dim] = size(int_prof_all);

% Get the R^2 values from fitPeaks and find all of the lanes that exceed
% the minimum R^2 value 
r2 = data_struct.R2;

fit_coeffs = data_struct.fit_coefficients;

if num_peaks == 1
    good_r2 = find(r2 >= r2_threshold);
else
    
    for i = 1:length(r2)
        r2_comparison(i) = min(r2(:,i));
    end
    good_r2 = find(r2_comparison >= r2_threshold);
end
% determine array position of high r2 value fit lanes
good_indices = data_struct.good_indices;
good_fits = good_indices(good_r2);

good_int_profiles = zeros(x_dim,y_dim,length(good_r2));

% Slice out the lanes that meet the R^2 value threshold
for i = 1:length(good_r2)
    good_int_profiles(:,:,i) = int_prof_all(:, :, good_fits(i));
    good_coeffs(:,:,i) = fit_coeffs(:,:,good_r2(i));
end

% set number of rows/columns of subplots to display in each figure window
n=5;
num_subplots = n * n;

% Calculate the number of plots that need to be displayed and the number of
% n x n panels that are required
plots_display=length(good_r2);
number_subplots = ceil(plots_display/(n*n));

% Preallocate a bitmask to indicate the lanes that pass QC
good_devices = ones(length(good_r2), 1);

% Preallocate a vector containing the indices of devices to analyze
dev_to_analyze=zeros(z_dim, 1);

% Preallocate a vector for the intensity profiles that pass QC
good_subplots = ones(plots_display,1);

disp(number_subplots);


% for loop to generate subplots for user inspection of the intensity profiles
for i=1:number_subplots
    
    disp(i);
    
    h = figure;
    if i==1
        devices_subplot = (1:(n*n));
        
    elseif i*n*n <= plots_display
        devices_subplot=((i*n*n)-(n*n)+1):((i*n*n));
    
    else 
        
        devices_subplot=((i*n*n)-(n*n)):(plots_display);
    end
    
    % Plot all of the intensity profiles
	for j=1:length(devices_subplot)
              
            % Get the index for the intensity profile subplot
            dev_number = devices_subplot(j);
            
            % Get the intensity profiles
            device = good_int_profiles(:,:,dev_number);
            xval=device(:,1);
            yval=device(:,2);
            
            device_coeffs = good_coeffs(:,:,dev_number);
            % Make the subplot
            subplot(n, n, j);
            
            % Plot. The buttondownfucn Toggles the subplot selection.
            % Function at the bottom of the script
            plot(xval,yval,'LineWidth',2,'Tag', sprintf('%d', dev_number),...
                'buttondownfcn', @clickTest);   
                    
         
            hold on
            for k = 1:num_peaks
                coeff = device_coeffs(k,:);
                a = coeff(1);
                b = coeff(2);
                c = coeff(3);
                yGauss = a*exp(-((xval-b)/c).^2);
                plot(xval, yGauss,'r','LineWidth',2);
                hold on
            end
                
            
        
    end
    next=0;
    
    % Add a button to go to the next set of intensity profiles. Function 
    % at the bottom of the script
    btn = uicontrol('Style', 'pushbutton', 'String', 'Next',...
        'Position', [500 15 50 30],...
        'Callback',@continueButton2);
    
    zoom_off = 0;
    btn = uicontrol('Style', 'pushbutton', 'String', 'Zoom Off',...
    'Position', [250 15 50 30],...
    'Callback',@zoomOff);

       hManager = uigetmodemanager(h);
            try
                set(hManager.WindowListenerHandles, 'Enable', 'off');  % HG1
            catch
                [hManager.WindowListenerHandles.Enabled] = deal(false);  % HG2
            end
            set(h, 'WindowKeyPressFcn', []);
            set(h, 'buttondownfcn', []);
            set(h, 'buttondownfcn', @clickTest);
            set(h, 'KeyPressFcn', @myKeyPressCallback);
              
 
     % Wait for the user to be done selecting lanes before moving to the
    % next panel
    while next == 0   

        pause(0.01);
    end
    
    %good_devices(devices_subplot(1):devices_subplot(end))=good_subplots;
    
% Close the current window before creating the next one 
close(gcf)    
end 

% Save the a bitmask of the good devices The indices refer to
% the lanes that were fit (struct.good_indices) 
data_struct.dev_to_analyze = good_subplots;

% Save a vector of the indices to the good devices. The indices refer to
% the lanes that were fit (struct.good_indices)
good_subplot_ind = find(good_subplots == 1);
data_struct.index_dev_to_analyze = good_r2(good_subplot_ind);
end

function [exitZoom] = zoomOff(zoom_off,event)
zoom_off = 1;
if zoom_off ==1
    zoom off
end
end


function [next]=continueButton2(qstring,title,str1,str2,default)
% This function is called with the button the QC GUI is pressed. It
% determines if the is satisified with their selection and moves to the
% next set of intensity profiles if they are.

% Create a dialog to ask the user if they are done.
qstring = 'Are you done selecting devices to throw out?';
title = 'Device Quality Control';
str1 = 'Yes';
str2 = 'No';
default ='Yes';
choice = questdlg(qstring,title,str1,str2,default);
    
    % If the user is done, move to the next set of intensity profiles. If
    % they are not, allow them to continue selecting intensity profiles    
        switch choice
            case 'Yes';
                disp([choice 'Great, let''s keep going then!'])
                next=1;
            case 'No';
                disp([choice 'Okay, please finish selecting devices to throw out'])
                next=0;
        end
                  
% Set the "next" variable to true to allow the while loop to terminate and
% the next set of intensity profiles to be displayed
assignin('caller', 'next', next);

end


function clickTest(line_handle, event)
  % This function handles toggling the subplot state when it is clicked

  % Get the vector of subplot states from the calling workspace.
  good_subplots = evalin('caller', 'good_subplots');
  
  % Determine the subplot number of the selected plot.
  current_tag = get(line_handle, 'Tag');
  
  % Get the current state of the selected intensity profile
  subplot_num = str2num(current_tag);
  subplot_state = good_subplots(subplot_num);
  
  % Uncomment the line below to display which subplot is selected and its
  % current state (for debugging)
  %disp(sprintf('%d, %d', subplot_num, subplot_state));

  % Toggle the selection based on the last character in the tag 
  % (0 = off, 1 = on)
  
  if (subplot_state)
      
     set(line_handle, 'Color', [1, 0, 0]);
 
     good_subplots(subplot_num) = 0;
     
  else
      
     set(line_handle, 'Color', [0, 0, 1]);

     
     good_subplots(subplot_num) = 1;
     
      
  end
  
%  disp(good_subplots);
  assignin('caller', 'good_subplots', good_subplots);
  
    
end


