function [struct] = roiGeneration2(filename,horzspacing,vertspacing,bidirectional,struct)
%This function rotates/aligns the raw image and segments the image into ROIs
%   Outputs:  
%Struct [structure]: A data structure containing objects (3D matrix
%   with each ROI contained in a different z, and coordinates of the ROIs)
%   Inputs:
%filename [string]: A string containing the name of the fluorescence image
%                   to be processed
%horzspacing [num]: Well-to-well spacing (horizontal)
%vertspacing [num]: Well-to-well spacing (vertical)
%bidirectional [num]: Input the number "1" if the separation is
%bidirectional. This will offset the ROIs so the microwell is in the
%center of the lane.

%%versions
% 0.1-Created April, 2016
% 0.2-(2017) added bidirectional option
% 0.3 (3.3.20): Updated to apply same transform for ROI generation if user
% inputs a struct with the fields "angle" and "rotate".
%% Check input arguments
switch nargin
    % If the user only provides the image, horizontal and vertical spacing
    case 4
        transform = 0;
    case 5
        transform = 1;
        
        tf = isstruct(struct);
            if tf == 0
                 error('Input argument "struct" is not a structure.');
            
            return
            
            end
            
        % retrieve previously determined angle for transformation of image
        angle = struct.angle;
        
        % retrieve previously determined array boundaries
        array_bounds = struct.array_bounds;
        
        % extract the individual x an y coordinates of the array boundaries
        x_upperleftwell = array_bounds(1, 1);
        y_upperleftwell = array_bounds(1, 2);
        
        x_upperrightwell = array_bounds(2, 1);
        y_upperrightwell = array_bounds(2, 2);
        
        x_lowerrightwell = array_bounds(3, 1);
        y_lowerrightwell = array_bounds(3, 2);
        
    
    otherwise
        
        error('Invalid number of input arguments');
            
        return
        
end

%Load the image file in Matlab
img=imread(filename);
    if transform == 0
%Display more contrasted img in window
contrasted_img=histeq(img);
imshow(contrasted_img);
title('Take a look at the array and determine if the wells are oriented left of the bands or right of the bands. Then press any key');
pause()
% Construct a questdlg to ask the user how the image is currently oriented
% for coarse rotation
    choice = questdlg('Are the wells currently left of the bands or right of the bands?', ...
	'Current array orientation', ...
	'Wells are left of bands','Wells are right of bands','Wells are right of bands');
        % Handle response
        switch choice
            case 'Wells are left of bands';
                disp([choice 'Okay, the image will be rotated to the right!'])
                rotate = -90;
            case 'Wells are right of bands';
                disp([choice 'Okay, the image will be rotated to the left!'])
                rotate = 90;
        end
   % Store the course rotation angle to orient the array vertically to
        % the struct
        
        struct.rotate = rotate;
    else
         
        rotate = struct.rotate;
    end
    
  imgrotated=imrotate(img,rotate);
  contrasted_img_r=histeq(imgrotated);
  imshow(contrasted_img_r);
  
    %If struct was not an input argument (and there is no previous
  %angle/array boundary values to draw from), the user will now manually
  %select the array boundaries.
  while transform == 0
    test=1;
    while test==1
    %Prompt user to select the upper right well of the array. 
    title('Please zoom in on the the middle of the upper left well and press any key.');
    zoom on;   % use mouse button to zoom in or out
    pause()
    zoom off;
    
 % preallocate array bounds matrix
    array_bounds = zeros(3, 2);
    
    title('Please click on the middle of the upper left well.');
    [x_click,y_click]=ginput(1);
    x_upperleftwell=x_click;
    y_upperleftwell=y_click;
    zoom out;
    
     array_bounds(1,:) = [x_upperleftwell, y_upperleftwell];

    %Prompt user to select the upper right well 
    title('Please zoom in on the middle of the upper right well and press any key.')
    zoom on;   % use mouse button to zoom in or out
    pause()
    zoom off;
    title('Please click on the middle of the upper right well.');
    [x_click,y_click]=ginput(1);
    x_upperrightwell=x_click;
    y_upperrightwell=y_click;
    zoom out;
    
    array_bounds(2,:) = [x_upperrightwell, y_upperrightwell];
    
    
    %Prompt user to select the lower right well of the array 
    title('Please zoom in on the middle of the lower right well and press any key.')
    zoom on;   % use mouse button to zoom in or out
    pause()
    zoom off;
    title('Please click on the middle of the lower right well.');
    [x_click,y_click]=ginput(1);
    x_lowerrightwell=x_click;
    y_lowerrightwell=y_click;
     
    array_bounds(3,:) = [x_lowerrightwell, y_lowerrightwell];
    
    % store all of the coordinates of the array bounds to the struct
     
    struct.array_bounds = array_bounds;
       
    % Construct a questdlg to ask the user if they are happy with their
       % well selection
    choice = questdlg('Are you happy with your well selections?', ...
	'Well selections for array boundaries', ...
	'Yes','No','Yes');
        % Handle response
        switch choice
            case 'Yes';
                disp([choice 'Great, let''s keep going then!'])
                test = 0;
            case 'No';
                disp([choice 'That''s okay, try again!'])
                test = 1;
        end
    if (x_upperrightwell<x_upperleftwell || y_upperrightwell>y_lowerrightwell)        
        test=1;
        title('Oh no! We detected you selected the wells in the wrong order. Please try again. Press any key to continue')
        pause()
    else
        test=0;
    end
end
%store the coordinates of the direction vector that extends from the upper left well to the right most point of the array
dir_vector1=[x_upperrightwell,y_upperleftwell]-[x_upperleftwell,y_upperleftwell];
%store the coordinates of the direction vector that extends from the upper left well to the upper right well 
dir_vector2=[x_upperrightwell,y_upperrightwell]-[x_upperleftwell,y_upperleftwell];
%Find angle between the two direction vectors [angle in degrees]
cosangle=dot(dir_vector1,dir_vector2)/(norm(dir_vector1)*norm(dir_vector2));
angle=acosd(cosangle);
    if (y_upperrightwell<y_upperleftwell)
        angle=-angle;
    end
    
% store the angle used to straigten the image in the struct
    struct.angle=angle;  
    transform=1;
  end
  %Rotate the image so the array is aligned
b=imrotate(imgrotated,angle,'nearest','crop');
b_contrasted=histeq(b);
imshow(b_contrasted);
hold on
sz=size(b)/2;
rotation_matrix=[cosd(-angle),-sind(-angle);sind(-angle),cosd(-angle)];
new_upper_left=rotation_matrix*[(x_upperleftwell-(sz(2)));(y_upperleftwell-sz(1))];
new_upper_right=rotation_matrix*[(x_upperrightwell-sz(2));(y_upperrightwell-sz(1))];
new_lower_right=rotation_matrix*[(x_lowerrightwell-sz(2));(y_lowerrightwell-sz(1))];
x_new_upper_left=new_upper_left(1)+sz(2);
y_new_upper_left=new_upper_left(2)+sz(1);
x_new_upper_right=new_upper_right(1)+sz(2);
y_new_upper_right=new_upper_right(2)+sz(1);
x_new_lower_right=new_lower_right(1)+sz(2);
y_new_lower_right=new_lower_right(2)+sz(1);
%generate matrix to store ROIs
%Determine number of wells per row
wells_per_row=round((x_new_upper_right-x_new_upper_left)/horzspacing);
%Determine number of rows
rows=round((y_new_lower_right-y_new_upper_right)/vertspacing);
%Determine total number of wells
total_wells=wells_per_row*rows;
%for loop to fill in the 3D matrix with ROIs from the image (proceeds row by row of the microwell array from left to right)
mat=zeros(vertspacing,horzspacing,total_wells);

%the start location of the column array will be offset for bidirectional
%separations (so the well is centered in the ROI)
if bidirectional==1
    well_offset=vertspacing/2;
else
    well_offset=0;
end 

for i=1:rows
    for j=1:wells_per_row
        z=(wells_per_row)*(i-1)+j;
        row_start=(round(x_new_upper_left)-horzspacing/2)+((j-1)*horzspacing);
        row_end=row_start+horzspacing;
        col_start=(round(y_new_upper_left)-well_offset+((i-1)*vertspacing));
        col_end=col_start+vertspacing;
        x=row_start:1:(row_end-1);
        y=repmat(col_start,1,length(x));
        y2=col_start:1:(col_end-1);
        x2=repmat((row_end-1),1,length(y2));
        mat(:,:,z)=b(col_start:(col_end-1),row_start:(row_end-1));
        plot(x',y','Color','w','LineStyle','-');
        plot(x',y','Color','k','LineStyle',':');
        plot(x2',y2','Color','w','LineStyle','-');
        plot(x2',y2','Color','k','LineStyle',':');
    end
end
struct.rois=mat;
end

