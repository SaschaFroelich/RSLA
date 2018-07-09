function [ROISCOURSE, INTENS, numberofrois] = coursematrix(BOLDfile,mask)
% Creates ROISCOURSE matrix by averaging the BOLD timeseries for one ROI
% over all involved voxels.
%
% Input:
%   BOLDfile:   .nii-file complete with path.
%   mask:       .nii-file specifying ROIs. Voxels with identical intensity
%               values are considered as one ROI. Voxels with intensity
%               value 0 are ignored.
%
% Output:
%   ROISCOURSE: ROIsx(timesteps+1) matrix. Each row vector in ROISCOURSE is
%               the time-series of a ROI. NB: First column in ROISCOURSE
%               contains intensity values of ROIs. Timecourses follow in
%               columns 2:end.
%   INTENS:     A vector containing the intensity values of all ROIs
%               (except 0).
%   numberofrois:   no. of ROIs.

% read BOLD-data from nifti-file
[data,~,~,~,~] = readnifti(['' BOLDfile '']);
info = nii_read_header(['' BOLDfile '']);

% read mask 
[data2,~,~,~,~] = readnifti(['' mask '']);
info2 = nii_read_header(['' mask '']);

% Check if data is being reshaped in readnifti.m or not
if size(data,2) == 1
    error('Must reshape data!');
end

% no. of voxels per image. Every voxel will be uniquely identified by its
% position in the array
voxels = size(data,1)*size(data,2)*size(data,3);

% Check if same no. of voxels
if voxels ~= size(data2,1)*size(data2,2)*size(data2,3)
    error('Number of voxels in imaging data does not correspond to no. of voxels in mask file!');
end

% reshape data to VOXELSxTIME matrix COURSES
% For the voxels in data, which is an axbxcxd matrix,
% data(i,j,k,l) = data(m), where m=(l-1)*abc+(k-1)*ab+(j-1)*a+i
COURSES = reshape(data,[voxels,size(data,4)]);

ROIS=reshape(data2,[voxels,1]);

%un stores the unique values of ROIS (in mask file) in ascending order. The first value has to be 0 (which must not be a ROI in the brain)!
un = unique(ROIS);
unlen = length(un);

%check if lowest value in ROIS is zero
if un(1) ~= 0
    msgbox('The lowest intensity in the mask file is different from zero.','Error','error');
end

%Create Matrix of size ROISxTIME (unlen-1 because first entry of un is 0 which does not correspond to a ROI)
ROISCOURSE = zeros(unlen-1,size(data,4));

%The vector INTENS stores the intensitiy values of the rois in the mask corresponding
%to the timecourses in ROISCOURSE. 
INTENS = zeros(unlen-1,1);
INTENS = un(2:unlen);

%Next, average the timecourses of all voxels inside one ROI
%start from 2 as unlen(1) is zero and does not correspond to a ROI!
for i=2:unlen
    %The elements of I describe the position where in the array ROIS an
    %element has the value of un(i)
    I = find(ROIS == un(i));
    %len is equal to the number of voxels in that specific ROI
    len = length(I);
    %for every timestep, calculate the average BOLD signal value
    for j=1:size(data,4) %iterates over time steps
        for k=1:len %iterates over all voxels in ROI
            %Calculate the sum over all ROIS in COURSES at that time point
            ROISCOURSE(i-1,j) =  ROISCOURSE(i-1,j) + COURSES(I(k),j);
        end
        %Divide the sum by the number of voxels in that ROI to obtain the
        %average
        ROISCOURSE(i-1,j) = ROISCOURSE(i-1,j)/len;
    end
end

%Now ROISCOURSE(i,:) is the average time course of all the voxels in the
%ROI that has the intensity value un(i+1)!

numberofrois = unlen-1;

end