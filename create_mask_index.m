% create_mask_index(input1,input2)
%   INPUT:
%   input1:             .nii-file, complete with path
%   input2 (optional):  min_vox: Ignore all ROIs with min_vox voxels or less
%
% Script creates .tst-file and .mat-file in same directory as input-file
function []=create_mask_index(varargin)

mask = varargin{1};

if nargin==2
   min_vox=varargin{2}; 
else
   min_vox=0;
end

filename=mask(1:end-4) %nifti-files

[data2,~,~,~,~] = readnifti(['' mask '']);
info2 = nii_read_header(['' mask '']);
SIZE2 = size(data2);
voxels = SIZE2(1)*SIZE2(2)*SIZE2(3);
ROIS=reshape(data2,[voxels,1]);

% un stores the unique values of ROIS (i.e. the intensity values in mask file) in ascending order. The first value has to be 0 (which must not be a ROI in the brain)!
un = unique(ROIS);
unlen = length(un);

% For the voxels in data, which is an axbxc matrix,
% data(i,j,k) = data(m), where m=(k-1)*ab+(j-1)*a+i

rois=[];
fiD=fopen(['' filename '.txt'],'w');
for iter=2:unlen
    % The elements of I describe the position where in the array ROIS an
    % element has the value of un(i).
    I = find(ROIS == un(iter));

    vox_in_ROI=length(I);
    
    %use only first entry of I as exemplatory. I(1)=m=(k-1)*ab+(j-1)*a+i.
    %Determine i,j,k.
    
    %determine how many columns (there are SIZE2(1) rows per column)
    %col=I(floor(length(I)/2))/SIZE2(1);
    col=I(1)/SIZE2(1);
    %determine how many pages that corresponds to (there are SIZE2(2) columns per page)
    pag=ceil(col/SIZE2(2));
    
    %It is thus in which column on that page?
    col = col-(pag-1)*SIZE2(2);
    
    %Thus it is on the pag page in column ceil(col). Thus, the row is
    row = I(1)-((pag-1)*SIZE2(1)*SIZE2(2)+(ceil(col)-1)*SIZE2(1));
    
    %These coordinates are correct in matlab, but but necessarily correct in fsl, as fsl
    %starts with 0!
    i = row;
    j = ceil(col);
    k = pag;
    
    intensity=data2(i,j,k);
    
    roi = iter-1;

    if vox_in_ROI > min_vox
        rois=[rois, intensity];
        fprintf(fiD,['Roi no. ' num2str(roi) ': i,j,k (MATLAB)=' num2str(i) ',' num2str(j) ',' num2str(k) '. Intensity=' num2str(intensity) ', #voxels: ' num2str(vox_in_ROI) ' \n']);
    end
   
end

intensities_in_mask=un(2:end);
save(['' filename '_atlas'],'intensities_in_mask');

fprintf(fiD,['\n \nAll ROIs with ' num2str(min_vox) ' voxels or less were ignored. \n \nROIs=\n']);

for iter=1:length(rois)-1
   fprintf(fiD,['' num2str(rois(iter)) ',']);
end
if length(rois)>0
    fprintf(fiD,['' num2str(rois(end)) '']);
end

fclose(fiD);

%No. of non-zero voxels in .nii-file:
disp(['' num2str(length(find(data2~=0))) ' non-zero voxels.']);