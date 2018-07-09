function []=groupnifti(filename,filepath)
% groupnifti(filename,filepath)
%   INPUT:
%   filename:   cell array containing the filenames of the files from which
%               to create the group-level file.
%   filepath:   Directory containing all files in "filename". With trailing
%               file seperator.
%
% Script creates a group-level nifti file of all files listed in "filename"

disp(['Computing Group-Level Nifti-File of ' num2str(length(filename)) ' files.']);

[data,pixdim,~,~,~]=readnifti(['' strcat(filepath,filename{1}) '']);
[rows,columns,pages]=size(data);
groupdata=zeros(rows,columns,pages);

file1=strcat(filepath,filename{1});
[filedata,~,~,~,~]=readnifti(['' file1 '']);
DivMat=filedata;
DivMat(:)=0;

for fl=1:length(filename)
    file=strcat(filepath,filename{fl});
    [filedata,~,~,~,~]=readnifti(['' file '']);
    
    
    % DivMat is "divider matrix". For every entry in the niftidata-matrix
    % it contains the no. through whiich to divide for a
    % nanmean-group-level computation.
    
    % Convert every entry in filedata that is NaN to 0 and every other
    % entry to 1
    filedata(find(~isnan(filedata)))=1;
    filedata(find(isnan(filedata)))=0;
    
    DivMat=DivMat+filedata;
    
end

    % If DivMat contains zeros, that means that every file has NaN at that
    % position. In that case, filedata will have entry 0 in this position
    % due to manipulation of filedata in following loop. Therefore, set
    % DivMat to 1, as 0/1=0.
    disp(['' num2str(length(find(DivMat==0))) ' NaN-voxels set to zero.']); 
    DivMat(find(DivMat==0))=1;
   

for fl=1:length(filename)
    file=strcat(filepath,filename{fl});
    [filedata,~,~,~,~]=readnifti(['' file '']);
    
    
    % DivMat is "divider matrix". For every entry in the niftidata-matrix
    % it contains the no. through whiich to divide for a
    % nanmean-group-level computation.
    
    filedata(find(isnan(filedata)))=0;
    
    groupdata=groupdata+filedata./DivMat;
end

newfile=strcat(filepath,filename{1});
newfile=strsplit(newfile,filesep);
newfile=newfile{end};
newfile=strcat(filepath,'GROUPFILE_',newfile(1:end-4),'.nii')

strcat(filepath,filename{1}),newfile
write_nifti(strcat(filepath,filename{1}),newfile,groupdata);

commandwindow
disp('Computation of Group-Level Nifti-File terminated.');

end