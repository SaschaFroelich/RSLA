function []=flipfile()
% This script "flips" the header of a nifti-file from left to right. 
% Example:
% If the rightmost voxel of a nifti-file is at x-value 0 and mm-value 90,
% after flipping that same, rightmost voxel will have the largest possible
% x-value with decreasing x-values going from right to left, but the same
% mm-value. Thus the matrix-entries in the nifti-files are flipped, while
% the actual positions in space of the anatomical structures are not
% changetocd. This is sometimes needed for reasons of compatibility.

[flipfilename, flipfilepath]=uigetfile(['' pwd '/*.nii'],'MultiSelect','On','Select nifti-file(s) which to flip');

if ~iscell(flipfilename)
    datafiles=cell(1);
    datafiles{1}=flipfilename;
else
    datafiles=flipfilename;
end

clear datafilename

for file=1:length(datafiles)
flipfilename=datafiles{file};
    
[data,~,~,dtype,~] = readnifti(['' strcat(flipfilepath,flipfilename) '']);

newfile=strcat(flipfilepath,'f',flipfilename(1:end-4),'.nii');
[data2,~,~,~,~] = readnifti(['' strcat(flipfilepath,flipfilename) '']);

%THE NIFTI HEADER IS OF SIZE 348BYTES
%%%READS HEADER INFO. DO NOT TOUCH--------------------------------------%%%
fid=fopen(['' strcat(flipfilepath,flipfilename) ''],'r');
sizeof_hdr=fread(fid,1,'int32');
notused=fread(fid,35,'int8');
dim_info=fread(fid,1,'char');
dim=fread(fid,8,'short'); %stores the dimensions (1st entry) and no of voxels (ie matrix entries)
intent_p1=fread(fid,1,'float');
intent_p2=fread(fid,1,'float');
intent_p3=fread(fid,1,'float');
intent_code=fread(fid,1,'short');
datatype=fread(fid,1,'short');
bitpix=fread(fid,1,'short');
slice_start=fread(fid,1,'short');
pixdim=fread(fid,8,'float'); %dimensions of the voxels (ie matrix entries). The first value in pixdim has a special meaning, it should always be either -1 or 1
% Change header
pixdim
pixdim(1)=pixdim(1)*(-1);
vox_offset=fread(fid,1,'float'); %Offset in bytes of the data inside the file
scl_slope=fread(fid,1,'float');
scl_inter=fread(fid,1,'float');
slice_end=fread(fid,1,'short');
slice_code=fread(fid,1,'char');
xyzt_units=fread(fid,1,'char');
cal_max=fread(fid,1,'float'); %intended max display intensity when image is opened
cal_min=fread(fid,1,'float'); %intended min display intensity when image is opened
slice_duration=fread(fid,1,'float');
toffset=fread(fid,1,'float');
notused2=fread(fid,8,'int8');
descrip=fread(fid,80,'char');
aux_file=fread(fid,24,'char');
qform_code=fread(fid,1,'short');
sform_code=fread(fid,1,'short');
quatern_b=fread(fid,1,'float');
quatern_c=fread(fid,1,'float');
quatern_d=fread(fid,1,'float');
qoffset_x=fread(fid,1,'float');
% Change header
qoffset_x=qoffset_x*(-1);
qoffset_y=fread(fid,1,'float');
qoffset_z=fread(fid,1,'float');
srow_x=fread(fid,4,'float');
% Change header
srow_x
srow_x=srow_x*(-1);
srow_y=fread(fid,4,'float');
srow_z=fread(fid,4,'float');
intent_name=fread(fid,16,'char');
magic=fread(fid,4,'char');
%if fseek is pointed to 76, then fread() will start reading out at 76Bytes
%fread(fid,8,'float=>int8'); will read 8 floats from the file, then store
%them in an array of length 8, where every entry is an int8. Default:
%stores to double
offset_contents=fread(fid,vox_offset-348,'int8');
fclose(fid);
%%%-----------------------------------------------------------------------%

%%%REDEFINE HEADER INFO HERE----------------------------------------------%
%%%PARTS OF HEADER THAT ARE NOT CALLED HERE REMAIN UNCHANGED--------------%

%%%-----------------------------------------------------------------------%

%REDEFINE data to write here----------------------------------------------%

for i=1:size(data,1)
   for j=1:size(data,2)
       for k=1:size(data,3)
            data_new(i,j,k)=data2(end-i+1,j,k);
       end
   end
end
%%%-----------------------------------------------------------------------%

%%%WRITE TO NEW FILE-----------------------------------------------------%
fileID = fopen(['' newfile ''],'w');
fwrite(fileID,sizeof_hdr,'int32');
fwrite(fileID,notused,'int8');
fwrite(fileID,dim_info,'char');
fwrite(fileID,dim,'short');
fwrite(fileID,intent_p1,'float');
fwrite(fileID,intent_p2,'float');
fwrite(fileID,intent_p3,'float');
fwrite(fileID,intent_code,'short');
fwrite(fileID,datatype,'short');
fwrite(fileID,bitpix,'short');
fwrite(fileID,slice_start,'short');
fwrite(fileID,pixdim,'float');
fwrite(fileID,vox_offset,'float');
fwrite(fileID,scl_slope,'float');
fwrite(fileID,scl_inter,'float');
fwrite(fileID,slice_end,'short');
fwrite(fileID,slice_code,'char');
fwrite(fileID,xyzt_units,'char');
fwrite(fileID,cal_max,'float');
fwrite(fileID,cal_min,'float');
fwrite(fileID,slice_duration,'float');
fwrite(fileID,toffset,'float');
fwrite(fileID,notused2,'int8');
fwrite(fileID,descrip,'char');
fwrite(fileID,aux_file,'char');
fwrite(fileID,qform_code,'short');
fwrite(fileID,sform_code,'short');
fwrite(fileID,quatern_b,'float');
fwrite(fileID,quatern_c,'float');
fwrite(fileID,quatern_d,'float');
fwrite(fileID,qoffset_x,'float');
fwrite(fileID,qoffset_y,'float');
fwrite(fileID,qoffset_z,'float');
fwrite(fileID,srow_x,'float');
fwrite(fileID,srow_y,'float');
fwrite(fileID,srow_z,'float');
fwrite(fileID,intent_name,'char');
fwrite(fileID,magic,'char');

%AS HEADER IS OF SIZE 348BYTES AND DATA IS SUPPOSED TO START AS VOX_OFFSET,
%INSERT ADDITIONAL BYTES.
fwrite(fileID,offset_contents,'int8');

fwrite(fileID,data_new,dtype);
fclose(fileID);
disp('Writing to file.');
end

%%%-----------------------------------------------------------------------%
disp('Done.');


