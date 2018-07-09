%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               READ NIFTI FILE HEADERS                %
%             Only for Nifti-1 file format             %
%                                                      %
% Author: Sascha Frölich B.Sc.          Dec. 2016      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [header] = readniftifileheader(file)

%Script is used by createplots.m and lagmapsanal.m

[data,pixdim,orient,dtype,slice_code] = readnifti(['' file '']);

% NIFTI: HEADER SIZE IS 348B
% NIFTI2: HEADER SIZE IS 540B
%%%READS HEADER INFO. DO NOT TOUCH--------------------------------------%%%
fid=fopen(['' file ''],'r');
header.sizeof_hdr=fread(fid,1,'int32');
if header.sizeof_hdr==540
   error('This script is only for Nifti-1 files. Your file has file-format nifti-2') 
end
header.notused=fread(fid,35,'int8');
header.dim_info=fread(fid,1,'char');
header.dim=fread(fid,8,'short'); %stores the dimensions (1st entry) and no of voxels (ie matrix entries)
header.intent_p1=fread(fid,1,'float');
header.intent_p2=fread(fid,1,'float');
header.intent_p3=fread(fid,1,'float');
header.intent_code=fread(fid,1,'short');
header.datatype=fread(fid,1,'short');
header.bitpix=fread(fid,1,'short');
header.slice_start=fread(fid,1,'short');
header.pixdim=fread(fid,8,'float'); %dimensions of the voxels (ie matrix entries). The first value in pixdim has a special meaning, it should always be either -1 or 1
header.vox_offset=fread(fid,1,'float'); %Offset in bytes of the data inside the file
header.scl_slope=fread(fid,1,'float');
header.scl_inter=fread(fid,1,'float');
header.slice_end=fread(fid,1,'short');
header.slice_code=fread(fid,1,'char');
header.xyzt_units=fread(fid,1,'char');
header.cal_max=fread(fid,1,'float'); %intended max display intensity when image is opened
header.cal_min=fread(fid,1,'float'); %intended min display intensity when image is opened
header.slice_duration=fread(fid,1,'float');
header.toffset=fread(fid,1,'float');
header.notused2=fread(fid,8,'int8');
header.descrip=fread(fid,80,'char');
header.aux_file=fread(fid,24,'char');
header.qform_code=fread(fid,1,'short');
header.sform_code=fread(fid,1,'short');
header.quatern_b=fread(fid,1,'float');
header.quatern_c=fread(fid,1,'float');
header.quatern_d=fread(fid,1,'float');
header.qoffset_x=fread(fid,1,'float');
header.qoffset_y=fread(fid,1,'float');
header.qoffset_z=fread(fid,1,'float');
header.srow_x=fread(fid,4,'float');
header.srow_y=fread(fid,4,'float');
header.srow_z=fread(fid,4,'float');
header.intent_name=fread(fid,16,'char');
header.magic=fread(fid,4,'char');
%if fseek is pointed to 76, then fread() will start reading out at 76Bytes
%fread(fid,8,'float=>int8'); will read 8 floats from the file, then store
%them in an array of length 8, where every entry is an int8. Default:
%stores to double
fclose(fid);
%%%-----------------------------------------------------------------------%

end
