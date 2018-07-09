function [fileinfo] = readniftiheader(filename,filenametype)
%This fucntion is used by gsr.m and readnifti.m
%output "fileinfo" fields:
%   dimension
%   pixdim
%   datatype
%   bitpix
%   voxoffset
%   byteswap
%   formcodes
%   qmatrix
%   smatrix
%   slicecode
%
% Open Header
%    filename is the name of the file
%    filenametype is either 'img' or 'nii'
%

if isempty(filenametype)
  error('not implemented yet');
else
  if strcmp(filenametype,'img')
    fid = fopen([filename,'.hdr'],'r');
    if fid < 0
      error('Error opening header file (check file directory!)');
    end
  else
    fid = fopen([filename,'.nii'],'r');
    if fid < 0 %if fopen cannot open the file, then fileID is -1
      error('Error opening header file (check file directory!)');
    end
  end
end
%fseek moves to a specified position in the file. The dimension information
%has offset 39Bytes in the nifti header and is one Byte long. 'bof'
%specifies that the offset (here 40) is counted from the Beginning Of the
%File
fseek(fid,40,'bof');
%read the output into an array of length 8 (is stored in header as array of length 8, each entry is a int16)
dimension = fread(fid,8,'int16');
%fread(fid,8,'int16') says that the read content is to be stored in an
%array of length 8 with each entry being an int16
byteswap = 'n';

% The next part tries to figure out if we have a byte swapping problem by
% looking at the dimension bit in the header file. If it's bigger than 10
% (we probably aren't dealing with more than 10-dimensional data) then it
% tries byte swapping. If it still fails then it outputs an error.
% dimension
if (dimension(1) > 10)
  byteswap = 'b';
  fclose(fid);
  if isempty(filenametype)
    error('not implemented yet');
  else
    if strcmp(filenametype,'img')
      fid = fopen([filename,'.hdr'],'r',byteswap);
      if fid < 0
        error('Error opening header file');
      end
    else
      fid = fopen([filename,'.nii'],'r',byteswap);
      if fid < 0
        error('Error opening header file');
      end
    end
  end
  fseek(fid,40,'bof');
  dimension = fread(fid,8,'int16');
  if (dimension(1) > 10)
    byteswap = 'l';
    fclose(fid);
    if isempty(filenametype)
      error('not implemented yet');
    else
      if strcmp(filenametype,'img')
        fid = fopen([filename,'.hdr'],'r',byteswap);
        if fid < 0
          error('Error opening header file');
        end
      else
        fid = fopen([filename,'.nii'],'r',byteswap);
        if fid < 0
          error('Error opening header file');
        end
      end
    end
    fseek(fid,40,'bof');
    dimension = fread(fid,8,'int16');
    if (dimension(1) > 10)
      fclose(fid);
      error('Error opening file. Dimension argument is not valid');
    end
  end
end

% Read information from header file
% Note: there is information in the header that is not currently read
% like the intent codes.

% datatype, possible values range from 0 to 2304 which gives a max of
% 16bit=2byte
fseek(fid,40+30,'bof');
datatype = fread(fid,1,'int16');

% bits per pixel
fseek(fid,40+32,'bof');
bitpix = fread(fid,1,'int16');

% pixel dimension
fseek(fid,40+36,'bof');
pixdim = fread(fid,8,'float');

% data offset. The field vox_offset indicates, for single files, the byte
% offset before the imaging data starts
fseek(fid,108,'bof');
voxoffset = fread(fid,1,'float');

%Slice acquisition order
fseek(fid,122,'bof');
slicecode = fread(fid,1,'char');

% orientation 
fseek(fid,252,'bof');
%read qform_code and sform_code
formcodes = fread(fid,2,'int16');
fseek(fid,256,'bof');
qmatrix = fread(fid,6,'float'); %float is 32 bits (=4 byte), float('double') is 64 bits
fseek(fid,280,'bof');
smatrix = fread(fid,12,'float');

% Close header file
fclose(fid);

fileinfo = struct('dimension',dimension,...
                  'pixdim',pixdim,...
                  'datatype',datatype,...
                  'bitpix',bitpix,...
                  'voxoffset',voxoffset,...
                  'byteswap',byteswap,...
                  'formcodes',formcodes,...
                  'qmatrix',qmatrix,...
                  'smatrix',smatrix,...
                  'slicecode',slicecode);
end