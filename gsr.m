function [GS,GSregressed_data] = gsr(varargin)
% performs GSR on data. Every voxel in the brainmask-file whose intensity is
% larger than 0 is considered part of the mask.
% --------------------------
% Input (in this order):
%   file:       The nifti-file from which you want to regress the global
%               signal. May be a single file or a cell array for a number of files.
%   brainmask:  Defines the space over which the GS is computed. All voxels
%               with a value of >0 will be interpreted as one large region over which
%               to average the voxe-values to obtain the global signal.
%               May be a single file or a cell array for a numbver of
%               files.
%   strg_dir:   Storage directory, where the regressed nifti-files shall be
%               stored.
%   quitTF (optional):  Boolean function to determine whether to close matlab 
%                       after function terminates.
%                       DEFAULT=1
%
% Output:
%   GS :                GS is the global signal.
%   GSregressed_data :  The data with the GS regressed out.
%
% --------------------------

files=varargin{1};
brainmasks=varargin{2};

if ~iscell(files)
    filesTEMP=cell(1,1);
    filesTEMP{1}=files;
    files=filesTEMP;
    clear filesTEMP
    
    bmTEMP=cell(1,1);
    bmTEMP{1}=brainmasks;
    brinmasks=bmTEMP;
    clear bmTEMP
end

if nargin>2
    strg_dir=varargin{3};
end

if nargin==4
   quitTF=varargin{4};
else
   quitTF=1; 
end

disp('PERFORMING GSR...');


for fl=1:length(files)
    file=files{fl};
    brainmask=brainmasks{fl};

    disp(['File: ' file '']);
    disp(['Brainmask: ' brainmask '']);

    % This script is to perform a GSR on a nii-file and store the data with GSR as a
    % separate file
    newfile='';

    splitfile = strsplit(file,filesep);
    if nargin>2
        newfile=strg_dir;
    else
        if isunix
            newfile = strcat(newfile,filesep);
        end

        for i = 1:(length(splitfile)-1)
            if ~isempty(splitfile{i})
                newfile = strcat(newfile,splitfile(i)); 
                newfile = strcat(newfile,filesep);
            end
        end
    end

    % Contains the global signal (is mat-file)
    GSfile=strcat(newfile,'GS_');
    GSfile=strcat(GSfile,splitfile{end});

    % nifti file with global signal regressed out
    GSregr=strcat(newfile,'GSR');
    GSregr=strcat(GSregr,splitfile{end});


    filetype=file(end-2:end);
    filename=file(1:end-4);

    if filetype~='nii'
       error('Files must be .nii-files'); 
    end

    filetype_brain=brainmask(end-2:end);

    if filetype_brain~='nii'
       error('Brainmask must be .nii-file'); 
    end


    [data,pixdim,orient,dtype,slice_code] = readnifti(['' file '']);

    [GS,GSregressed_data]=performgsr(data,brainmask);

    fileinfo=readniftiheader(['' filename ''],'nii');
    %%%READS HEADER INFO. DO NOT TOUCH--------------------------------------%%%
    fid=fopen(file,'r');
    sizeof_hdr=fread(fid,1,'int32');
    if sizeof_hdr==540
       error('This script is only for Nifti-1 files. Your file has file-format nifti-2') 
    end
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
    qoffset_y=fread(fid,1,'float');
    qoffset_z=fread(fid,1,'float');
    srow_x=fread(fid,4,'float');
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

    %%%%WRITE to GSRregressed_* FILE WHICH CONTAINS REGRESSED DATA
    %First, recreate the file
    if iscell(GSregr)
        fileID = fopen(['' GSregr{1} ''],'w');
    else
        fileID = fopen(['' GSregr ''],'w');
    end
    %----------TRANSFERS HEADER INFORMATION TO NEWFILE----------------------%%%
    %Create new file and write to new file
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
    fwrite(fileID,magic,'char'); %reached bit is bit 348
    %%%%--------------------------------------------------------------------%%%
    %VoxOffset is greater than 348, so imaging data does not immediately follow the header
    fwrite(fileID,offset_contents,'int8');

    %%%%TRANSFER REGRESSGS ED DATA MATRIX TO NEW FILE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fwrite(fileID,GSregressed_data,dtype);
    %%%%--------------------------------------------------------------------%%%
    fclose(fileID);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Save file containing global signal
    if iscell(GSfile)
        GSfile=GSfile{1};
    end
    GSfile=strcat(GSfile(1:end-3),'mat');
    disp(['Global Signal is stored in .mat- file ' GSfile '.']);
    fprintf('\n');
    save(GSfile,'GS');

end


disp('Done!');

if quitTF
	quit
end
end
