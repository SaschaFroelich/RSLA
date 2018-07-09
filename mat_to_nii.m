function []=mat_to_nii(varargin)

file=varargin{1}; % matfile complete with path
plotting=varargin{2}; % 0=lagthread 1=lagmaps 2=lag projection

if plotting~=2
   seed=varargin{3}; 
   
   if plotting==1
       %reference map contains the grouplevel file for the given files. If
       %a single subject contains a NaN-entry at voxel i,j,k, its value in
       %the niftifile will be replaced by its value in the grouplevelfile
      reference_map=varargin{4};
       
      refmap=load(['' reference_map '']);
      
   end
   
end

m=load(file);

if ~strcmp(m.mask1,m.mask2)
    
    error('File uses two different masks','Error','error');
    
end

if plotting==0
    
    disp('Made sure that lag threads are correct? Will use entries of L as they are.');
    
    lag_vec=m.L(:,seed);
    
    newfile=strcat(file(1:end-4),['_lag_thread' num2str(seed) '.nii'])
    
elseif plotting==1
    
    I=find(seed==m.TDrois_cols);
    I2=find(seed==refmap.TDrois_cols);
    
    lag_vec=m.TD(:,I);
    refmap_vec=refmap.TD(:,I2);
    
    for i=1:length(lag_vec);
       if isnan(lag_vec(i))
           lag_vec(i)=refmap_vec(i);
       end
    end
    
    newfile=strcat(file(1:end-4),['_LagMap' num2str(seed) '.nii'])
    
    
elseif plotting==2
    if isfield(m,'GroupLagProjection')
        lag_vec=m.GroupLagProjection;
    else
        lag_vec=nanmean(m.TD);
    end
    
    newfile=strcat(file(1:end-4),'_LagProjection.nii')
    
end

[data,pixdim,orient,dtype,slice_code] = readnifti(m.mask1);   
[rows,columns,pages]=size(data);

data_new=data;

for row=1:rows
    for col=1:columns
        for page=1:pages
            %Check which ROI we are handling
            roi_intens=data_new(row,col,page);
            
            if roi_intens~=0
            
                %Check which position in the lag_vec this ROI occupies
                roi_pos=find(roi_intens==m.rois_in_mask1);

                %Replace the value in data_new by the lag value in lag_vec
                data_new(row,col,page)=lag_vec(roi_pos);
                
            end
            
        end
    end
end
    
    
%THE NIFTI HEADER IS OF SIZE 348BYTES
%%%READS HEADER INFO. DO NOT TOUCH--------------------------------------%%%
fid=fopen(['' m.mask1 ''],'r');
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
fclose(fid);
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
additional=vox_offset-348;
empty=0;
empty=uint8(empty);

%Each iteration writes one byte
for iter=1:additional
   fwrite (fileID,empty,'uint8');
end

fwrite(fileID,data_new,dtype);
fclose(fileID);


%%%-----------------------------------------------------------------------%
    
    
disp('Done!');

quit

end