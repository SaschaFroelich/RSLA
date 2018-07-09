cd /data_august/sascha/cer/trunk/lagged_covar/lag_anal_GUI_2masks/masks

[status,out]=system('ls | grep .mat');
mitramasks=strsplit(out,'\n');

if isempty(mitramasks{end})
   mitramasks=mitramasks{1:end-1}; 
end

no_masks=length(mitramasks);

%matrix ROIS: each column is a mask, with the intensities of voxels in that
%mask
ROIS=[];

for iter=1:no_masks
    m=load(['' mitramasks{iter} '']);
    
    if iter>1
       [rows,columns]
    end
    
end