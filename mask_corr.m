function []=mask_corr(groupfile,mask)
%Takes as input a groupfile.mat "file" (complete with path!) and a ROI-mask (complete with path!). File can be nifti or .mat containing a field "m.BOLDfile" It then
%computes the average zerolag-correlation across the ROIs of given mask.

if strcmp(groupfile(end-2:end),'nii')
    
elseif strcmp(groupfile(end-2:end),'mat')
    m=load(['' groupfile '']);
    
    if isfield(m,'filenames')
       
        filenames=m.filenames;
        mask_subject_average=zeros(1,length(filenames));
        
        for file_no=1:length(filenames)
           
            disp(['Processing file no' num2str(file_no) '.']);
            single_subject_file=strcat(m.pathtofiles,filenames{file_no});
            
            m2=load(['' single_subject_file '']);
            
            [ROISCOURSE, INTENS, numberofrois] = coursematrix(m2.BOLDfile,mask);
            
            
            zerolagmask=mask;
            ZL_zerolagmask=zeros(numberofrois,numberofrois);
            
            for roi1=1:numberofrois
                for roi2=1:numberofrois
                    ZL_zerolagmask(roi1,roi2)=corr2(ROISCOURSE(roi1,:),ROISCOURSE(roi2,:));
                end
            end
            
            ZL_zerolagmask_average=mean(mean(ZL_zerolagmask));

            mask_subject_average(file_no)=ZL_zerolagmask_average;
            
        end
        
    end
    
    zerolagmask_average=mean(mask_subject_average);
    
    zerolagmask=strsplit(zerolagmask,filesep);
    zerolagmask=zerolagmask{end};
    
    newfile=strcat(groupfile(1:end-4),'_AvgCorr_',zerolagmask(1:end-4),'.mat');
    
    
    
    save(newfile,'zerolagmask');
    save(newfile,'mask_subject_average','-append');
    save(newfile,'zerolagmask_average','-append');
    
    
end




end
