function []=TD_corr(varargin)
% Function TD_corr computes the convergence of random group level TD
% matrices of subsets of the whole group and corresponds those matrices to
% the actual matrix TD.
%
% INPUT (in this order):
%       groupfile:      matlab-file containing a group lebvel time delay
%                       matrix TD and fields "filepath" and "filenames".
%                       Those fields indicate which first level matlab
%                       files were used for computation of group level TD.
%       no_iterations:  The number of iterations per subset size.
%       savenewfile:    Boolean function whether or not to save the TD
%                       convergence statistics in a new matlab file.
% OUTPUT:
%       none

disp('FUNCTION TD_corr');

groupfile=varargin{1};
no_iterations=varargin{2};
savenewfile=varargin{3}; %tf toggle

m=load(['' groupfile '']);

disp(['' num2str(length(m.filenames)) ' files.'])

max_subset=length(m.filenames);

% crr is a vector of correlation values (between subset groupTD and reference TD) per subset size
% p_vals is a vector of p-values for given correlations
crr=zeros(1,max_subset);
p_vals=zeros(1,max_subset);

fast=1;
   
% Very time-consuming version below
if fast
    GB=m.noofROIS_mask1*m.noofROIS_mask2*64*length(m.filenames)/(8*1024*1024*1024);
    prompt=strcat('This will need ',num2str(GB),' GB of memory. What do you want to do?');
    choice=questdlg(prompt,'Necessary Memory Capacity','Continue','Convert to single (halfs the requirement)','Cancel','Cancel');
    
    switch choice
        case 'Continue'
            AllTDs=zeros(m.noofROIS_mask1,m.noofROIS_mask2,length(m.filenames));
            grpTD=zeros(m.noofROIS_mask1,m.noofROIS_mask2);
            for fl=1:length(m.filenames)
                file=strcat(m.pathtofiles,m.filenames{fl});
                tdfile=load(file);
                tdfile.TD
                AllTDs(:,:,fl)=tdfile.TD;
            end  
        case 'Convert to single (halfs the requirement)'
            AllTDs=single(zeros(m.noofROIS_mask1,m.noofROIS_mask2,length(m.filenames)));
            grpTD=single(m.noofROIS_mask1,m.noofROIS_mask2);
            for fl=1:length(m.filenames)
                file=strcat(m.pathtofiles,m.filenames{fl});
                tdfile=load(file);
                AllTDs(:,:,fl)=single(tdfile.TD);
            end  
        case 'Cancel'
            error('Process terminated by user.','error','error');
    end
   
    if ~isempty(choice)
            commandwindow
           for subset_size=1:max_subset

               subset_size

               for iter=1:no_iterations

                   iter_vec=randperm(length(m.filenames),subset_size);

                   grpTD=nanmean(AllTDs(:,:,iter_vec),3);
                   
                   [crrltn,pval]=corrcoef(grpTD,m.TD,'rows','complete');

                   crr(subset_size)=crr(subset_size)+crrltn(1,2)/no_iterations;
                   p_vals(subset_size)=p_vals(subset_size)+pval(1,2)/no_iterations;

               end   

           end

           if savenewfile
                TDConvMaxSubset=max_subset;
                TDConvCorr=crr;
                TDConvP=p_vals;
                TDConvNoIt=no_iterations;
                save(groupfile,'TDConvMaxSubset','-append');
                save(groupfile,'TDConvCorr','-append');
                save(groupfile,'TDConvP','-append');
                save(groupfile,'TDConvNoIt','-append');
           end

           figure;
           plot(1:max_subset,crr);
           xlabel('subset size');
           ylabel('Pearson r');

           disp('Done!');
    end
    
else
    commandwindow
   for subset_size=1:max_subset
       
       subset_size
       
       for iter=1:no_iterations
       
           iter_vec=randperm(length(m.filenames),subset_size);
           
           grpTD=zeros(m.noofROIS_mask1,m.noofROIS_mask2);
           
           for fl=1:subset_size
            
                file=strcat(m.pathtofiles,m.filenames{iter_vec(fl)});
                tdfile=load(['' file '']);
                grpTD=grpTD+tdfile.TD/subset_size;
            
           end
          
           [crrltn,pval]=corrcoef(grpTD,m.TD,'rows','complete');
          
           crr(subset_size)=crr(subset_size)+crrltn(1,2)/no_iterations;
           p_vals(subset_size)=p_vals(subset_size)+pval(1,2)/no_iterations;
           
       end   
       
   end
    
   if savenewfile

        save(newfile,'groupfile');
        save(newfile,'max_subset','-append');
        save(newfile,'crr','-append');
        save(newfile,'p_vals','-append');
        save(newfile,'no_iterations','-append');
       
        
        save(groupfile,'statsfilepath','-append');
        save(groupfile,'statsfilename','-append');

   end

   figure;
   plot(1:max_subset,crr);
   xlabel('subset size');
   ylabel('Pearson r');

    disp('Done!');
end

if 0
%%%DO THE SAME WITH RANDOM MATRICES
if m.noofROIS_mask1 < 1000 || m.noofROIS_mask2 < 1000
    
    for method=1:3
        clear TD;
        %method 1 create a no. of completely random TDs.
        if method==1
                for i=1:max_subset
                    TD(:,:,i)=rand(m.noofROIS_mask1,m.noofROIS_mask2);
                end
        %method 2: force TD to be skew-symmetric and have zeroes on
        %diagonal (random values only positive)
        elseif method==2
                TD=zeros(m.noofROIS_mask1,m.noofROIS_mask2);
                for i=1:max_subset
                    for row=1:m.noofROIS_mask1
                        for col=row+1:m.noofROIS_mask2
                            TD(row,col,i)=rand;
                            TD(col,row,i)=-TD(row,col,i);
                        end
                    end
                end
        %as method 2, but with random values being negative or positive
        elseif method==3
                TD=zeros(m.noofROIS_mask1,m.noofROIS_mask2);
                for i=1:max_subset
                    for row=1:m.noofROIS_mask1
                        for col=row+1:m.noofROIS_mask2
                            sign=rand-0.5;
                            if sign<0
                                sign=-1;
                            else
                                sign=1;
                            end

                            TD(row,col,i)=sign*rand;
                            TD(col,row,i)=-TD(row,col,i);
                        end
                    end
                end
        end





        grp_TD=mean(TD,3);

       crr_rand=zeros(1,max_subset);
       p_vals_rand=zeros(1,max_subset);

       for subset_size=1:max_subset

           subset_size

           for iter=1:no_iterations

               iter_vec=randperm(max_subset,subset_size);

               for fl=1:subset_size

                    TDtemp(:,:,fl)=TD(:,:,iter_vec(fl));

               end

               grpTDtemp=nanmean(TDtemp,3);

               [crrltn,pval]=corrcoef(grpTDtemp,grp_TD);

               crr_rand(subset_size)=crr_rand(subset_size)+crrltn(1,2)/no_iterations;
               p_vals_rand(subset_size)=p_vals_rand(subset_size)+pval(1,2)/no_iterations;

           end   

       end

       figure;
       plot(1:max_subset,crr_rand,'red',1:max_subset,crr,'green')
       title(['green=real,red=random (method ' num2str(method) ')']);

       if savenewfile

           if method==1
                crr_rand_method1=crr_rand;
                p_vals_rand_method1=p_vals_rand;
                save(newfile,['crr_rand_method1'],'-append');
                save(newfile,['p_vals_rand_method1'],'-append');
           elseif method==2
                crr_rand_method2=crr_rand;
                p_vals_rand_method2=p_vals_rand;
                save(newfile,['crr_rand_method2'],'-append');
                save(newfile,['p_vals_rand_method2'],'-append');
           elseif method==3
                crr_rand_method3=crr_rand;
                p_vals_rand_method3=p_vals_rand;
                save(newfile,['crr_rand_method3'],'-append');
                save(newfile,['p_vals_rand_method3'],'-append');
           end

       end


    end
   
    
end

end

    
end