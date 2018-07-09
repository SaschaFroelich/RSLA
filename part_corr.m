clc;
clear all;
close all;

% This script can do two things, depending on the selection ("choice", see below):

% Option 1
% This file stores the desired lags as columns to the psychometric
% observations in psychometr_msrs

% Option 2
% this file takes the selected file with psychometric msrs and performs
% partial correlation and multiple regression with the chosen regressors.

choice=questdlg('Choose option','What do you want to do?','Add lag to psychometrics file','Perform partial correlation','partcorr');

switch choice
    case 'lag'

        % sourcefile is the file from which to take all psychometric measures (stored as a cell)
        [sourceFILENAME, sourcePATHNAME] = uigetfile(['' pwd filesep '*.mat'],'Choose sourcefile containing cell table "psychometr_msrs"');
        sourcefile=strcat(sourcePATHNAME,sourceFILENAME)
        m=load(sourcefile);

        % newfile is wehere the new table will be stored. The new data corresponds
        % to the table in sourcefile plus some appended columns
        newfile=(strcat(sourcePATHNAME,sourceFILENAME(1:end-4),'_plus'));

        % Now this is important: Choose a nifti-file, or a mat-file that stands represantative of
        % all the nifti-files from which you want to extract a specific lag.
        [lagFILENAME, lagPATHNAME] = uigetfile(['' pwd filesep '*.mat;*.nii'],'Choose lagfile (nii or mat) single-subject file, represantative of all included files. (Can be from either grp1 or grp2)');
        niifile=strcat(lagPATHNAME,lagFILENAME)

        % In its name, niifile should have the subject's ID as it is stored exactly
        % in the first column of the psychometrics file (the "name" column).
        % Possibly, niifile has a longer name than that. Determine suffix & prefix:
        for ii=1:size(m.psychometr_msrs,1);
            k=strfind(lagFILENAME,m.psychometr_msrs{ii,1});
            if ~isempty(k)
                prefix=lagFILENAME(1:k-1);
                suffix=lagFILENAME(k+length(m.psychometr_msrs{ii,1}):end);
            end
        end

        if strcmp(lagFILENAME(end-3:end),'.nii')
            prompt={ 'i (MATLAB)', 'j (MATLAB)', 'k (MATLAB)' };
            dlg_title='Enter voxel coordinates in MATLAB space';
            num_lines=1;
            answer=inputdlg(prompt,dlg_title,num_lines);
            i=str2double(answer{1});
            j=str2double(answer{2});
            k=str2double(answer{3});
        elseif strcmp(lagFILENAME(end-3:end),'.mat')
            prompt='Enter ROIs (intensity values as given in column list) over which to average lag value. Sperate by comma.';
            dlg_title='Input';
            num_lines=1;
            rois_temp=inputdlg(prompt,dlg_title,num_lines);
            rois_temp=strsplit(rois_temp{1},',');

            rois=[];
            for entry=1:length(rois_temp)
                if ~isempty(rois_temp(entry))
                    rois=[rois, str2double(rois_temp(entry))];
                end
            end
        end


        prompt={ 'Enter name of new Column of the psychometr_msrs table' };
        dlg_title='Enter name of new Column:';
        num_lines=1;
        column_name=inputdlg(prompt,dlg_title,num_lines);


        if strcmp(lagFILENAME(end-3:end),'.nii')
             psychometr_msrs=m.psychometr_msrs;
             cluster_lags=cell(length(psychometr_msrs),1);
             cluster_lags{1,1}=column_name{1};

             psychometr_msrs=[psychometr_msrs, cluster_lags];

            for line=2:size(psychometr_msrs,1)
                line

                lagPATHNAMEparent=strsplit(lagPATHNAME,'grp');
               if psychometr_msrs{line,2}==1
                   file=strcat(lagPATHNAMEparent{1},'grp1/',prefix,psychometr_msrs{line,1},suffix);
               elseif psychometr_msrs{line,2}==2
                   temp=[];
                   file=strcat(lagPATHNAMEparent{1},'grp2/',prefix,psychometr_msrs{line,1},suffix);
               end

                if exist(['' file ''])
                    [data,~,~,~,~]=readnifti(['' file '']);
                    psychometr_msrs{line,end}=data(i,j,k);
                end

            end

        elseif strcmp(lagFILENAME(end-3:end),'.mat')

             psychometr_msrs=m.psychometr_msrs;
             cluster_lags=cell(length(psychometr_msrs),1);
             cluster_lags{1,1}=column_name{1};
             psychometr_msrs=[psychometr_msrs, cluster_lags];

            for line=2:length(psychometr_msrs)
                line=line
                file=strcat(lagPATHNAME,psychometr_msrs{line,1},'.mat');

                rois_idx=[];
                if exist(['' file ''])
                    subj=load(['' file '']);

                    for jj=rois
                        rois_idx=[rois_idx find(jj==subj.TDrois_cols)];
                    end

                    avrg_lag=mean(subj.ColumnMeans(rois_idx));
                    psychometr_msrs{line,end}=avrg_lag;

                end


            end
        end

         save(newfile,'psychometr_msrs');
         disp('Done.');

     case 'Perform partial correlation'
         
        % sourcefile is the file from which to take all psychometric measures (stored as a cell)
        [sourceFILENAME, sourcePATHNAME] = uigetfile(['' pwd filesep '*.mat'],'Choose sourcefile containing cell table "psychometr_msrs"');
        sourcefile=strcat(sourcePATHNAME,sourceFILENAME)
        m=load(sourcefile);

         % The following vector stores the columns of the psychometrics
         % file that serve as partial regressors.
         regressors= [4 10 11]; % GA, BW, Opti
                
         
         % The following vector contains the columns which serve as
         % dependent variables in the partial correlation
         DepVars=[14 15]; %Left Frontop, Right Frontop
         
         
         for group=1:3 %1=pre, 2=term, 3=both
         
             if group~=3
                 if group==1
                     disp('Preterm');
                 else
                     disp('Term');
                 end

                 I=find([m.psychometr_msrs{2:end,2}]'==group);
                 
                Z=[];
                for ii=regressors
                    Z=[Z, [m.psychometr_msrs{min(I)+1:max(I)+1,ii}]'];
                end
                 
                 for jj=DepVars
                    if jj==14
                        disp('Left Frontop');
                    else
                        disp('Right Frontop');
                    end
                    X=[];
                    X=[[m.psychometr_msrs{min(I)+1:max(I)+1,jj}]',Z];

                    % RHO = partialcorr(X) returns the sample linear partial correlation
                    % coefficients between pairs of variables in X,
                    % controlling for the remaining variables in X.
                    
                    disp('Matrix: jj,GA,BW,Opti')
                    disp('Partial Correlation');
                    [RHO,PVAL]=partialcorr(X) 

                    % regstats(RESPONSES,DATA,MODEL) performs multiple
                    % regression
                    disp('Multiple regression');
                    RESPONSES=[m.psychometr_msrs{min(I)+1:max(I)+1,jj}]';
                    DATA=Z;
                    stats=regstats(RESPONSES,DATA);
                    pval=stats.tstat.pval
                    beta=stats.tstat.beta
                    model_p=stats.fstat.pval
                    
                    % Perform outlier diagnostics with regress
                    disp('Outlier Diagnostics');
                    
                    % Returns n-by-2 matrix rint of intervals that can be
                    % used to diagnose outliers. Id the interval rint(i,:)
                    % for observation i does not contain zero, the
                    % corresponding residual is larger than epected 95% of
                    % new observations, suggesting an outlier.
                    
                    % regress(y,X) on responses in y on predictions in X. y
                    % is an n-by-1 vector of observed responses
                    
                    y=[m.psychometr_msrs{min(I)+1:max(I)+1,jj}]';
                    X=[ones(size(Z,1),1) Z];
                    
                    [~,~,~,rint,~]=regress(y,X);
                    
                    outliers=all(rint'<0) + all(rint'>0);
                    
                    RESPONSES(find(outliers==1),:)=[];
                    DATA(find(outliers==1),:)=[];
                    
                    disp('After removal of outliers:')
                    disp('Multiple regression');
                    stats=regstats(RESPONSES,DATA);
                    pval=stats.tstat.pval
                    beta=stats.tstat.beta
                    model_p=stats.fstat.pval
                    
                 end
                 
             else
                 
                 disp('Large group');
                 
                 Z=[];
                 for ii=regressors
                    Z=[Z, [m.psychometr_msrs{2:end,ii}]'];
                 end

                 for jj=DepVars
                    if jj==14
                        disp('Left Frontop');
                    else
                        disp('Right Frontop');
                    end
                    X=[];
                    X=[[m.psychometr_msrs{2:end,jj}]',Z];

                    % RHO = partialcorr(X) returns the sample linear partial correlation
                    % coefficients between pairs of variables in X,
                    % controlling for the remaining variables in X.
                    
                    disp('Matrix: jj,GA,BW,Opti')
                    disp('Partial Correlation');
                    [RHO,PVAL]=partialcorr(X) 

                    % regstats(RESPONSES,DATA,MODEL) performs multiple
                    % regression
                    disp('Multiple regression');
                    RESPONSES=[m.psychometr_msrs{2:end,jj}]';
                    DATA=Z;
                    stats=regstats(RESPONSES,DATA);
                    pval=stats.tstat.pval
                    beta=stats.tstat.beta
                    model_p=stats.fstat.pval
                    
                    
                    y=RESPONSES;
                    X=[ones(size(DATA,1),1) DATA];
                    
                    [~,~,~,rint,~]=regress(y,X);
                    
                    outliers=all(rint'<0) + all(rint'>0);
                    
                    RESPONSES(find(outliers==1),:)=[];
                    DATA(find(outliers==1),:)=[];
                    
                    disp('After removal of outliers:')
                    disp('Multiple regression');
                    stats=regstats(RESPONSES,DATA);
                    pval=stats.tstat.pval
                    beta=stats.tstat.beta
                    model_p=stats.fstat.pval
                    
                 end
                 
             end
             
         end
         
end
