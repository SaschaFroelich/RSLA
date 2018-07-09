function [filename] = local_extrema_xcorr(varargin)
% This function creates a time delay matrix TD from a nifti file and stores it as a field in a matlab file (file format .mat).
% INPUT (in this order):
%   BOLDfile:   The nifti file containing 4D imaging data.
%   mask1:      Brainmask determining the rows of TD. Voxels with the same intensity will constitute one ROI.
%   mask2:      Brainmask determining the columns of TD. Voxels with the same intensity will constitute one ROI.
%   maxlag:     The maximum number of time-shifts around t=0 in the cross-correlation function. This defines a window [-maxlag,+maxlag] around t=0 in which the maximum of the cross-correlation function is looked for.
%   r_min:      The minimum absolute correlation value at t=0 for a pair of ROIs. If this criterion is not met, the corresponding entry in TD will contain NaN.
%   m_min:      The minimum absolute correlation value of the chosen maximum. If this criterion is not met, the corresponding entry in TD will contain NaN.
%   TR:         The TR value of the imaging sequence, in seconds.
%   minmaxdiscr:    Determines how the maximum in the cross-correlation function is defined. minmaxdiscr==2->zero-lag correlation, minmaxdiscr==3->most signficant
%   storage_dir (optional): Directory where to store new matlab file. If not specified, file will be stored in current directory.
%   tag (optional): Specific tag to add to name of new matlab file.
%   quitTF:     Boolean function to determine whether to close matlab after function terminates
%
% OUTPUT:
%   filename:   Name of the matlab-file containing time delay matrix TD.

BOLDfile=varargin{1};
mask1=varargin{2};
mask2=varargin{3};
maxlag=varargin{4};
r_min=varargin{5};
m_min=varargin{6};
TR=varargin{7};
minmaxdiscr=varargin{8}; %2->zero-lag correlation, minmaxdiscr==3->most signficant
quitTF=varargin{end};

if nargin==10
	storage_dir=varargin{9}; %where to store newly generated files. 
elseif nargin==11
	storage_dir=varargin{9}; %where to store newly generated files. 
	tag=varargin{10}; %Give the generated .mat-files simpler names to incorporate their analysis in a batch-script analysis with a group analysis following. tag is a string preceding the filenames
end

disp('Beginning Analysis.');
tic
c=clock;
disp(['Time is ' num2str(c(4)) 'hrs:' num2str(c(5)) 'mins on ' num2str(c(3)) 'th of ' num2str(c(2)) '.']);
%--------------------------------------------------------------------------%
%repetition time in ms
TR=TR*1000;
%--------------------------------------------------------------------------%

% coursematrix creates ROISxTIME MATRIX by averaging the BOLD timeseries for one ROI
% over all involved voxels. INTENS is a ROISx1 vector containing the
% intensities of the ROIS (from mask file). So INTENS(2,1) is the (mask-file)-intensity
% of the ROI which exhibits the signal timecourse ROISCOURSE(2,:) 

if 0
    %Check that ROIs in maskfile are of the form [2,3,4,5,6,...]
    if INTENS~=linspace(2,max(INTENS),length(INTENS))'
        error('ROIs in maskfile are not of the correct format (has to be 2,3,4,5....).');
    end
end

[ROISCOURSE1, INTENS1, numberofrois1] = coursematrix(BOLDfile, mask1);
[ROISCOURSE2, INTENS2, numberofrois2] = coursematrix(BOLDfile, mask2);

%mask no. 1 must have more ROIs than mask no.2
if numberofrois2>numberofrois1
    mask_temp1=mask1;
    mask_temp2=mask2;
    mask1=mask_temp2;
    mask2=mask_temp1;

    [ROISCOURSE1, INTENS1, numberofrois1] = coursematrix(BOLDfile, mask1);
    [ROISCOURSE2, INTENS2, numberofrois2] = coursematrix(BOLDfile, mask2);
end


disp(['Number of individual ROIS in first mask-file: ' num2str(numberofrois1) '']);
disp(['Number of individual ROIS in second mask-file: ' num2str(numberofrois2) '']);

if 0
    %Check that ROIs in maskfile are of the form [2,3,4,5,6,...]
    if isequal(INTENS1,linspace(2,max(INTENS1),length(INTENS1))') || isequal(INTENS2,linspace(2,max(INTENS2),length(INTENS2))')
        msgbox('ROIs in maskfile are not of the correct format (has to be 2,3,4,5....).','Error','error');
    end
end


disp('Create matrix LAGCORR');
toc
[LAGCORR,lags]=create_xcorr(ROISCOURSE1,ROISCOURSE2,r_min,maxlag);


%Now, compute the maximum of the lagged correlation by parabolic interpolation.
[rows, ~, pages] = size(LAGCORR);

%Create TD where the time delay will eventually be stored
TD=zeros(pages,rows);

%Using the same mask twice
if strcmp(mask1,mask2)
    disp('Generate TD matrix from one mask file.');
    %minmaxdiscr==2->zero-lag correlation, minmaxdiscr==3->most signficant
    if minmaxdiscr==2
        for page=1:pages
            for row=page:rows
                zerocorr=LAGCORR(row,maxlag+1,page);
                %If zero-lag is negative, look for minimun, else look for maximum
                if zerocorr>=0
                    [pks,locs]=findpeaks(LAGCORR(row,:,page));
                    [m,pos]=max(pks);
                elseif zerocorr<0
                    [~,locs]=findpeaks(-LAGCORR(row,:,page));
                    pks = LAGCORR(row,locs,page);
                    [m,pos]=min(pks);
                elseif isnan(zerocorr)
                    pks=[];
                end
                if isempty(pks)
                    TD(page,row)=NaN;
                    TD(row,page)=NaN;
                else
                   I=locs(pos); 

                   %prepare input for parabolic interpolation qint. If maximum is at
                   %timeshift maxlag or -maxlag the actual maximum cannot be
                   %deter_mined. Store NaN in TD.
                   if isnan(m)
                        TD(page,row)=NaN; 
                        TD(row,page)=NaN;
                   elseif abs(m) < m_min %abs(r) has to be larger than or equal to 0.116 in order to be p=0.05-significant
                        TD(page,row)=NaN;
                        TD(row,page)=NaN;
                   else
                    xm1=lags(I-1);
                    ym1=LAGCORR(row,I-1,page);
                    x0=lags(I);
                    y0=m;
                    xp1=lags(I+1);
                    yp1=LAGCORR(row,I+1,page);

                   %peak position in terms of TR (=timestep interval)
                    peak=qint(xm1,ym1,x0,y0,xp1,yp1);

                   %Enter the time-delay in terms ofxcorr(A,B) milliseconds into TD:
                   %Round the entries to the nearest multiple of 10ms
                   %TD(page,row) = round(peak*TR,-1);

                    %Round the entries to the nearest multiple of 50ms
                    temp=peak*TR-50*floor(peak*TR/50);
                    if temp>=25
                        temp=50;
                    elseif temp<25
                        temp=0;
                    end
                    new=temp+50*floor(peak*TR/50);
                    TD(page,row)=new;
                    TD(row,page)=-new;

                    %TD(1,2) > 0 means ROI 1 is BEHIND ROI 2 (2->1)

                   end
                end
            end
        end
    elseif minmaxdiscr==3 %Most significant
        for page=1:pages
            for row=page:rows
                %Look for minimum AND maximum, then determine which one is more
                %significant
                [pks_maximum,locs_maximum]=findpeaks(LAGCORR(row,:,page));
                [m_maximum,pos_maximum]=max(pks_maximum);

                [~,locs_minimum]=findpeaks(-LAGCORR(row,:,page));
                pks_minimum = LAGCORR(row,locs_minimum,page);
                [m_minimum,pos_minimum]=min(pks_minimum);

                if ~isempty(m_minimum) && ~isempty(m_maximum)
                    if abs(m_maximum) >= abs(m_minimum)
                        m=m_maximum;
                        locs=locs_maximum;
                        pos=pos_maximum;
                    elseif abs(m_maximum) < abs(m_minimum)
                        m=m_minimum;
                        locs=locs_minimum;
                        pos=pos_minimum;
                    end
                elseif isempty(m_minimum) && ~isempty(m_maximum)
                    m=m_maximum;
                    locs=locs_maximum;
                    pos=pos_maximum;
                elseif ~isempty(m_minimum) && isempty(m_maximum)
                    m=m_minimum;
                    locs=locs_minimum;
                    pos=pos_minimum;
                elseif isempty(m_minimum) && isempty(m_maximum)
                    m=[];
                    locs=[];
                    pos=[];
                end

                if isempty(m)
                    TD(page,row)=NaN;
                    TD(row,page)=NaN;
                else
                   I=locs(pos); 

                   %prepare input for parabolic interpolation qint. If maximum is at
                   %timeshift maxlag or -maxlag the actual maximum cannot be
                   %deter_mined. Store NaN in TD.
                   if isnan(m)
                        TD(page,row)=NaN; 
                        TD(row,page)=NaN;
                   elseif abs(m)<m_min %abs(r) has to be larger than or equal to 0.116 in order to be p=0.05-significant
                        TD(page,row)=NaN;
                        TD(row,page)=NaN;
                   else
                    xm1=lags(I-1);
                    ym1=LAGCORR(row,I-1,page);
                    x0=lags(I);
                    y0=m;
                    xp1=lags(I+1);
                    yp1=LAGCORR(row,I+1,page);

                   %peak position in terms of TR (=timestep interval)
                    peak=qint(xm1,ym1,x0,y0,xp1,yp1);

                   %Enter the time-delay in terms of milliseconds into TD:
                   %Round the entries to the nearest multiple of 10ms
                   %TD(page,row) = round(peak*TR,-1);

                    %Round the entries to the nearest multiple of 50ms
                    temp=peak*TR-50*floor(peak*TR/50);
                    if temp>=25
                        temp=50;
                    elseif temp<25
                        temp=0;
                    end
                    new=temp+50*floor(peak*TR/50);
                    TD(page,row)=new;
                    TD(row,page)=-new;

                    %TD(1,2) > 0 means ROI 1 is BEHIND ROI 2 (2->1)

                   end
                end
            end
        end
    end
else %Using two different masks
    disp('Generate TD matrix from two different mask-files.');
    %minmaxdiscr==2->zero-lag correlation, minmaxdiscr==3->most signficant
    if minmaxdiscr==2
        for page=1:pages
            for row=1:rows
                zerocorr = LAGCORR(row,maxlag+1,page);
                %If zero-lag is negative, look for minimun, else look for maximum
                if zerocorr>=0
                    [pks,locs]=findpeaks(LAGCORR(row,:,page));
                    [m,pos]=max(pks);
                elseif zerocorr<0
                    [~,locs]=findpeaks(-LAGCORR(row,:,page));
                    pks = LAGCORR(row,locs,page);
                    [m,pos]=min(pks);
                elseif isnan(zerocorr)
                    pks=[];
                end
                if isempty(pks)
                    TD(page,row)=NaN;
                else
                   I=locs(pos); 

                   %prepare input for parabolic interpolation qint. If maximum is at
                   %timeshift maxlag or -maxlag the actual maximum cannot be
                   %deter_mined. Store NaN in TD.
                   if isnan(m)
                        TD(page,row)=NaN; 
                   elseif abs(m) < m_min %abs(r) has to be larger than or equal to 0.116 in order to be p=0.05-significant
                        TD(page,row)=NaN;
                   else
                    xm1=lags(I-1);
                    ym1=LAGCORR(row,I-1,page);
                    x0=lags(I);
                    y0=m;
                    xp1=lags(I+1);
                    yp1=LAGCORR(row,I+1,page);

                   %peak position in terms of TR (=timestep interval)
                    peak = qint(xm1,ym1,x0,y0,xp1,yp1);

                   %Enter the time-delay in terms ofxcorr(A,B) milliseconds into TD:
                   %Round the entries to the nearest multiple of 10ms
                   %TD(page,row) = round(peak*TR,-1);

                    %Round the entries to the nearest multiple of 50ms
                    temp = peak*TR-50*floor(peak*TR/50);
                    if temp>=25
                        temp=50;
                    elseif temp<25
                        temp=0;
                    end
                    new = temp + 50*floor(peak*TR/50);
                    TD(page,row) = new;

                    %TD(1,2) > 0 means ROI 1 is BEHIND ROI 2 (2->1)

                   end
                end
            end
        end
    elseif minmaxdiscr==3 %Most significant
        for page=1:pages
            for row=1:rows
                %Look for minimum AND maximum, then determine which one is more
                %significant
                [pks_maximum,locs_maximum]=findpeaks(LAGCORR(row,:,page));
                [m_maximum,pos_maximum]=max(pks_maximum);

                [~,locs_minimum]=findpeaks(-LAGCORR(row,:,page));
                pks_minimum = LAGCORR(row,locs_minimum,page);
                [m_minimum,pos_minimum]=min(pks_minimum);

                if ~isempty(m_minimum) && ~isempty(m_maximum)
                    if abs(m_maximum) >= abs(m_minimum)
                        m=m_maximum;
                        locs=locs_maximum;
                        pos=pos_maximum;
                    elseif abs(m_maximum) < abs(m_minimum)
                        m=m_minimum;
                        locs=locs_minimum;
                        pos=pos_minimum;
                    end
                elseif isempty(m_minimum) && ~isempty(m_maximum)
                    m=m_maximum;
                    locs=locs_maximum;
                    pos=pos_maximum;
                elseif ~isempty(m_minimum) && isempty(m_maximum)
                    m=m_minimum;
                    locs=locs_minimum;
                    pos=pos_minimum;
                elseif isempty(m_minimum) && isempty(m_maximum)
                    m=[];
                    locs=[];
                    pos=[];
                end

                if isempty(m)
                    TD(page,row)=NaN;
                else
                   I=locs(pos); 

                   %prepare input for parabolic interpolation qint. If maximum is at
                   %timeshift maxlag or -maxlag the actual maximum cannot be
                   %deter_mined. Store NaN in TD.
                   if isnan(m)
                        TD(page,row)=NaN; 
                   elseif abs(m) < m_min %abs(r) has to be larger than or equal to 0.116 in order to be p=0.05-significant
                        TD(page,row)=NaN;
                   else
                    xm1=lags(I-1);
                    ym1=LAGCORR(row,I-1,page);
                    x0=lags(I);
                    y0=m;
                    xp1=lags(I+1);
                    yp1=LAGCORR(row,I+1,page);

                   %peak position in terms of TR (=timestep interval)
                    peak = qint(xm1,ym1,x0,y0,xp1,yp1);

                   %Enter the time-delay in terms of milliseconds into TD:
                   %Round the entries to the nearest multiple of 10ms
                   %TD(page,row) = round(peak*TR,-1);

                    %Round the entries to the nearest multiple of 50ms
                    temp = peak*TR-50*floor(peak*TR/50);
                    if temp>=25
                        temp=50;
                    elseif temp<25
                        temp=0;
                    end
                    new=temp+50*floor(peak*TR/50);
                    TD(page,row)=new;

                    %TD(1,2) > 0 means ROI 1 is BEHIND ROI 2 (2->1)

                   end
                end
            end
        end
    end
end

fprintf('\n \n \n');

noofNaN=length(find(isnan(TD))); %stores number of NaN entries in TD

disp('Write to file.');
%%%WRITE TD TO MAT-FILE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Determine Filename of generated file containing TD
splitfile = strsplit(BOLDfile,filesep);
filename=splitfile{end};
filename=filename(1:end-4); %-4 for nifti-files
if exist('tag')
    filename=strcat(filename,'_',tag);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BOLDfile=['' BOLDfile ''];

rois_in_mask1=INTENS1;
rois_in_mask2=INTENS2;
noofROIS_mask1=numberofrois1;
noofROIS_mask2=numberofrois2;

TDrois_rows=rois_in_mask1;
TDrois_cols=rois_in_mask2;

%Create a ROISx(1+TIMEPOINTS)-matrix to store the timecourses of the ROIS as
%well. The first column entries are the intensity values assigned by the
%mask, the following columns are the time-series:
ROIS_TIMECOURSES1=[INTENS1, ROISCOURSE1];
ROIS_TIMECOURSES2=[INTENS2, ROISCOURSE2];

Method='Local Extrema. abs(r(t=0))>r_min. ';

if minmaxdiscr==2
    Method=strcat(Method,'Chosen extremum determined by zero-lag correlation. ');
elseif minmaxdiscr==3
    Method=strcat(Method,'Chosen extremum determined by significance. ');
end

Method=strcat(Method,'Used algorithm: xcorr ');

%Make TD=-TD, such that TD(2,1)<0 if 1->2. This way, Column Means is the actual lag projection map. 
TD=TD*(-1);

comment='TD*(-1) in order to make ColumnMeans the actual lag projection map.';

disp(['No of NaN in TD:' num2str(length(find(isnan(TD)))) '']);

ColumnMeans=nanmean(TD);
RowMeans=nanmean(TD');

%   If *Means only has NaN in a column, nanmean of that column will return NaN -> replace by 0
ColumnMeans(find(isnan(ColumnMeans)))=0;
RowMeans(find(isnan(RowMeans)))=0;

%   Compute ZL (zero-lag correlation matrix)
%   ZL=zeros(noofROIS_mask1,noofROIS_mask2);

for row=1:noofROIS_mask1
	for col=1:noofROIS_mask2
        ZL(row,col)=corr2(ROIS_TIMECOURSES1(row,2:end),ROIS_TIMECOURSES2(col,2:end));
	end
end

% First level (single subject) or second level (group)
level=1;

if ~exist('storage_dir') %no storage_dir
	save(['' pwd filename '.mat'],'BOLDfile','mask1','mask2','rois_in_mask1','rois_in_mask2','lags','noofNaN','noofROIS_mask1','noofROIS_mask2','Method','TD','ROIS_TIMECOURSES1','ROIS_TIMECOURSES2','r_min','m_min','TDrois_rows','TDrois_cols','comment','ColumnMeans','RowMeans','level');
    fisp(['File stored ad ' pwd filename '.mat']);
elseif exist('storage_dir')
	save(['' storage_dir '' filename '.mat'],'BOLDfile','mask1','mask2','rois_in_mask1','rois_in_mask2','lags','noofNaN','noofROIS_mask1','noofROIS_mask2','Method','TD','ROIS_TIMECOURSES1','ROIS_TIMECOURSES2','r_min','m_min','TDrois_rows','TDrois_cols','comment','ColumnMeans','RowMeans','level');
    disp(['File saved as ' storage_dir '' filename '.mat']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
disp('Done!');

if quitTF
    quit
end

end
