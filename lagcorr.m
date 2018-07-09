function [LAGCORR] = lagcorr(file,maxlag,fig,usexcorr,roi1,roi2)
% Function lagcorr() computes the lagged cross-correlation function of BOLD
% timecourses of two given ROIs.
% INPUT:
%       file:       matlab file including matrices ROIS_TIMECOURSES1 and 
%                   ROIS_TIMECOURSES2 (see below for explanation).
%       maxlag:     The maximum number of timesteps (in terms of TR) around 
%                   t=0 within which a maximum for time delay determination 
%                   is sought. 
%       fig:        Boolean function whether or not to open a figure window 
%                   for with plots of the correlations between roi1 and roi2.
%       usexcor:    Boolean function whether or not to use MATLAB function 
%                   xcorr(). Should always be 1.
%       roi1:       First roi which is used to compute the cross-correlation.
%       roi2:       Second roi which is used to compute the cross-correlation.
% OUTPUT:
%       LAGCORR:    2x(2*maxlag+1) matrix. The first row contains the lagged
%                   cross-correlation of roi1 with roi2. The second row 
%                   contains the lagged cross-correlation of roi2 with roi1.
%
% FURTHER INFO
%   ROIS_TIMECOURSES1 is a nx(m+1) matrix. Each row corresponds to the n
%   ROIs as defined by mask1. The first column of each row is that ROIs
%   BOLD timecourse (m timesteps). Same for ROIS_TIMECOURSES2, which corresponds to the
%   ROIs as defined by mask2.
%
m=load(['' file '']);
if isfield(m,'mask')
    INTENS=m.ROIS_TIMECOURSES(:,1);
    TIMECOURSES=m.ROIS_TIMECOURSES(:,2:end);
elseif isfield(m,'mask2')
    INTENS1=m.ROIS_TIMECOURSES1(:,1);
    TIMECOURSES1=m.ROIS_TIMECOURSES1(:,2:end);
    INTENS2=m.ROIS_TIMECOURSES2(:,1);
    TIMECOURSES2=m.ROIS_TIMECOURSES2(:,2:end);
end

%sort varargin with the ROI with the smallest index standing at varargin(1)
ROIs=[roi1, roi2];

LAGCORR = zeros(2,2*maxlag+1);


%when the correlation is computed at lags(lag)=-5 it means that the
%second timecourse was shifted five timesteps to the left.
%a,b as in corr2(a,b)
if usexcorr
    if ~isfield(m,'mask2')
        pos1=find(INTENS==roi1);
        pos2=find(INTENS==roi2);
        [LAGCORR_temp1,lags_temp]=xcorr(TIMECOURSES(pos1,:),TIMECOURSES(pos2,:));
        LAGCORR_temp2=xcorr(TIMECOURSES(pos2,:),TIMECOURSES(pos1,:));
    elseif isfield(m,'mask2')
        pos1=find(INTENS1==roi1);
        pos2=find(INTENS2==roi2);
        [LAGCORR_temp1,lags_temp]=xcorr(TIMECOURSES1(pos1,:),TIMECOURSES2(pos2,:));
        LAGCORR_temp2=xcorr(TIMECOURSES2(pos2,:),TIMECOURSES1(pos1,:));
    end
    I=find(lags_temp==0);
    LAGCORR(1,:)=LAGCORR_temp1(I-maxlag:I+maxlag);
    LAGCORR(2,:)=LAGCORR_temp2(I-maxlag:I+maxlag);
    lags=lags_temp(I-maxlag:I+maxlag);
else
    lags=linspace(-maxlag,maxlag,2*maxlag+1);
    for lag = 1:length(lags)
        %Now we have the intensity of the ROIs to be corrrelated,
        %determine which position a and b have in TIMECOURSES
        if ~isfield(m,'mask2')
            pos1=find(INTENS==roi1);
            pos2=find(INTENS==roi2);
            LAGCORR(1,lag)= corr2(TIMECOURSES(pos1,maxlag+1:end-maxlag),TIMECOURSES(pos2,2*maxlag+2-lag:end+1-lag));
            LAGCORR(2,lag)= corr2(TIMECOURSES(pos2,maxlag+1:end-maxlag),TIMECOURSES(pos1,2*maxlag+2-lag:end+1-lag));
        elseif isfield(m,'mask2')
            pos1=find(INTENS1==roi1);
            pos2=find(INTENS2==roi2);
            LAGCORR(1,lag)= corr2(TIMECOURSES1(pos1,maxlag+1:end-maxlag),TIMECOURSES2(pos2,2*maxlag+2-lag:end+1-lag));
            LAGCORR(2,lag)= corr2(TIMECOURSES2(pos2,maxlag+1:end-maxlag),TIMECOURSES1(pos1,2*maxlag+2-lag:end+1-lag));
        end
    end
end
if fig
    figure;
    plot(lags,LAGCORR(1,:));
    xlabel('lag in TIMESHIFTS (not sec)');
    set(gca,'xtick',lags(1):1:lags(end)); 
    title(['lagged corr between roi1 (intens= ' num2str(roi1) ') and roi2 (intens= ' num2str(roi2) ')']);
    grid on;

    figure;
    plot(lags,LAGCORR(2,:));
    xlabel('lag in TIMESHIFTS (not sec)');
    set(gca,'xtick',lags(1):1:lags(end)); 
    title(['lagged corr between roi1 (intens= ' num2str(roi1) ') and roi2 (intens= ' num2str(roi2) ')']);
    grid on;
end

%i,j stand for the ROIs in LAGCORR.
%varargin(i) stands for the row, in which the ROI's timecourse was saved in
%the huge ROIS_TIMECOURSES matrix. 
%INTENS(i) stands for the intensity of the ROI in the mask file.
    
rois=ROIs;
%disp('Done.');

end