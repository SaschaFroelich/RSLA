function [LAGCORR,lags] = create_xcorr(varargin)
% Input:
%       ROISCOURSE1: ROIsx(timesteps+1) matrix. Each row vector in ROISCOURSE1 is
%                    the time-series of a ROI in the first maskfile. 
%                    NB: First column in ROISCOURSE contains intensity values 
%                    of ROIs. Timecourses follow in columns 2:end.
%       ROISCOURSE2: Same as above, for second maskfile.
%       r_min:       The minimum absolute correlation value at t=0 for a pair of ROIs. If this criterion is not met, the corresponding entry in TD will contain NaN.
%       maxlag:      The maximum number of time-shifts around t=0 in the cross-correlation function. This defines a window [-maxlag,+maxlag] around t=0 in which the maximum of the cross-correlation function is looked for.
% Output:
%       LAGCORR:     A threedimendional matrix. Every row contains the
%                    cross-correlation values of the BOLD timeseries
%                    in ROISCOURSE1 with the ones in ROISCOURSE2
%       lags:        Vector of same length as rows in LAGCORR. Specifies
%                    the timeshifts for the cross-correlation values in
%                    LAGCORR.
%
ROISCOURSE1=varargin{1};
ROISCOURSE2=varargin{2};
r_min=varargin{3};
maxlag=varargin{4};

[rois1, timepoints1]=size(ROISCOURSE1);
[rois2, timepoints2]=size(ROISCOURSE2);

if timepoints1~=timepoints2
    msgbox('Number of timepoints are not the same. Implementation error.','Error','error');
end

%Make sure that ROISCOURSE1 is the matrix with most rois.
if rois1<rois2
   mtr1=ROISCOURSE1;
   mtr2=ROISCOURSE2;

   clear ROISCOURSE1
   clear ROISCOURSE2

   ROISCOURSE2=mtr1;
   ROISCOURSE1=mtr2;

   clear mtr1
   clear mtr2

[rois1, timepoints1]=size(ROISCOURSE1);
[rois2, timepoints2]=size(ROISCOURSE2);

end

LAGCORR_temp=zeros(1,2*maxlag-1);
LAGCORR = zeros(rois2,2*maxlag+1,rois1);

[throwaway,lags]=xcorr(ROISCOURSE1(1,:),ROISCOURSE2(1,:));
I=find(lags==0);

for page=1:rois1
    for row = 1:rois2
        %Compute correlation at lag zero
        zerocorr = corr2(ROISCOURSE1(page,:),ROISCOURSE2(row,:));
        if abs(zerocorr) > r_min
            LAGCORR_temp=xcorr(ROISCOURSE1(page,:),ROISCOURSE2(row,:));
            LAGCORR(row,:,page)=LAGCORR_temp(I-maxlag:I+maxlag);
        else
            LAGCORR(row,:,page)=NaN;
        end
    end
end

lags=lags(I-maxlag:I+maxlag);


end