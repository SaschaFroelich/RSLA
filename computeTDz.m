function [TDz,ColumnMeans,RowMeans] = computeTDz(TD)
% computeTDz() conputes TDz from TD.
%   INPUT:
%       TD:             Matrix TD
%   OUTPUT:
%       TDz:            Zero-centered version of TD.
%       ColumnMeans:    Vector containing the arithmetic means of the 
%                       columns in matrix TD.
%       RowMeans:       Vector containing the arithmetic means of the 
%                       rows in matrix TD.

	ColumnMeans=nanmean(TD);
    ColumnMeans(find(isnan(ColumnMeans)))=0;
    RowMeans=nanmean(TD');
    RowMeans(find(isnan(RowMeans)))=0;
    TDz=bsxfun(@minus,TD,nanmean(TD));
    % If TD has only NaN in a column, nanmean of that column will return NaN
    % -> replace by 0
    TDz(find(isnan(TDz)))=0;
    % If TD contains NaN-entries, so will now TDz, as NaN-numbrt=NaN

end