function [] = grouplagproj(varargin)
% grouplagproj computes the 1D group level lag projection, which corresponds
% to the average ColumnMeans vector of all involved TD matrixes. 1D
% group-level lag projection is NOT the ColumnMeans vector of the groupTD
% matrix (according to Mitra).
% INPUT:
%       filename:           The group level matlab file.
%       filepath:           Path to the group level matlab file (with trailing file
%                           separator)
%       quitTF (optional):  Boolean function to determine whether to close matlab 
%                           after function terminates.
%                           DEFAULT=1
disp('FUNCTION GROUPLAGPROJ()');
disp('Preparing for group-level Lag Projection computation.');

%"matfile" contains groupTD
matfile=strcat(varargin{2},varargin{1});
m=load(matfile);
    
if ~isfield(m,'pathtofiles')
	error('Path to Files not specified');
end

if nargin==3
    quitTF=varargin{3};
else
   quitTF=1; 
end
   
%Create a filesxrois matrix CLMNS, containing as rows the ColumnMeans of each
%file while each ColumnMeans only contains the ROIs that are contained
%inside groupTD.


%Only so far implemented for com-processed files
[rows,columns]=size(m.TD);
%Both RWS and CLMN are stored as row vectors
CLMNS=zeros(length(m.filenames),columns);
RWS=zeros(length(m.filenames),rows);

for iter=1:length(m.filenames)
	file=strcat(m.pathtofiles,m.filenames{iter});
	m2=load(file);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Check if all involved single-subject matrices have same number of
    %ROIs in rows and columns
    
    if iter==1
        TDrois_rows=m2.TDrois_rows;
        TDrois_cols=m2.TDrois_cols;
    else
        TDrois_rows_temp=m2.TDrois_rows;
        TDrois_cols_temp=m2.TDrois_cols;
        
        if length(TDrois_rows_temp)~=length(TDrois_rows) || length(TDrois_cols_temp)~=length(TDrois_cols)
            error('TD-matrices do not have same dimensions across all files. Cannot perform grouplagproj.m'); 
        end
        
        TDrois_rows=TDrois_rows_temp;
        TDrois_cols=TDrois_cols_temp;
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Both RWS and CLMN are stored as row vectors
    if isfield(m2,'ColumnMeans')
        CLMNS(iter,:)=m2.ColumnMeans;
        RWS(iter,:)=m2.RowMeans;
    else
        CLMNS(iter,:)=nanmean(m2.TD);
        RWS(iter,:)=nanmean(m2.TD');  
    end
end

    
lag_projection_clmns=nanmean(CLMNS);
lag_projection_rws=nanmean(RWS);

GroupLagProjection=lag_projection_clmns;

disp('Saving group level lag projection as ''GroupLagProjection''');
save(matfile,'GroupLagProjection','-append');

disp('Computation of group-level lag projection terminated.');

if quitTF
    quit
end

end
