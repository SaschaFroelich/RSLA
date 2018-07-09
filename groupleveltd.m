function [] = groupleveltd(varargin)
% groupleveltd creates a .mat-file containing a group level time delay
% matrix TD by averaging over all TDs in the input files.
% INPUT (in this order):
%       filenames:          cell array containing filenames. Each entry corresponds to a mat-file containing a time delay matrix TD.
%       filepath:           Path to files in cell array "filenames".
%       name:               name of mat file containing group level TD.
%       strg_dir:           storage directory of mat file containing group level TD. The new file will thus be stored as strcat(strg_dir,name)
%       quitTF (optional):  Boolean function to determine whether to close matlab 
%                           after function terminates.
%                           DEFAULT=1

fprintf('\n');
disp('COMPUTING group-level TD (groupleveltd.m)');

filenames=varargin{1};
filepath=varargin{2};
name=varargin{3};
strg_dir=varargin{4};
if nargin==5
    quitTF=varargin{5};
else
    quitTF=1;
end
    
if ~strcmp(strg_dir(end),filesep)
  strg_dir=strcat(strg_dir,filesep); 
end

newfile=strcat(strg_dir,name);
        
for iter=1:length(filenames)
   file = strcat(filepath,filenames{iter});
   m=load(['' file '']);

    if iter==1
        TDrois_rows=m.TDrois_rows;
        TDrois_cols=m.TDrois_cols;
    else
        TDrois_rows_temp=m.TDrois_rows;
        TDrois_cols_temp=m.TDrois_cols;


        if (length(TDrois_rows_temp)~=length(TDrois_rows)) || (length(TDrois_cols_temp)~=length(TDrois_cols))
            error('TD-matrices do not have same dimensions across all files. Cannot perform grouplagproj.m'); 
        end

        TDrois_rows=TDrois_rows_temp;
        TDrois_cols=TDrois_cols_temp;
    end

    %Dimensions of groupTD
    groupTDrois_rows=m.TDrois_rows; 
    groupTDrois_cols=m.TDrois_cols;

end
        
% divider counts the no. of non-NaN-entries for every TD-entry
% across all files
divider=zeros(length(groupTDrois_rows),length(groupTDrois_cols));
disp('Now creating nanmean-divider.');
for iter=1:length(filenames)
    disp(['File ' num2str(iter) '']);
    file = strcat(filepath,filenames{iter});
    m=load(file);

    % Retrieve information for the save-function at the end
    r_min=m.r_min;
    m_min=m.m_min;
    mask1=m.mask1;
    mask2=m.mask2;

    % Improvise a nanmean method without storing all TD-matrices
    for row=1:length(m.TDrois_rows);
        for col=1:length(m.TDrois_cols);
            if ~isnan(m.TD(row,col))
                divider(row,col)=divider(row,col)+1;
            end
        end
    end

end

groupTD=zeros(length(groupTDrois_rows),length(groupTDrois_cols));
% groupZL=zeros(length(groupTDrois_rows),length(groupTDrois_cols));
disp('Now computing groupTD.');
for iter=1:length(filenames)
    disp(['File ' num2str(iter) '']);
    file = strcat(filepath,filenames{iter});
    m=load(['' file '']);

   % testzl(:,:,iter)=m.ZL;
    for row=1:length(m.TDrois_rows);
        for col=1:length(m.TDrois_cols);
            if ~isnan(m.TD(row,col)) && divider(row,col)~=0
                groupTD(row,col)=groupTD(row,col)+m.TD(row,col)/divider(row,col);
            elseif isnan(m.TD(row,col)) && divider(row,col)==0 %This has to be an && statement, as divider==0 means that ALL matrices contain a NaN-entry for that particular position
                groupTD(row,col)=0;
            end
        end
    end
end

TD=groupTD;
noofNaN=0;
noofROIS_mask1=m.noofROIS_mask1;
noofROIS_mask2=m.noofROIS_mask2;
% First level (single subject) or second level (group)
level=2;

filenames=filenames';
pathtofiles=filepath;

description='File contains group-level TD from 2-mask files.';
lags=m.lags;

save(newfile,'TD','pathtofiles','filenames','r_min','m_min','mask1','mask2','description','TDrois_rows','noofNaN','noofROIS_mask1','noofROIS_mask2','lags','level');

disp([' Storing group level file as ' newfile '']);
disp('Done!');

if quitTF
    quit
end

end
