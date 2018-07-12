function []=combineCSV()
% Combines CSV-files into one by appending the second csv-file to the
% first. IMPORTANT: Matching of the rows is done by comparing the 
% contents of the first column (rows 2 to end) of the first file with the
% contenst of the first column (rows 2 to end) of the second file. If two 
% cells of the first column contain the same entry in both files, then
% they are matched and appended. If no cell matching can be done for a cell
% in first column of second file, a new row will be created for the
% corresponding cell.
%   INPUT:
%       none
%   OUTPUT:
%       none

[csvfile1Name,csvfile1Path]=uigetfile(['' pwd filesep '*.csv'],'MultiSelect','Off','Select .csv-file containing behavioral scores (each subject occupies one row.).');
csvfile1=strcat(csvfile1Path,csvfile1Name)

fid = fopen(['' csvfile1 '']); % 'Data' worksheet

clear contents1
contents1=cell(1,1);
contents1Temp=cell(1,1);
EscLoop=0;
line=0;
while EscLoop~=1
    line=line+1;
    if line==1
        contents1Temp{line,1}=fgetl(fid);
    else
        contents1Temp=[contents1Temp; fgetl(fid)];
    end
    if contents1Temp{line,1}==-1
       EscLoop=1;
       contents1=contents1Temp(1:end-1,:);
       clear contents1Temp
    end
end

acols = length(find(contents1{1,1}==','))+1; % number of columns
aformat = repmat('%s ', 1, acols); % create format
fclose(fid);

fid = fopen(['' csvfile1 '']);
ColumnsInFile1 = textscan(fid,aformat,'Delimiter',',');
fclose(fid);      

no_lines1=size(contents1,1);

[csvfile2Name,csvfile2Path]=uigetfile(['' pwd filesep '*.csv'],'MultiSelect','Off','Select second .csv-file, to be added to first one (each subject occupies one row.).');
csvfile2=strcat(csvfile2Path,csvfile2Name)

fid = fopen(['' csvfile2 '']); % 'Data' worksheet

clear contents2
contents2=cell(1,1);
contents2Temp=cell(1,1);
EscLoop=0;
line2=0;
while EscLoop~=1
    line2=line2+1;
    if line2==1
        contents2Temp{line2,1}=fgetl(fid);
    else
        contents2Temp=[contents2Temp; fgetl(fid)];
    end
    if contents2Temp{line2,1}==-1
       EscLoop=1;
       contents2=contents2Temp(1:end-1,:);
       clear contents2Temp
    end
end

acols = length(find(contents2{1,1}==','))+1; % number of columns
aformat = repmat('%s ', 1, acols); % create format
fclose(fid);

fid = fopen(['' csvfile2 '']);
ColumnsInFile2 = textscan(fid,aformat,'Delimiter',',');
fclose(fid);   

no_lines2=size(contents2,1);


if no_lines1 ~= no_lines2
   % error('The two csv-files differ in number of rows.');
end


SubjMatch=zeros(no_lines1,2);

FirstColFl1=ColumnsInFile1{1};
FirstColFl2=ColumnsInFile2{1};
for ii=2:no_lines1
    SubjMatch(ii,1)=ii;
    if isempty(find(~cellfun(@isempty,strfind(FirstColFl2,FirstColFl1{ii}))))   
        disp(['' FirstColFl1{ii} ' not found in second csv-file.']);
        SubjMatch(ii,2)=-1;
    else
        SubjMatch(ii,2)=find(~cellfun(@isempty,strfind(FirstColFl2,FirstColFl1{ii})));
    end
end


DIRECTORYNAME=uigetdir(pwd,'Where to store combined csv-file (will be called ''combi.csv'')?');

fid=fopen(['' DIRECTORYNAME '' filesep 'combinecsv_' csvfile1Name csvfile2Name],'w');
% Print first row
fprintf(fid,'%s,%s\n',contents1{1,1},contents2{1,1});

% Determine no. of cells added by second file
NoCells=length(strsplit(contents2{1,:},',','CollapseDelimiters',false));
Placeholder='NaN';

for cl=2:NoCells
    Placeholder=strcat(Placeholder,',NaN');
end

for ii=2:no_lines1
    if SubjMatch(ii,2) ~= -1
        fprintf(fid,'%s,%s\n',contents1{ii,1},contents2{SubjMatch(ii,2),1});
    else
        fprintf(fid,'%s,%s\n',contents1{ii,1},Placeholder);
    end
end

fclose(fid);

disp('csv combination Done.');
end