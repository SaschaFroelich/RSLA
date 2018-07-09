function varargout = rslastarter(varargin)
% RSLASTARTER MATLAB code for rslastarter.fig
%      RSLASTARTER, by itself, creates a new RSLASTARTER or raises the existing
%      singleton*.
%
%      H = RSLASTARTER returns the handle to a new RSLASTARTER or the handle to
%      the existing singleton*.
%
%      RSLASTARTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RSLASTARTER.M with the given input arguments.
%
%      RSLASTARTER('Property','Value',...) creates a new RSLASTARTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rslastarter_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rslastarter_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rslastarter

% Last Modified by GUIDE v2.5 23-Apr-2018 12:35:45

% AUTHOR: Sascha Frölich, sascha.froelich@gmail.com

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rslastarter_OpeningFcn, ...
                   'gui_OutputFcn',  @rslastarter_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before rslastarter is made visible.
function rslastarter_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rslastarter (see VARARGIN)

% Choose default command line output for rslastarter
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

handles=guidata(hObject);
handles.rmin=str2double(get(handles.r_min,'String')); %handles.r_min is the text field, handles.rmin is the value

contents=cellstr(get(handles.m_min,'String'));
handles.mmin=contents{get(handles.m_min,'Value')}; %handles.m_min is the drop-down menu, handles.mmin is the value
guidata(hObject,handles);

% UIWAIT makes rslastarter wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = rslastarter_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function boldfile_descr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to boldfile_descr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes during object creation, after setting all properties.
function text4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% --- Executes on button press in RUN.
function RUN_Callback(hObject, eventdata, handles)
% hObject    handle to RUN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set(handles.RUN,'Enable','off');
contents=cellstr(get(handles.m_min,'String'));
sign_level=contents{get(handles.m_min,'Value')};
usexcorr=1; %If set to 1, use Matlab xcorr method (as suggested by Mitra)

if ~get(handles.boldfile_descr,'value') || ~get(handles.text4,'value') || (get(handles.use_second_mask,'Value') && ~get(handles.second_mask,'Value'))
    msgbox('Choose Boldfile/Mask.','Error','error');
elseif ~isfield(handles,'TR') || (isfield(handles,'TR') && isnan(str2double(get(handles.TR,'string'))))
    msgbox('Parameters not valid (field TR).','Error','error'); 
elseif isnan(handles.lags) || isnan(handles.rmin)
    msgbox('Parameters not valid (Either field "Lags" of field "rmin").','Error','error'); 
elseif strcmp(sign_level,'Lagged Corr. p')
    msgbox('Parameters not valid.','Error','error'); 
else
    set(handles.chooselag,'Enable','Off');
    set(handles.r_min,'Enable','Off');
    set(handles.m_min,'Enable','Off');
    set(handles.RUN,'Enable','Off');
    set(handles.TR,'Enable','Off');
    drawnow;
    filename=strsplit(get(handles.boldfile_descr,'string'),' ');
    filename=filename{end};
    maskname=strsplit(get(handles.text4,'string'),' ');
    maskname=maskname{end};
    filepath=strsplit(get(handles.boldfile_descr,'Tooltipstring'),' ');
    filepath=filepath{end};
    maskpath=strsplit(get(handles.text4,'Tooltipstring'),' ');
    maskpath=maskpath{end};
    file=strcat(filepath,filename);
    mask=strcat(maskpath,maskname);

    fileinfo=nii_read_header(file);
    maskinfo=nii_read_header(mask);
    maskpixdim=maskinfo.PixelDimensions;
    filepixdim=fileinfo.PixelDimensions;

    %maskqform=maskinfo.QformCode
    %masksform=maskinfo.SformCode
    %fileqform=fileinfo.QformCode
    %filesform=fileinfo.QformCode
    
    BOLDfile=file;
    maxlag=str2double(get(handles.chooselag,'string'));
    r_min=str2double(get(handles.r_min,'string'));

    contents=cellstr(get(handles.m_min,'String'));
    sign_level=contents{get(handles.m_min,'Value')};
    sign_level=str2double(num2str(sign_level));
    if sign_level == 0.01
            m_min = 0.15;
    elseif sign_level== 0.05
            m_min= 0.113;
    elseif sign_level== 0.1
            m_min=0.095;
    elseif sign_level== 0.2
            m_min=0.075;
    elseif sign_level== 1
            m_min=0;
    end
    
    if maskpixdim(1) ~= filepixdim(1)
        msgbox('File and Mask are of conflicting orientations.','Error','error'); 
    else
        minmaxdiscr=2;
        maxlag
        r_min
        m_min
        minmaxdiscr
        usexcorr
        if get(handles.use_second_mask,'Value')
            newfile=local_extrema_xcorr(BOLDfile,mask,handles.secondmask,maxlag,r_min,m_min,str2double(get(handles.TR,'string')),minmaxdiscr,handles.strg_dir,0);
        else
            newfile=local_extrema_xcorr(BOLDfile,mask,mask,maxlag,r_min,m_min,str2double(get(handles.TR,'string')),minmaxdiscr,handles.strg_dir,0);
        end
        msgbox(['Saved new file as ' newfile '']);        
        disp(['Saved new file as ' newfile '']);
        set(handles.newanalysis,'Visible','On');
        set(handles.examinefile,'Visible','On');
    end
    
end
set(handles.chooselag,'Enable','On');
set(handles.r_min,'Enable','On');
set(handles.m_min,'Enable','On');
set(handles.RUN,'Enable','On');
drawnow;

function chooselag_Callback(hObject, eventdata, handles)
% hObject    handle to chooselag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of chooselag as text
%        str2double(get(hObject,'String')) returns contents of chooselag as a double

lag=str2double(get(handles.chooselag,'string'));

handles=guidata(hObject);
handles.lags=lag;
guidata(hObject,handles);

if ~isnan(lag)
    set(handles.LagsText,'string',get(handles.chooselag,'string'));
else
    set(handles.LagsText,'string','NaN');
end

% --- Executes during object creation, after setting all properties.
function chooselag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chooselag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadfile.
function loadfile_Callback(hObject, eventdata, handles)
% hObject    handle to loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.boldfile_descr,'Value',0.0);

if isfield(handles,'batch') && handles.batch % Create batch for a number of files
        if get(handles.tdgen,'Value') || get(handles.tdgen_group,'Value') || get(handles.gsr,'Value')
            if isfield(handles,'Recentfilepath')
                [BOLDfilename,BOLDfilepath]=uigetfile(['' handles.Recentfilepath '*.nii'],'MultiSelect','On','Select nifti-file(s)');
            else
                [BOLDfilename,BOLDfilepath]=uigetfile(['' pwd filesep '*.nii'],'MultiSelect','On','Select nifti-file(s)');
            end

            if strfind(BOLDfilepath,' ');
               msgbox('Make sure there are no space characters in the BOLDfile directory path.');
               error('Make sure there are no space characters in the BOLDfile directory path.'); 
            end
            
            if ~isequal(BOLDfilepath,0)
                handles=guidata(hObject);
                handles.Recentfilepath=BOLDfilepath;
                guidata(hObject,handles);
            else
                return
            end
        else
            [BOLDfilename,BOLDfilepath]=uigetfile(['' pwd filesep '*.mat'],'MultiSelect','On','Select mat-file(s)');
        end
else ~isfield(handles,'batch') || ~handles.batch % one file to select
    [BOLDfilename,BOLDfilepath]=uigetfile(['' pwd filesep '*.nii'],'Select nifti-file');
end

if ~isempty(BOLDfilename)
    if iscell(BOLDfilename) %several files selected, BATCH creator active
        set(handles.boldfile_descr,'string',['' num2str(length(BOLDfilename)) ' files selected.']); 
        GSRInfoString=cell(length(BOLDfilename),1);
        for file=1:length(BOLDfilename)
            GSRInfoString{file}=strcat(BOLDfilename{file},' : no brainmask selected');
        end
        set(handles.GSRBrainmasksInfo,'string',GSRInfoString); 
    elseif ~iscell(BOLDfilename) %one file selected
        BOLDfile=strcat(BOLDfilepath,BOLDfilename);
        set(handles.boldfile_descr,'String',['File: ' BOLDfilename '']);
    end
    
    set(handles.parallelize,'Value',0.0);
    set(handles.boldfile_descr,'Value',1.0);

    set(handles.boldfile_descr,'TooltipString',['Path: ' BOLDfilepath '']);
    handles=guidata(hObject);
    handles.BOLDfilename=BOLDfilename;
    handles.BOLDfilepath=BOLDfilepath;
    guidata(hObject,handles);
end

% --- Executes on button press in loadmask.
function loadmask_Callback(hObject, eventdata, handles)
% hObject    handle to loadmask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text4,'Value',0.0)
ThisFile=mfilename('fullpath');
path=strcat(ThisFile(1:end-length(mfilename)),'masks',filesep);
[maskname,maskpath]=uigetfile(['' path '*.nii']);

if strfind(maskpath,' ');
   msgbox('Make sure there are no space characters in the mask''s directory.');
   error('Make sure there are no space characters in the mask''s directory.'); 
end

clear ThisFile path
mask=strcat(maskpath,maskname);
set(handles.text4,'String',['Mask: ' maskname ''])
set(handles.text4,'TooltipString',['Path: ' maskpath '']);
if exist(mask)
    set(handles.text4,'Value',1.0);
    
    handles=guidata(hObject);
    handles.mask=mask;
    guidata(hObject,handles);
    
end


function r_min_Callback(hObject, eventdata, handles)
% hObject    handle to r_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of r_min as text
%        str2double(get(hObject,'String')) returns contents of r_min as a double
r_min=str2double(get(handles.r_min,'String'));

handles=guidata(hObject);
handles.rmin=r_min;
guidata(hObject,handles);

if ~isnan(r_min)
    set(handles.RText,'string',get(handles.r_min,'string'));
else
    set(handles.RText,'string','NaN');
end


% --- Executes during object creation, after setting all properties.
function r_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in m_min.
function m_min_Callback(hObject, eventdata, handles)
% hObject    handle to m_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns m_min contents as cell array
%        contents{get(hObject,'Value')} returns selected item from m_min
contents=cellstr(get(handles.m_min,'String'));
sign_level=contents{get(handles.m_min,'Value')};

if strcmp(sign_level,'Lagged Corr. p')
    
else
    sign_level=str2num(num2str(sign_level));
    if sign_level == 0.01
            m_min = 0.15;
    elseif sign_level== 0.05
            m_min= 0.113;
    elseif sign_level== 0.1
            m_min=0.095;
    elseif sign_level== 0.2
            m_min=0.075;
    elseif sign_level== 1
            m_min=0;
    end
end

handles=guidata(hObject);
handles.mmin=m_min;
guidata(hObject,handles);

set(handles.MText,'string',sign_level);

% --- Executes during object creation, after setting all properties.
function m_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to m_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in newanalysis.
function newanalysis_Callback(hObject, eventdata, handles)
% hObject    handle to newanalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.newanalysis,'Visible','off');
set(handles.RUN,'visible','on');
set(handles.chooselag,'Enable','On');
set(handles.r_min,'Enable','On');
set(handles.m_min,'Enable','On');
set(handles.RUN,'Enable','On');
set(handles.examinefile,'Visible','Off');

% --- Executes on button press in examinefile.
function examinefile_Callback(hObject, eventdata, handles)
% hObject    handle to examinefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
CorrViewer;

% --- Executes on button press in activate_batch.
function activate_batch_Callback(hObject, eventdata, handles)
% hObject    handle to activate_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of activate_batch

if get(handles.activate_batch,'Value')
    set(handles.BashPanel,'Visible','On'); 
    set(handles.compute_lagthreads,'Visible','On');
    set(handles.text24,'Visible','On');
    set(handles.grp1_name,'Visible','On');
    set(handles.batchoptions,'Visible','On');
    set(handles.RUN,'Enable','Off');
    set(handles.parallelize,'Visible','On');

    gsr_Callback(hObject, eventdata, handles);

    handles=guidata(hObject);
    handles.batch=1;
    guidata(hObject,handles);

else
    set(handles.BashPanel,'Visible','Off');
    set(handles.compute_lagthreads,'Visible','Off');
    set(handles.text24,'Visible','Off');
    set(handles.RunBatch,'Visible','Off');
    set(handles.PerformGSR,'Visible','Off');
    set(handles.grp1_name,'Visible','Off');
    set(handles.batchoptions,'Visible','Off');
    set(handles.parallelize,'Visible','Off');
    set(handles.BrainmasksPanel,'Visible','Off');

    handles=guidata(hObject);
    handles.batch=2;
    guidata(hObject,handles);

    set(handles.loadmask,'Enable','On');
    set(handles.load_second_mask,'Enable','On');
    set(handles.chooselag,'Enable','On');
    set(handles.r_min,'Enable','On');
    set(handles.m_min,'Enable','On');
    set(handles.RUN,'Enable','On');
    
end

% --- Executes on button press in perform_pca.
function perform_pca_Callback(hObject, eventdata, handles)
% hObject    handle to perform_pca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of perform_pca


% --- Executes on button press in group_anal.
function group_anal_Callback(hObject, eventdata, handles)
% hObject    handle to group_anal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of group_anal
if get(hObject,'Value');
    set(handles.chooselag,'Enable','Off');
    set(handles.r_min,'Enable','Off');
    set(handles.m_min,'Enable','Off');
    set(handles.loadmask,'Enable','Off');
    set(handles.load_second_mask,'Enable','Off');
    set(handles.parallelize,'Enable','Off');
    set(handles.TR,'Enable','Off');
    set(handles.compute_lagthreads,'Visible','On');
end

% --- Executes on button press in create_batch.
function create_batch_Callback(hObject, eventdata, handles)
% hObject    handle to create_batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'strg_dir') || (isfield(handles,'strg_dir') && ~exist(['' handles.strg_dir '']))
    msgbox('Choose Storage Directory.','Error','error');
else
    %Create a tag
    [status,tag]=system('tag=$(date | md5sum); tag=$(echo "${tag: : -3}"); tag=$(echo "${tag: 0: 10}"); echo $tag');
    spaces=find(isspace(tag)==1);
    
    for spc=1:length(spaces)
        if spaces(spc)==1
            tag=tag(2:end);
        elseif spaces(spc)==length(tag)
            tag=tag(1:end-1);
        else
            tag=[tag(1:spaces(spc)-1) tag(spaces(spc)+1:end)];
        end
    end
    
    %Only Group Analysis Batch 
    if get(handles.group_anal,'Value') %Create Group-Level Analysis Batch script

            name=strcat(get(handles.grp1_name,'string'),'group_anal.mat');

            %%%CREATE BASH-SCRIPT (=BATCH FILE)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if get(handles.parallelize,'Value')
                Servers=handles.Servers;
                for server=1:length(Servers)
                   file=strcat(pwd,filesep,['batch_group_' tag '.sh']);
                end
                
            else
                file=strcat(pwd,filesep,['batch_group_' tag '.sh']);
                fiD=fopen(file,'w');
                fprintf(fiD,'#!/bin/bash \n');
                fprintf(fiD,'\n');
                fprintf(fiD,'SECONDS=0;\n\n');
                fprintf(fiD,'files="{');
                BOLDfilename=handles.BOLDfilename;
                for i=1:length(handles.BOLDfilename)
                    fprintf(fiD,['''' BOLDfilename{i} ''' ']);
                    if i~=length(handles.BOLDfilename)
                        fprintf(fiD,' ');
                    end
                end
                fprintf(fiD,'}";\n');
                fprintf(fiD,'\n');
                fprintf(fiD,['path=' handles.BOLDfilepath ';\n']);
                fprintf(fiD,['name=' name ';\n']);
                fprintf(fiD,['LogDir=' handles.strg_dir ';\n']);
                fprintf(fiD,['storage_dir=' handles.strg_dir ';\n']);
                fprintf(fiD,'newfile=$storage_dir$name;');
                fprintf(fiD,'\necho "COMPUTING ON HOSTNAME $HOSTNAME" > ${LogDir}"$name".txt\n\n');
                fprintf(fiD,'echo "COMPUTING ON HOST $HOST" >> ${LogDir}"$name".txt\n');
                fprintf(fiD,'\necho $(date) >> ${LogDir}"$name".txt\n');
                fprintf(fiD,'echo "Initiating group-level TD computation." >> ${LogDir}$name.txt;\n');
                ThisFile=mfilename('fullpath');
                ThisFolder=ThisFile(1:end-length(mfilename));
                fprintf(fiD,['cd ' ThisFolder '']);
                clear ThisFile ThisFolder
                fprintf(fiD,'\n');
                if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                    fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "groupleveltd($files,''$path'',''$name'',''$storage_dir'');" >> ${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_terminated.txt) &\n']);
                else
                    fprintf(fiD,['MATLABexec=']);
                    fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "groupleveltd($files,''$path'',''$name'',''$storage_dir'');" >> ${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_terminated.txt) &\n']);
                    commandwindow
                    disp('--------------------------------------------------');
                    fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                    disp('--------------------------------------------------');
                end
                fprintf(fiD,'\ncheck=${LogDir}"$name"_terminated.txt;\n');
                fprintf(fiD,'exists=0;');
                fprintf(fiD,'\nwhile [ "$exists" -ne "1" ]; do\n');
                fprintf(fiD,'if [ -e "$check" ]; then\n');
                fprintf(fiD,'exists=1;');
                fprintf(fiD,'\nfi\ndone\n');
                fprintf(fiD,'#Remove termination file as it only exists to indicate that group-level TD computation is terminated');
                fprintf(fiD,'\nrm ${LogDir}"$name"_terminated.txt;\n');
                fprintf(fiD,'#Perform PCA');
                fprintf(fiD,'\n');
                fprintf(fiD,'printf "\\n" >> ${LogDir}"$name".txt; printf "\\n" >> ${LogDir}"$name".txt; echo "Now performing PCA on groupTD" >> ${LogDir}"$name".txt;printf "\\n" >> ${LogDir}"$name".txt;');
                fprintf(fiD,'\necho "$SECONDS seconds elapsed" >> ${LogDir}$name.txt;\n\n');
                if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "performpca(''$newfile'',1,1);" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_PCA_terminated.txt) &']);
                else
                    fprintf(fiD,['MATLABexec=']);
                    fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "performpca(''$newfile'',1,1);" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_PCA_terminated.txt) &']);
                    commandwindow
                    disp('--------------------------------------------------');
                    fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                    disp('--------------------------------------------------');
                end
                fprintf(fiD,'\ncheck=${LogDir}"$name"_PCA_terminated.txt;');
                fprintf(fiD,'\nexists=0;\n');
                fprintf(fiD,'while [ "$exists" -ne "1" ]; do');
                fprintf(fiD,'\nif [ -e "$check" ]; then\n');
                fprintf(fiD,'exists=1;');
                fprintf(fiD,'\nfi\ndone\n');
                fprintf(fiD,'rm ${LogDir}"$name"_PCA_terminated.txt;');
                fprintf(fiD,'\n#Compute 1D lag projection\n');
                fprintf(fiD,'printf "\\n" >> ${LogDir}"$name".txt;printf "\\n" >> ${LogDir}"$name".txt;echo "Now computing group-level 1D lag projection" >> ${LogDir}"$name".txt;printf "\\n" >> ${LogDir}"$name".txt;');
                fprintf(fiD,'\n\necho "$SECONDS seconds elapsed" >> ${LogDir}$name.txt;\n');
                if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                    fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "grouplagproj(''$name'',''$storage_dir'');" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_1d_lagproj_terminated.txt) &']);
                else
                    fprintf(fiD,['MATLABexec=']);
                    fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "grouplagproj(''$name'',''$storage_dir'');" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_1d_lagproj_terminated.txt) &']);
                    commandwindow
                    disp('--------------------------------------------------');
                    fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                    disp('--------------------------------------------------');
                end
                fprintf(fiD,'\ncheck=${LogDir}"$name"_1d_lagproj_terminated.txt;\n');
                fprintf(fiD,'exists=0;');
                fprintf(fiD,'\nwhile [ "$exists" -ne "1" ]; do\n');
                fprintf(fiD,'if [ -e "$check" ]; then');
                fprintf(fiD,'\n');
                fprintf(fiD,'exists=1;');
                fprintf(fiD,'\nfi\ndone\n');
                fprintf(fiD,'printf "\\n" >> ${LogDir}"$name".txt;\n');
                fprintf(fiD,'echo "Batch Group - Level Analysis Terminated!" >> ${LogDir}"$name".txt;\n');
                fprintf(fiD,'rm ${LogDir}"$name"_1d_lagproj_terminated.txt;\n');
                fprintf(fiD,'echo "$SECONDS seconds elapsed" >> ${LogDir}$name.txt;\n');
                fclose(fiD);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                msgbox(['Batch-file created (batch_group_' tag '.sh) '])
                disp(['Batch-file created (batch_group_' tag '.sh) '])
            end

    %.nii-file analysis batch
    elseif get(handles.tdgen,'Value')
        
        sleeping=str2double(get(handles.SleepTime,'string'));
        if isnan(sleeping)
           error('Sleep-Time value not acceptable.'); 
        end

        contents=cellstr(get(handles.m_min,'String'));
        sign_level=contents{get(handles.m_min,'Value')};

        if get(handles.boldfile_descr,'value')==0 || get(handles.text4,'value')==0 
            msgbox('Choose Boldfile/Mask.','Error','error'); 
        elseif isnan(handles.lags) || isnan(handles.rmin)
            msgbox('Parameters not valid.','Error','error');    
        elseif strcmp(sign_level,'Lagged Corr. p')
            msgbox('Parameters not valid.','Error','error'); 
     	elseif isempty(get(handles.TR,'string')) || isnan(str2double(get(handles.TR,'string')))
            msgbox('Enter valid TR!','Error','error');
        else
            % handles.r_min is the text field, handles.rmin is the value
            % handles.m_min is the drop-down menu, handles.mmin is the value

                use_function='local_extrema_xcorr';

            %%%CREATE BASH-SCRIPT (=BATCH FILE)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                file=strcat(handles.strg_dir,['Batchfile_' tag '.sh']);
                fiD=fopen(file,'w');
                fprintf(fiD,'#!/bin/bash \n\n');
                fprintf(fiD,'#logfiles/$file_"$tag" is the logfile of the TD-analysis for each nifti-file.\n');
                
                if get(handles.use_second_mask,'Value')
                    fprintf(fiD,'\n');
                    fprintf(fiD,['mask1=' handles.mask ';\n']);
                    fprintf(fiD,['mask2=' handles.secondmask ';\n']);
                else
                    fprintf(fiD,'\n');
                    fprintf(fiD,['mask1=' handles.mask ';\n']);
                    fprintf(fiD,['mask2=' handles.mask ';\n']);
                end
                
                fprintf(fiD,['path=' handles.BOLDfilepath ';\n']);
                fprintf(fiD,'\n');
                fprintf(fiD,'files=(');
                % If BOLDfilename is no cell, only one file was chosen for batch
                % processing. Make into a cell array of length 1 to be compatible
                % with rest of code.
                BOLDfilename=handles.BOLDfilename;
                
                if ~iscell(BOLDfilename)
                    filename=BOLDfilename;
                    BOLDfilename=cell(1);
                    BOLDfilename{1}=filename;
                end
                
                for i=1:length(BOLDfilename)
                    fprintf(fiD,['''' BOLDfilename{i} ''' ']);
                end
                    fprintf(fiD,');\n\n');
                    fprintf(fiD,['LogDir=' handles.strg_dir ';\n']);
                    fprintf(fiD,['storage_dir=' handles.strg_dir ';\n']);
                    fprintf(fiD,'\n\n');
                    fprintf(fiD,'tag=$(date | md5sum);\n');
                    fprintf(fiD,'tag=$(echo "${tag: : -3}");');
                    fprintf(fiD,'tag=$(echo "${tag: 0: 10}");\n');
                    fprintf(fiD,'number_of_files=$(echo ${#files[@]});\n');
                    fprintf(fiD,'\n');
                    fprintf(fiD,'for file in "${files[@]}"; do\n');
                    fprintf(fiD,'echo "Calculating file $path$file"\n');
                    if get(handles.use_second_mask,'Value')
                        fprintf(fiD,'echo "Brainmask 1 is $mask1"\n');
                        fprintf(fiD,'echo "Brainmask 2 is $mask2"\n');
                    else
                            fprintf(fiD,'echo "Brainmask 1 is $mask1"\n');
                            fprintf(fiD,'echo "Brainmask 2 is $mask2"\n');
                    end
                    fprintf(fiD,'echo "COMPUTING ON HOSTNAME $HOSTNAME" > ${LogDir}"$file"_"$tag".txt\n');
                    fprintf(fiD,'echo "COMPUTING ON HOST $HOST" >> ${LogDir}"$file"_"$tag".txt\n');
                    fprintf(fiD,'echo "storage_dir=$storage_dir" >> ${LogDir}"$file"_"$tag".txt\n');
                    if get(handles.use_second_mask,'Value')
                        fprintf(fiD,['echo "USING MASKS $mask1 and $mask2 ON $file" >> ${LogDir}"$file"_"$tag".txt']);
                    else  
                        fprintf(fiD,['echo "USING MASK $mask ON $file" >> ${LogDir}"$file"_"$tag".txt']);
                    end
                    fprintf(fiD,'\n\n\n');
                    ThisFile=mfilename('fullpath');
                    ThisFolder=ThisFile(1:end-length(mfilename));
                    fprintf(fiD,['cd ' ThisFolder ';\n']);
                    clear ThisFile ThisFolder
                    minmaxdiscr=2;
                    if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                        fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "local_extrema_xcorr(''$path$file'',''$mask1'',''$mask2'',' num2str(handles.lags) ',' num2str(handles.rmin) ',' num2str(handles.mmin) ',' num2str(str2double(get(handles.TR,'string'))) ',' num2str(minmaxdiscr) ',''$storage_dir'',''$tag'',1);" >> ${LogDir}"$file"_"$tag".txt 2>&1; touch ${LogDir}"$file"_"$tag"_terminated.txt) &\n']);
                    else
                        fprintf(fiD,['MATLABexec=']);
                        fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "local_extrema_xcorr(''$path$file'',''$mask1'',''$mask2'',' num2str(handles.lags) ',' num2str(handles.rmin) ',' num2str(handles.mmin) ',' num2str(str2double(get(handles.TR,'string'))) ',' num2str(minmaxdiscr) ',''$storage_dir'',''$tag'',1);" >> ${LogDir}"$file"_"$tag".txt 2>&1; touch ${LogDir}"$file"_"$tag"_terminated.txt) &\n']);
                        commandwindow
                        disp('--------------------------------------------------');
                        fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                        disp('--------------------------------------------------');
                    end
                    fprintf(fiD,['sleep ' num2str(sleeping) ';\n']);
                    fprintf(fiD,'done\n\n');
                    fprintf(fiD,'files_terminated=0;\n');
                    fprintf(fiD,'#Remove backup files to be sure.\n');
                    fprintf(fiD,'rm ${LogDir}*~;\n');
                    fprintf(fiD,'while [ "$files_terminated" -ne "$number_of_files" ]; do\n');
                    fprintf(fiD,['files_terminated=$(ls ${LogDir} | grep $tag | grep terminated |wc -l);']);
                    fprintf(fiD,'\n');
                    fprintf(fiD,'done\n');
                    fprintf(fiD,['rm ${LogDir}*_"$tag"_terminated.txt;']);
                    fprintf(fiD,'\n\necho "TD-computation terminated."');
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    msgbox(['Batch-file created (' file ') ']);
                    disp(['Batch-file created (' file ') ']);
                    set(handles.create_batch,'Enable','On');
        end

    %TD generation + ensuing group analysis
    elseif get(handles.tdgen_group,'Value')

        sleeping=str2double(get(handles.SleepTime,'string'));
        if isnan(sleeping)
           error('Sleep-Time value not acceptable.'); 
        end
        
        name=get(handles.grp1_name,'string');

        contents=cellstr(get(handles.m_min,'String'));
        sign_level=contents{get(handles.m_min,'Value')};

        if get(handles.boldfile_descr,'value')==0 || get(handles.text4,'value')==0 
            msgbox('Choose Boldfile/Mask.','Error','error'); 
        elseif isnan(handles.lags) || isnan(handles.rmin)
            msgbox('Parameters not valid.','Error','error'); 
        elseif strcmp(sign_level,'Lagged Corr. p')
            msgbox('Parameters not valid.','Error','error'); 
        elseif isempty(get(handles.grp1_name,'string'))
            msgbox('Enter group-level filename!','Error','error');
     	elseif isempty(get(handles.TR,'string')) || isnan(str2double(get(handles.TR,'string')))
            msgbox('Enter valid TR!','Error','error');
        else
            % handles.r_min is the text field, handles.rmin is the value
            % handles.m_min is the drop-down menu, handles.mmin is the value

                use_function='local_extrema_xcorr';

            %%%CREATE BASH-SCRIPT (=BATCH FILE)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if get(handles.parallelize,'Value')
                Servers=handles.servers;
                answer=handles.parallel;
                
                total_no_files=0;
                for ans=1:length(answer)
                    total_no_files=total_no_files+str2double(answer{ans});
                end
                
                file=strcat(handles.strg_dir,['GroupBatchFile_' tag '.sh']);
                fiD=fopen(file,'w');
                fprintf(fiD,'#!/bin/bash \n\n');
                fprintf(fiD,['LogDir=' handles.strg_dir '\n']);
                fprintf(fiD,['tag=' tag '']);
                fprintf(fiD,'\n');

                for server=1:length(Servers) 
                    fprintf(fiD,['$(ssh ' handles.ServerUser{1} '@' Servers{server} ' ''bash -s'' < ./batch_TD_group_' Servers{server} '_' tag '.sh $tag>${LogDir}batch_TD_group_' Servers{server} '_"$tag".txt 2>&1 && exit) & \n']);
                end
                
                fprintf(fiD,'\n');
                fclose(fiD);

                for server=1:length(Servers)
                    
                    no_files=0;
                    for fl=1:(server-1)
                       % no_files is the number of preceding files in the files-list processed by other
                       % servers.
                        no_files=no_files+str2double(answer{fl});
                    end

                    file=strcat(handles.strg_dir,filesep,['batch_TD_group_' Servers{server} '_' tag '.sh']);
                    fiD=fopen(file,'w');
                    fprintf(fiD,'#!/bin/bash \n\n');
                    fprintf(fiD,['total_no_files=' num2str(total_no_files) '']);
                    fprintf(fiD,'\n\n#logfiles/$file_"$tag" is the logfile of the TD-analysis for each nifti-file.\n');
                    if get(handles.use_second_mask,'Value')
                        fprintf(fiD,'\n');
                        fprintf(fiD,['mask1=' handles.mask ';\n']);
                        fprintf(fiD,['mask2=' handles.secondmask ';\n']);
                    else
                        fprintf(fiD,'\n');
                        fprintf(fiD,['mask1=' handles.mask ';\n']);
                        fprintf(fiD,['mask2=' handles.mask ';\n']);
                    end
                    fprintf(fiD,['path=' handles.BOLDfilepath ';\n']);
                    fprintf(fiD,'\ntag=$1\n');
                    fprintf(fiD,'files=(');
                    % If BOLDfilename is no cell, only one file was chosen for batch
                    % processing. Make into a #!cell array of length 1 to be compatible
                    % with rest of code.
                    BOLDfilename=handles.BOLDfilename;
                    if ~iscell(BOLDfilename)
                        filename=BOLDfilename;
                        BOLDfilename=cell(1);
                        BOLDfilename{1}=filename;
                    end
                    for i=(no_files+1):(no_files+str2double(answer{server}))
                        fprintf(fiD,['''' BOLDfilename{i} ''' ']);
                    end
                        fprintf(fiD,');\n');
                        fprintf(fiD,['LogDir=' handles.strg_dir ';\n']);
                        fprintf(fiD,['storage_dir=' handles.strg_dir ';\n']);
                        fprintf(fiD,'number_of_files=$(echo ${#files[@]});\n');
                        fprintf(fiD,'for file in "${files[@]}"; do\n');
                        fprintf(fiD,'    echo "Calculating file $path$file"\n');
                        if get(handles.use_second_mask,'Value')
                            fprintf(fiD,'    echo "Brainmask 1 is $mask1"\n');
                            fprintf(fiD,'    echo "Brainmask 2 is $mask2"\n');
                        else
                            fprintf(fiD,'    echo "Brainmask 1 is $mask1"\n');
                            fprintf(fiD,'    echo "Brainmask 2 is $mask2"\n');
                        end
                        fprintf(fiD,'    echo "COMPUTING ON HOSTNAME $HOSTNAME" > ${LogDir}"$file"_"$tag".txt\n');
                        fprintf(fiD,'    echo "COMPUTING ON HOST $HOST" >> ${LogDir}"$file"_"$tag".txt\n');
                        if get(handles.use_second_mask,'Value')
                            fprintf(fiD,['    echo "USING MASKS $mask1 and $mask2 ON $file" >> ${LogDir}"$file"_"$tag".txt']);
                        else  
                            fprintf(fiD,['    echo "USING MASK $mask1 ON $file" >> ${LogDir}"$file"_"$tag".txt']);
                        end
                        fprintf(fiD,'    echo "storage_dir=$storage_dir" >> ${LogDir}"$file"_"$tag".txt\n\n');
                        ThisFile=mfilename('fullpath');
                        ThisFolder=ThisFile(1:end-length(mfilename));
                        fprintf(fiD,['    cd ' ThisFolder '']);
                        clear ThisFile ThisFolder
                        minmaxdiscr=2;
                        fprintf(fiD,'\n');
                        str2double(get(handles.TR,'string'))
                        if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                            fprintf(fiD,['    $(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "local_extrema_xcorr(''$path$file'',''$mask1'',''$mask2'',' num2str(handles.lags) ',' num2str(handles.rmin) ',' num2str(handles.mmin) ',' num2str(str2double(get(handles.TR,'string'))) ',' num2str(minmaxdiscr) ',''$storage_dir'',''$tag'',1);" >> ${LogDir}"$file"_"$tag".txt 2>&1; touch ${LogDir}"$file"_' Servers{server} '_"$tag"_terminated.txt) &\n']);
                        else
                            fprintf(fiD,['    MATLABexec=']);
                            fprintf(fiD,['    $(${MATLABexec} -nodesktop -nosplash -r "local_extrema_xcorr(''$path$file'',''$mask1'',''$mask2'',' num2str(handles.lags) ',' num2str(handles.rmin) ',' num2str(handles.mmin) ',' num2str(str2double(get(handles.TR,'string'))) ',' num2str(minmaxdiscr) ',''$storage_dir'',''$tag'',1);" >> ${LogDir}"$file"_"$tag".txt 2>&1; touch ${LogDir}"$file"_' Servers{server} '_"$tag"_terminated.txt) &\n']);
                            commandwindow
                            disp('--------------------------------------------------');
                            fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                            disp('--------------------------------------------------');
                        end
                        fprintf(fiD,['sleep ' num2str(sleeping) ';\n']);
                        fprintf(fiD,'done\n\n');
                        
                        % First server is the one for group analysis
                        if server==1
                            fprintf(fiD,'files_terminated=0;\n');
                            fprintf(fiD,'#Remove backup files to be sure.\n');
                            fprintf(fiD,'rm ${LogDir}*~;\n');
                            fprintf(fiD,'while [ "$files_terminated" -ne "$total_no_files" ]; do\n');
                            fprintf(fiD,['files_terminated=$(ls ${LogDir} | grep $tag | grep terminated |wc -l);']);
                            fprintf(fiD,'\n');
                            fprintf(fiD,'done\n');
                            fprintf(fiD,'\nrm ${LogDir}*_"$tag"_terminated.txt;\n');
                            fprintf(fiD,'#TD-computation terminated. Initiate group analysis. \n');
                            fprintf(fiD,'#logfiles/$name is the logfile of the group analysis. \n \n \n');
                            fprintf(fiD,['path=$storage_dir;\n']);
                            fprintf(fiD,['#name is the name of the group-file\n']);
                            fprintf(fiD,['name=' name '_"$tag".mat;\n']);
                            fprintf(fiD,'newfile=$path$name;\n');
                            fprintf(fiD,'echo "COMPUTING ON HOSTNAME $HOSTNAME" > ${LogDir}"$name".txt\n');
                            fprintf(fiD,'echo "COMPUTING ON HOST $HOST" >> ${LogDir}"$name".txt\n');
                            fprintf(fiD,'\necho $(date) >> ${LogDir}"$name".txt\n');
                            fprintf(fiD,'echo "Initiating group-level TD computation." >> ${LogDir}"$name".txt;\n');
                            ThisFile=mfilename('fullpath');
                            ThisFolder=ThisFile(1:end-length(mfilename));
                            fprintf(fiD,['cd ' ThisFolder ';\n']);
                            clear ThisFile ThisFolder
                            if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                                fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "groupleveltd($files,''$path'',''$name'',''$storage_dir'');" >> ${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_terminated.txt) &\n']);
                            else
                                fprintf(fiD,['MATLABexec=']);
                                fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "groupleveltd($files,''$path'',''$name'',''$storage_dir'');" >> ${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_terminated.txt) &\n']);
                                commandwindow
                                disp('--------------------------------------------------');
                                fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                                disp('--------------------------------------------------');
                            end
                            fprintf(fiD,'\ncheck=${LogDir}"$name"_terminated.txt;\n');
                            fprintf(fiD,'exists=0;\n');
                            fprintf(fiD,'while [ "$exists" -ne "1" ]; do\n');
                            fprintf(fiD,'    if [ -e "$check" ]; then\n');
                            fprintf(fiD,'        exists=1;\n');
                            fprintf(fiD,'    fi\ndone\n');
                            fprintf(fiD,'#Remove termination file as it only exists to indicate that group-level TD computation is terminated\n');
                            fprintf(fiD,'rm ${LogDir}"$name"_terminated.txt;\n');
                            fprintf(fiD,'#Perform PCA\n');
                            fprintf(fiD,'printf "\\n" >> ${LogDir}"$name".txt; printf "\\n" >> ${LogDir}"$name".txt; echo "Now performing PCA on groupTD" >> ${LogDir}"$name".txt;printf "\\n" >> ${LogDir}"$name".txt;\n');
                            if get(handles.compute_lagthreads,'Value')
                                if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                                    fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "performpca(''$newfile'',1,1,1);" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_PCA_terminated.txt) &']);
                                else
                                    fprintf(fiD,['MATLABexec=']);
                                    fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "performpca(''$newfile'',1,1,1);" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_PCA_terminated.txt) &']);
                                    commandwindow
                                    disp('--------------------------------------------------');
                                    fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                                    disp('--------------------------------------------------');
                                end
                            else
                                if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                                    fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "performpca(''$newfile'',1,1,0);" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_PCA_terminated.txt) &']);
                                else
                                    fprintf(fiD,['MATLABexec=']);
                                    fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "performpca(''$newfile'',1,1,0);" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_PCA_terminated.txt) &']);
                                    commandwindow
                                    disp('--------------------------------------------------');
                                    fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                                    disp('--------------------------------------------------');
                                end
                            end
                            fprintf(fiD,'\ncheck=${LogDir}"$name"_PCA_terminated.txt;\n');
                            fprintf(fiD,'exists=0;\n');
                            fprintf(fiD,'while [ "$exists" -ne "1" ]; do\n');
                            fprintf(fiD,'if [ -e "$check" ]; then\n');
                            fprintf(fiD,'exists=1;');
                            fprintf(fiD,'\nfi\ndone\n');
                            fprintf(fiD,'rm ${LogDir}"$name"_PCA_terminated.txt;\n');
                            fprintf(fiD,'#Compute 1D lag projection\n');
                            fprintf(fiD,'printf "\\n" >> ${LogDir}"$name".txt;printf "\\n" >> ${LogDir}"$name".txt;echo "Now computing group-level 1D lag projection" >> ${LogDir}"$name".txt;printf "\\n" >> ${LogDir}"$name".txt;\n\n');
                            if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                                fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "grouplagproj(''$name'',''$path'');" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_1d_lagproj_terminated.txt) &']);
                            else
                                fprintf(fiD,['MATLABexec=']);
                                fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "grouplagproj(''$name'',''$path'');" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_1d_lagproj_terminated.txt) &']);
                                commandwindow
                                disp('--------------------------------------------------');
                                fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                                disp('--------------------------------------------------');
                            end
                            fprintf(fiD,'\ncheck=${LogDir}"$name"_1d_lagproj_terminated.txt;\n');
                            fprintf(fiD,'exists=0;\n');
                            fprintf(fiD,'while [ "$exists" -ne "1" ]; do\n');
                            fprintf(fiD,'if [ -e "$check" ]; then\n');
                            fprintf(fiD,'exists=1;');
                            fprintf(fiD,'\nfi\ndone\n');
                            fprintf(fiD,'printf "\\n" >> ${LogDir}"$name".txt;\n\n\n');
                            fprintf(fiD,'echo "Batch Group - Level Analysis Terminated!" >> ${LogDir}"$name".txt;\n');
                            fprintf(fiD,'rm ${LogDir}"$name"_1d_lagproj_terminated.txt;\n\n');
                        else % Not the server for group analysis
                            
                            
                        end
                        
                        fclose(fiD);
                        
                end
                
                msgbox(['Batch-files created (' handles.strg_dir 'GroupBatchFile_' tag '.sh) ']);
                disp(['Batch-files created (' handles.strg_dir 'GroupBatchFile_' tag '.sh) ']);
                
            else %No parallelization
                file=strcat(handles.strg_dir,['batch_TD_group_' tag '.sh']);
                fiD=fopen(file,'w');
                fprintf(fiD,'#!/bin/bash \n\n');
                fprintf(fiD,'#logfiles/$file_"$tag" is the logfile of the TD-analysis for each nifti-file.\n');
                
                if get(handles.use_second_mask,'Value')
                    fprintf(fiD,'\n');
                    fprintf(fiD,['mask1=' handles.mask ';\n']);
                    fprintf(fiD,['mask2=' handles.secondmask ';\n']);
                else
                    fprintf(fiD,'\n');
                    fprintf(fiD,['mask1=' handles.mask ';\n']);
                    fprintf(fiD,['mask2=' handles.mask ';\n']);
                end
                
                fprintf(fiD,['path=' handles.BOLDfilepath ';\n']);
                fprintf(fiD,'\n');
                fprintf(fiD,'files=(');
                % If BOLDfilename is no cell, only one file was chosen for batch
                % processing. Make into a cell array of length 1 to be compatible
                % with rest of code.
                BOLDfilename=handles.BOLDfilename;
                
                if ~iscell(BOLDfilename)
                    filename=BOLDfilename;
                    BOLDfilename=cell(1);
                    BOLDfilename{1}=filename;
                end
                
                for i=1:length(BOLDfilename)
                    fprintf(fiD,['''' BOLDfilename{i} ''' ']);
                end
                    fprintf(fiD,');\n\n');
                    fprintf(fiD,['LogDir=' handles.strg_dir ';\n']);
                    fprintf(fiD,['storage_dir=' handles.strg_dir ';\n']);
                    fprintf(fiD,'\n\n');
                    fprintf(fiD,'tag=$(date | md5sum);\n');
                    fprintf(fiD,'tag=$(echo "${tag: : -3}");');
                    fprintf(fiD,'tag=$(echo "${tag: 0: 10}");\n');
                    %fprintf(fiD,'batchfile_id=$(date | md5sum);\n');
                    %fprintf(fiD,'batchfile_id=$(echo "${batchfile_id: : -3}");\n');
                    %fprintf(fiD,'touch ${LogDir}."$batchfile_id".txt;\n');
                    %fprintf(fiD,'batchlog=${LogDir}."$batchfile_id".txt;\n');
                    fprintf(fiD,'number_of_files=$(echo ${#files[@]});\n');
                    fprintf(fiD,'\n');
                    fprintf(fiD,'for file in "${files[@]}"; do\n');
                    fprintf(fiD,'echo "Calculating file $path$file"\n');
                    if get(handles.use_second_mask,'Value')
                        fprintf(fiD,'echo "Brainmask 1 is $mask1"\n');
                        fprintf(fiD,'echo "Brainmask 2 is $mask2"\n');
                    else
                            fprintf(fiD,'echo "Brainmask 1 is $mask1"\n');
                            fprintf(fiD,'echo "Brainmask 2 is $mask2"\n');
                    end
                    fprintf(fiD,'echo "COMPUTING ON HOSTNAME $HOSTNAME" > ${LogDir}"$file"_"$tag".txt\n');
                    fprintf(fiD,'echo "COMPUTING ON HOST $HOST" >> ${LogDir}"$file"_"$tag".txt\n');
                    fprintf(fiD,'echo "storage_dir=$storage_dir" >> ${LogDir}"$file"_"$tag".txt\n');
                    if get(handles.use_second_mask,'Value')
                        fprintf(fiD,['echo "USING MASKS $mask1 and $mask2 ON $file" >> ${LogDir}"$file"_"$tag".txt']);
                    else  
                        fprintf(fiD,['echo "USING MASK $mask1 ON $file" >> ${LogDir}"$file"_"$tag".txt']);
                    end
                    fprintf(fiD,'\n\n\n');
                    ThisFile=mfilename('fullpath');
                    ThisFolder=ThisFile(1:end-length(mfilename));
                    fprintf(fiD,['cd ' ThisFolder ';\n']);
                    clear ThisFile ThisFolder
                    minmaxdiscr=2;
                    if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                        fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "local_extrema_xcorr(''$path$file'',''$mask1'',''$mask2'',' num2str(handles.lags) ',' num2str(handles.rmin) ',' num2str(handles.mmin) ',' num2str(str2double(get(handles.TR,'string'))) ',' num2str(minmaxdiscr) ',''$storage_dir'',''$tag'',1);" >> ${LogDir}"$file"_"$tag".txt 2>&1; touch ${LogDir}"$file"_"$tag"_terminated.txt) &\n']);
                    else
                        fprintf(fiD,['MATLABexec=']);
                        fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "local_extrema_xcorr(''$path$file'',''$mask1'',''$mask2'',' num2str(handles.lags) ',' num2str(handles.rmin) ',' num2str(handles.mmin) ',' num2str(str2double(get(handles.TR,'string'))) ',' num2str(minmaxdiscr) ',''$storage_dir'',''$tag'',1);" >> ${LogDir}"$file"_"$tag".txt 2>&1; touch ${LogDir}"$file"_"$tag"_terminated.txt) &\n']);
                        commandwindow
                        disp('--------------------------------------------------');
                        fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                        disp('--------------------------------------------------');
                    end
                    fprintf(fiD,['sleep ' num2str(sleeping) ';\n']);
                    fprintf(fiD,'done\n\n');
                    fprintf(fiD,'files_terminated=0;\n');
                    fprintf(fiD,'#Remove backup files to be sure.\n');
                    fprintf(fiD,'rm ${LogDir}*~;\n');
                    fprintf(fiD,'while [ "$files_terminated" -ne "$number_of_files" ]; do\n');
                    fprintf(fiD,['files_terminated=$(ls ${LogDir} | grep $tag | grep terminated |wc -l);']);
                    fprintf(fiD,'\n');
                    fprintf(fiD,'done\n');
                    fprintf(fiD,['rm ${LogDir}*_"$tag"_terminated.txt;']);
                    fprintf(fiD,'\n\n\n\n');
                    fprintf(fiD,'#TD-computation terminated. Initiate group analysis.\n\n');
                    fprintf(fiD,'#logfiles/$name is the logfile of the group analysis.\n\n');
                    fprintf(fiD,'SECONDS=0;\n');
                    fprintf(fiD,'files="{');
                    BOLDfilename=handles.BOLDfilename;
                    for i=1:length(handles.BOLDfilename)
                        filename=BOLDfilename{i};
                        fprintf(fiD,[' ''' filename(1:end-4) '_"$tag".mat'' ']);
                        if i~=length(handles.BOLDfilename)
                            fprintf(fiD,',');
                        end
                    end
                    fprintf(fiD,'}";\n');
                    fprintf(fiD,['path=$storage_dir;\n']);
                    fprintf(fiD,['name=' name '_group_anal.mat;\n']);
                    fprintf(fiD,'newfile=$path$name;\n');
                    fprintf(fiD,'echo "COMPUTING ON HOSTNAME $HOSTNAME" > ${LogDir}"$name".txt\n');
                    fprintf(fiD,'echo "COMPUTING ON HOST $HOST" >> ${LogDir}"$name".txt\n');
                    fprintf(fiD,'echo "storage_dir=$storage_dir" >> ${LogDir}"$name".txt\n\n\n');
                    fprintf(fiD,'echo "Initiating group-level TD computation." >> ${LogDir}$name.txt;\n');
                    ThisFile=mfilename('fullpath');
                    ThisFolder=ThisFile(1:end-length(mfilename));
                    fprintf(fiD,['cd ' ThisFolder ';\n']);
                    clear ThisFile ThisFolder
                    fprintf(fiD,'\n\n\n\n');
                    if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                        fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "groupleveltd($files,''$path'',''$name'',''$storage_dir'');" >> ${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_terminated.txt) &\n']);
                    else
                        fprintf(fiD,['MATLABexec=']);
                        fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "groupleveltd($files,''$path'',''$name'',''$storage_dir'');" >> ${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_terminated.txt) &\n']);
                        commandwindow
                        disp('--------------------------------------------------');
                        fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                        disp('--------------------------------------------------');
                    end
                    fprintf(fiD,'\ncheck=${LogDir}"$name"_terminated.txt;');
                    fprintf(fiD,'\nexists=0;\n');
                    fprintf(fiD,'while [ "$exists" -ne "1" ]; do\n');
                    fprintf(fiD,'if [ -e "$check" ]; then');
                    fprintf(fiD,'\nexists=1;');
                    fprintf(fiD,'\nfi\nfone\n');
                    fprintf(fiD,'#Remove termination file as it only exists to indicate that group-level TD computation is terminated');
                    fprintf(fiD,'\nrm ${LogDir}"$name"_terminated.txt;\n');
                    fprintf(fiD,'#Perform PCA\n');
                    fprintf(fiD,'printf "\\n" >> ${LogDir}"$name".txt; printf "\\n" >> ${LogDir}"$name".txt; echo "Now performing PCA on groupTD" >> ${LogDir}"$name".txt;printf "\\n" >> ${LogDir}"$name".txt;');
                    fprintf(fiD,'\necho "$SECONDS seconds elapsed" >> ${LogDir}$name.txt;\n\n');
                    if get(handles.compute_lagthreads,'Value');
                        if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                            fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "performpca(''$newfile'',1,1,1);" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_PCA_terminated.txt) &']);
                        else
                            fprintf(fiD,['MATLABexec=']);
                            fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "performpca(''$newfile'',1,1,1);" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_PCA_terminated.txt) &']);
                            commandwindow
                            disp('--------------------------------------------------');
                            fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                            disp('--------------------------------------------------');
                        end
                    else
                        if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                            fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "performpca(''$newfile'',1,1,0);" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_PCA_terminated.txt) &']);
                        else
                            fprintf(fiD,['MATLABexec=']);
                            fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "performpca(''$newfile'',1,1,0);" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_PCA_terminated.txt) &']);
                            commandwindow
                            disp('--------------------------------------------------');
                            fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                            disp('--------------------------------------------------');
                        end
                    end
                    fprintf(fiD,'\ncheck=${LogDir}"$name"_PCA_terminated.txt;\n');
                    fprintf(fiD,'exists=0;');
                    fprintf(fiD,'\nwhile [ "$exists" -ne "1" ]; do');
                    fprintf(fiD,'\nif [ -e "$check" ]; then\n');
                    fprintf(fiD,'exists=1;');
                    fprintf(fiD,'\nfi\ndone\n');
                    fprintf(fiD,'rm ${LogDir}"$name"_PCA_terminated.txt;');
                    fprintf(fiD,'\n#Compute 1D lag projection\n');
                    fprintf(fiD,'printf "\\n" >> ${LogDir}"$name".txt;printf "\\n" >> ${LogDir}"$name".txt;echo "Now computing group-level 1D lag projection" >> ${LogDir}"$name".txt;printf "\\n" >> ${LogDir}"$name".txt;');
                    fprintf(fiD,'\n\necho "$SECONDS seconds elapsed" >> ${LogDir}$name.txt;');
                    if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                        fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "grouplagproj(''$name'',''$path'');" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_1d_lagproj_terminated.txt) &']);
                    else
                        fprintf(fiD,['MATLABexec=']);
                        fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "grouplagproj(''$name'',''$path'');" >>${LogDir}"$name".txt 2>&1; touch ${LogDir}"$name"_1d_lagproj_terminated.txt) &']);
                        commandwindow
                        disp('--------------------------------------------------');
                        fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                        disp('--------------------------------------------------');
                    end
                    fprintf(fiD,'\ncheck=${LogDir}"$name"_1d_lagproj_terminated.txt;');
                    fprintf(fiD,'\nexists=0;\n');
                    fprintf(fiD,'while [ "$exists" -ne "1" ]; do');
                    fprintf(fiD,'\nif [ -e "$check" ]; then\n');
                    fprintf(fiD,'exists=1;');
                    fprintf(fiD,'\nfi\ndone\n');
                    fprintf(fiD,'printf "\\n" >> ${LogDir}"$name".txt;');
                    fprintf(fiD,'\necho "Batch Group - Level Analysis Terminated!" >> ${LogDir}"$name".txt;\n');
                    fprintf(fiD,'rm ${LogDir}"$name"_1d_lagproj_terminated.txt;');
                    fprintf(fiD,'\necho "$SECONDS seconds elapsed" >> ${LogDir}$name.txt;\n');
                    fclose(fiD);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    msgbox(['Batch-file created (' file ') ']);
                    disp(['Batch-file created (' file ') ']);
                    set(handles.create_batch,'Enable','On');
            end

        end
            
    elseif get(handles.gsr,'Value')
        if ~isfield(handles,'gsr_brainmask')
            error('Please specify GSR brainmask!');
        end
            %%%CREATE BASH-SCRIPT (=BATCH FILE)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            file=strcat(handles.strg_dir,'gsr.sh');
            fiD=fopen(file,'w');
            fprintf(fiD,'#!/bin/bash \n');
            fprintf(fiD,['LogDir=' handles.strg_dir ';\n']);
            fprintf(fiD,['storagedir=' handles.strg_dir ';\n']);
            fprintf(fiD,'\nfiles="{');
            BOLDfilename=handles.BOLDfilename;
            if iscell(handles.BOLDfilename)
                for i=1:length(handles.BOLDfilename)
                    filename=BOLDfilename{i};
                    fprintf(fiD,['''' strcat(handles.BOLDfilepath,filename) '''']);
                    if i~=length(handles.BOLDfilename)
                        fprintf(fiD,',');
                    end
                end
            else
                fprintf(fiD,['''' strcat(handles.BOLDfilepath,handles.BOLDfilename) '''']);
            end
            fprintf(fiD,'}";\n');
            fprintf(fiD,'\nbrainmasks="{');
            MaskFilenames=handles.gsr_brainmasks_filenames;
            if iscell(handles.gsr_brainmasks_filenames)
                for i=1:length(handles.gsr_brainmasks_filenames)
                    filename=MaskFilenames{i};
                    fprintf(fiD,['''' strcat(handles.gsr_brainmasks_path,filename) '''']);
                    if i~=length(handles.gsr_brainmasks_filenames)
                        fprintf(fiD,',');
                    end
                end
            else
                fprintf(fiD,['''' strcat(handles.gsr_brainmasks_path,handles.gsr_brainmasks_filenames) '''']);
            end
            fprintf(fiD,'}";\n');
            
            ThisFile=mfilename('fullpath');
            ThisPath=ThisFile(1:end-length(mfilename));
            fprintf(fiD,['cd ' ThisPath ';\n']);
            if exist(strcat(matlabroot,filesep,'bin',filesep,'matlab'),'file')==2
                fprintf(fiD,['$(' strcat(matlabroot,filesep) 'bin' filesep 'matlab -nodesktop -nosplash -r "gsr($files,$brainmasks,''$storagedir'');" > ${LogDir}gsrLog.txt 2>&1;  )']);
            else
                fprintf(fiD,['MATLABexec=']);
                fprintf(fiD,['$(${MATLABexec} -nodesktop -nosplash -r "gsr($files,$brainmasks,''$storagedir'');" > ${LogDir}gsrLog.txt 2>&1;  )']);
                commandwindow
                disp('--------------------------------------------------');
                fprintf('\n\nNB: Could not find MATLAB executable. Insert MATLAB exec directory manually in bash-file in line that says ''MATLABexec=''.\n\n');
                disp('--------------------------------------------------');
            end
            fprintf(fiD,'\n');
            fclose(fiD);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            msgbox(['Batch-file created (' strcat(handles.strg_dir,'gsr.sh') ')']);
            disp(['Batch-file created (' strcat(handles.strg_dir,'gsr.sh') ')']);
            
    end
  
end


handles=guidata(hObject);
handles.strg_dir=0;
set(handles.storage_dir,'String','Choose Storage Directory.');
guidata(hObject,handles);


function grp1_name_Callback(hObject, eventdata, handles)
% hObject    handle to grp1_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of grp1_name as text
%        str2double(get(hObject,'String')) returns contents of grp1_name as a double


% --- Executes during object creation, after setting all properties.
function grp1_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to grp1_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in choosedir.
function choosedir_Callback(hObject, eventdata, handles)
% hObject    handle to choosedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'RecentStrgDir')
	strg_dir=uigetdir(handles.RecentStrgDir);
else
	strg_dir=uigetdir(pwd);
end

if strfind(strg_dir,' ');
   msgbox('Make sure there are no space characters in the storage directory path.');
   error('Make sure there are no space characters in the storage directory path.'); 
end

if ~isequal(strg_dir,0)
	handles=guidata(hObject);
	handles.RecentStrgDir=strg_dir;
	guidata(hObject,handles);
else
    return
end


if ~isempty(strg_dir)
    % handles.strg_dir should always end in filesep
    if ~strcmp(strg_dir(end),filesep);
        strg_dir=strcat(strg_dir,filesep);
    end

    set(handles.storage_dir,'String',strg_dir); 
    drawnow;
   
    handles=guidata(hObject);
    handles.strg_dir=strg_dir;
    guidata(hObject,handles);
end


% --- Executes on button press in tdgen.
function tdgen_Callback(hObject, eventdata, handles)
% hObject    handle to tdgen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdgen
if get(hObject,'Value');
    set(handles.chooselag,'Enable','On');
    set(handles.r_min,'Enable','On');
    set(handles.m_min,'Enable','On');
    set(handles.loadmask,'Enable','On');
    set(handles.load_second_mask,'Enable','On');
    set(handles.parallelize,'Visible','Off');
    set(handles.gsr_brainmasks,'Visible','Off');
    set(handles.BrainmasksPanel,'Visible','Off');
    set(handles.text24,'Visible','On');
    set(handles.grp1_name,'Visible','On');
    set(handles.compute_lagthreads,'Visible','On');
    set(handles.TR,'Enable','On');
    set(handles.SleepTime,'Visible','On');
    set(handles.text33,'Visible','On');
    set(handles.compute_lagthreads,'Visible','Off');
    set(handles.use_second_mask,'Enable','On');
    set(handles.PerformGSR,'Visible','Off');
    set(handles.RunBatch,'Visible','On');
    set(handles.text24,'Visible','Off');
    set(handles.grp1_name,'Visible','Off');
end


% --- Executes on button press in tdgen_group.
function tdgen_group_Callback(hObject, eventdata, handles)
% hObject    handle to tdgen_group (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tdgen_group


if get(hObject,'Value');
    set(handles.chooselag,'Enable','On');
    set(handles.r_min,'Enable','On');
    set(handles.m_min,'Enable','On');
    set(handles.parallelize,'Enable','On');
    set(handles.TR,'Enable','On');
    set(handles.compute_lagthreads,'Visible','On');
    set(handles.loadmask,'Enable','On');
    set(handles.load_second_mask,'Enable','On');
    set(handles.text24,'Visible','On');
    set(handles.grp1_name,'Visible','On');
    set(handles.gsr_brainmasks,'Visible','Off');
    set(handles.BrainmasksPanel,'Visible','Off');
    set(handles.parallelize,'Visible','On');
    set(handles.SleepTime,'Visible','On');
    set(handles.text33,'Visible','On');
    set(handles.use_second_mask,'Enable','On');
    set(handles.PerformGSR,'Visible','Off');
    set(handles.RunBatch,'Visible','On');
end


% --- Executes on button press in use_second_mask.
function use_second_mask_Callback(hObject, eventdata, handles)
% hObject    handle to use_second_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_second_mask
if get(handles.use_second_mask,'Value')
   set(handles.second_mask,'Visible','On'); 
   set(handles.load_second_mask,'Visible','On'); 
else
   set(handles.second_mask,'Visible','Off'); 
   set(handles.load_second_mask,'Visible','Off'); 
end

% --- Executes on button press in load_second_mask.
function load_second_mask_Callback(hObject, eventdata, handles)
% hObject    handle to load_second_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.second_mask,'Value',0.0)
ThisFile=mfilename('fullpath');
path=strcat(ThisFile(1:end-length(mfilename)),'masks',filesep);
[maskname,maskpath]=uigetfile(['' path '*.nii']);

if strfind(maskpath,' ');
   msgbox('Make sure there are no space characters in the mask''s directory.');
   error('Make sure there are no space characters in the mask''s directory.'); 
end

clear ThisFile path
mask=strcat(maskpath,maskname);
set(handles.second_mask,'String',['Mask: ' maskname ''])
set(handles.second_mask,'TooltipString',['Path: ' maskpath '']);
if exist(mask)
    set(handles.second_mask,'Value',1.0);
    
    handles=guidata(hObject);
    handles.secondmask=mask;
    guidata(hObject,handles);
    
end

% --- Executes on button press in parallelize.
function parallelize_Callback(hObject, eventdata, handles)
% hObject    handle to parallelize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parallelize

set(handles.create_batch,'enable','off');
if get(handles.parallelize,'Value')
    NoServers=inputdlg('Enter the number of servers','Servers',1);

    ServersPrompt=cell(1,str2double(NoServers));
    [ServersPrompt{1,:}]=deal('Specify Server name:'); 
    Servers=inputdlg(ServersPrompt,'Server names',1);        

    prompt={ strcat('Enter no. of files to be processed by ''',Servers{1},'''') };
    for server=2:length(Servers)
       prompt=[prompt, { strcat('Enter no. of files to be processed by ''',Servers{server},'''') } ];
    end

    ServerUser=inputdlg('Enter the username for the servers (needs to be the same username for each server).','User',1);

    dlg_title='Input';
    num_lines=1;
    answer=inputdlg(prompt,dlg_title,num_lines);
    handles=guidata(hObject);
    handles.servers=Servers;
    handles.parallel=answer;
    handles.ServerUser=ServerUser;
    guidata(hObject,handles);

    if 0
    if get(handles.tdgen_group,'Value')

       str_group={['' str{selection(1)} '']};
       for server=2:length(selection)
           str_group=[str_group, { Servers{server} } ];
       end

        [selection_group,ok_group]=listdlg('PromptString','Server for group analysis','SelectionMode','Single','ListString',str_group);

        if ok_group
            handles=guidata(hObject);

            handles.group_server=str_group{selection_group};
            guidata(hObject,handles);
        else
           set(handles.parallelize,'Value',0.0); 
        end

    end
    end
end
set(handles.create_batch,'enable','on');

function TR_Callback(hObject, eventdata, handles)
% hObject    handle to TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TR as text
%        str2double(get(hObject,'String')) returns contents of TR as a double


if ~isnan(str2double(get(handles.TR,'string')))
    set(handles.TRText,'string',get(handles.TR,'string'));
else
    set(handles.TRText,'string','NaN');
end


% --- Executes during object creation, after setting all properties.
function TR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gsr.
function gsr_Callback(hObject, eventdata, handles)
% hObject    handle to gsr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gsr

if get(handles.gsr,'Value')
    set(handles.loadmask,'Enable','Off');
    set(handles.load_second_mask,'Visible','Off');
    set(handles.parallelize,'Visible','Off');
    set(handles.chooselag,'Enable','Off');
    set(handles.m_min,'Enable','Off');
    set(handles.r_min,'Enable','Off');
    set(handles.TR,'Enable','Off');
    set(handles.gsr_brainmasks,'Visible','On');
    set(handles.BrainmasksPanel,'Visible','On');
    set(handles.compute_lagthreads,'Visible','Off');
    set(handles.text24,'Visible','Off');
    set(handles.grp1_name,'Visible','Off');
    set(handles.SleepTime,'Visible','Off');
    set(handles.text33,'Visible','Off');
    set(handles.use_second_mask,'Enable','Off');
    set(handles.PerformGSR,'Visible','On');
    set(handles.RunBatch,'Visible','Off');
else
   set(handles.gsr_brainmasks,'Visible','Off'); 
   set(handles.BrainmasksPanel,'Visible','Off');
end


% --- Executes on button press in gsr_brainmasks.
function gsr_brainmasks_Callback(hObject, eventdata, handles)
% hObject    handle to gsr_brainmasks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'RecentBrainmaskPath')
    [gsr_brainmasks_filenames,gsr_brainmasks_path]=uigetfile(['' handles.RecentBrainmaskPath '*.nii'],'MultiSelect','On','Select brainmasks');
else
    [gsr_brainmasks_filenames,gsr_brainmasks_path]=uigetfile(['' pwd filesep '*.nii'],'MultiSelect','On','Select brainmasks');
end

if strfind(gsr_brainmasks_path,' ');
   msgbox('Make sure there are no space characters in the GSR mask''s directory.');
   error('Make sure there are no space characters in the GSR mask''s directory.'); 
end

if ~isequal(gsr_brainmasks_path,0)
    handles=guidata(hObject);
    handles.RecentBrainmaskPath=gsr_brainmasks_path;
    guidata(hObject,handles);
else
    return
end

handles=guidata(hObject);
handles.gsr_brainmasks_filenames=gsr_brainmasks_filenames;
handles.gsr_brainmasks_path=gsr_brainmasks_path;
handles.gsr_brainmask=1; % TF function whether or not brainmasks are selected correctly.
guidata(hObject,handles);

GSRInfoString=cell(length(handles.BOLDfilename),1);
for file=1:length(handles.BOLDfilename)
    if length(gsr_brainmasks_filenames)==length(handles.BOLDfilename)
        GSRInfoString{file}=strcat(handles.BOLDfilename{file},' : ',gsr_brainmasks_path,gsr_brainmasks_filenames{file});
    else
        GSRInfoString{file}=strcat(handles.BOLDfilename{file},' : No brainmask selected ');
    end
end
set(handles.GSRBrainmasksInfo,'string',GSRInfoString); 

if length(gsr_brainmasks_filenames)~=length(handles.BOLDfilename)
    handles=guidata(hObject);
    handles.gsr_brainmask=0; % TF function whether or not brainmasks are selected correctly.
    guidata(hObject,handles);
    error('You need to select the same number of brainmasks as nifti-files.'); 
end

% set(handles.gsrbrainmaskText,'string',strcat(gsr_brainmask_path,gsr_brainmask_filename));




% --- Executes on button press in compute_lagthreads.
function compute_lagthreads_Callback(hObject, eventdata, handles)
% hObject    handle to compute_lagthreads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of compute_lagthreads

% --- Executes during object creation, after setting all properties.
function RUN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chooselag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.


% --- Executes on selection change in GSRBrainmasksInfo.
function GSRBrainmasksInfo_Callback(hObject, eventdata, handles)
% hObject    handle to GSRBrainmasksInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GSRBrainmasksInfo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GSRBrainmasksInfo


% --- Executes during object creation, after setting all properties.
function GSRBrainmasksInfo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GSRBrainmasksInfo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SleepTime_Callback(hObject, eventdata, handles)
% hObject    handle to SleepTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SleepTime as text
%        str2double(get(hObject,'String')) returns contents of SleepTime as a double


% --- Executes during object creation, after setting all properties.
function SleepTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SleepTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in PerformGSR.
function PerformGSR_Callback(hObject, eventdata, handles)
% hObject    handle to PerformGSR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'gsr_brainmask')
    error('Please specify GSR brainmask!');
end

gsr(strcat(handles.BOLDfilepath,handles.BOLDfilename),strcat(handles.gsr_brainmasks_path,handles.gsr_brainmasks_filenames),handles.strg_dir,0);


% --- Executes on button press in RunBatch.
function RunBatch_Callback(hObject, eventdata, handles)
% hObject    handle to RunBatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'strg_dir') || (isfield(handles,'strg_dir') && ~exist(['' handles.strg_dir '']))
    msgbox('Choose Storage Directory.','Error','error');
else
    commandwindow;
    fprintf('\n');
    disp(['Storage Directory: ' handles.strg_dir '']);
    fprintf('\n');

    for file=1:length(handles.BOLDfilename)
        fprintf(['\n------------------\nComputing on file ' handles.BOLDfilename{file} '\n']);
        if get(handles.use_second_mask,'Value')
            local_extrema_xcorr(strcat(handles.BOLDfilepath,handles.BOLDfilename{file}),handles.mask,handles.secondmask,handles.lags,handles.rmin,handles.mmin,str2double(get(handles.TR,'string')),2,handles.strg_dir,0);
        else
            local_extrema_xcorr(strcat(handles.BOLDfilepath,handles.BOLDfilename{file}),handles.mask,handles.mask,handles.lags,handles.rmin,handles.mmin,str2double(get(handles.TR,'string')),2,handles.strg_dir,0);
        end
        fprintf('------------------\n');
    end
    fprintf('\n\nBatch Processing terminated.\n');
end
