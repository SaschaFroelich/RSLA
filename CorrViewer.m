function varargout = CorrViewer(varargin)
% CORRVIEWER MATLAB code for CorrViewer.fig
%      CORRVIEWER, by itself, creates a new CORRVIEWER or raises the existing
%      singleton*.
%
%      H = CORRVIEWER returns the handle to a new CORRVIEWER or the handle to
%      the existing singleton*.
%
%      CORRVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRVIEWER.M with the given input arguments.
%
%      CORRVIEWER('Property','Value',...) creates a new CORRVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CorrViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CorrViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CorrViewer

% Last Modified by GUIDE v2.5 17-Apr-2018 16:51:55

% AUTHOR: Sascha Frölich, sascha.froelich@gmail.com

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CorrViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @CorrViewer_OutputFcn, ...
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


% --- Executes just before CorrViewer is made visible.
function CorrViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CorrViewer (see VARARGIN)

% Choose default command line output for CorrViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if length(varargin)
    handles=guidata(hObject);
    handles.filepath=varargin{1};
    handles.filename=varargin{2};
    if nargin
       handles.skipuiget=1; 
    end
    guidata(hObject,handles);

    loadfile_Callback(hObject, eventdata, handles);
end

% UIWAIT makes CorrViewer wait for user response (see UIRESUME)
% uiwait(handles.crrltns);


% --- Outputs from this function are returned to the command line.
function varargout = CorrViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
plot(handles.axes1,linspace(-5,-5,11),[0,0,0,0,0,0,0,0,0,0,0]); 
plot(handles.axes2,linspace(-5,-5,11),[0,0,0,0,0,0,0,0,0,0,0]); 
drawnow;

roi1=str2double(get(handles.roi1,'string'));
roi2=str2double(get(handles.roi2,'string'));
file=strcat(handles.filepath,handles.filename);
maxlag=str2double(get(handles.maxlag,'string'));
fig=get(handles.checkbox1,'Value');
m=load(['' file '']);

%In the new version of LAGCORR, we only correlate two rois. So LAGCORR will
%be a matrix with 2 rows and as many columns as timesteps
LAGCORR=lagcorr(file,maxlag,fig,1,roi1,roi2);
x=linspace(-maxlag,maxlag,2*maxlag+1);

if ~isfield(m,'mask2')
    I1=find(m.ROIS_TIMECOURSES(1:m.noofROIS)==roi1);
    I2=find(m.ROIS_TIMECOURSES(1:m.noofROIS)==roi2);
    set(handles.TD,'string',['TD(' num2str(roi1) ',' num2str(roi2) ')=' num2str(m.TD(I1,I2)) ' TD(' num2str(roi2) ',' num2str(roi1) ')=' num2str(m.TD(I2,I1)) '']);
elseif isfield(m,'mask2')
    I1=find(m.ROIS_TIMECOURSES1(1:m.noofROIS_mask1)==roi1);
    I2=find(m.ROIS_TIMECOURSES2(1:m.noofROIS_mask2)==roi2);
    set(handles.TD,'string',['TD(' num2str(roi1) ',' num2str(roi2) ')=' num2str(m.TD(I1,I2)) '']);
end

set(handles.r_min,'string',['r(t=0)_min=' num2str(m.r_min) '']);
set(handles.m_min,'string',['m_min=' num2str(m.m_min) '']);
set(handles.maxlag_info,'string',['maxlag=' num2str(max(m.lags)) '']);


xtick=linspace(-maxlag,maxlag,2*maxlag+1);
plot(handles.axes1,x,LAGCORR(1,:)); 
title(handles.axes1,['xcorr(' num2str(roi1) ',' num2str(roi2) ')']);
set(handles.axes1,'XTick',xtick);
set(handles.axes1,'XGrid','On');
set(handles.axes1,'YGrid','On');

plot(handles.axes2,x,LAGCORR(2,:));
title(handles.axes2,['xcorr(' num2str(roi2) ',' num2str(roi1) ')']);
set(handles.axes2,'XTick',xtick);
set(handles.axes2,'XGrid','On');
set(handles.axes2,'YGrid','On');

% --- Executes during object creation, after setting all properties.
function file_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in loadfile.
function loadfile_Callback(hObject, eventdata, handles)
% hObject    handle to loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%   If the GUI was called via another GUI (e.g. GUI lagmapsanal), then a field
%   handles.skipuiget was created and set to 1 (see opening function).
if isfield(handles,'skipuiget')
    if handles.skipuiget==1
        handles=guidata(hObject);
        handles.skipuiget=0; 
        guidata(hObject,handles);
        
        filename=handles.filename;
        filepath=handles.filepath;
    else
        if isfield(handles,'RecentFilePath')
            [filename,filepath]=uigetfile(['' handles.RecentFilePath '*.mat'],'Choose mat-file containing TD.');
            handles=guidata(hObject);
            handles.RecentFilePath=filepath;
            guidata(hObject,handles);
        else
            [filename,filepath]=uigetfile(['' pwd filesep '*.mat'],'Choose mat-file containing TD.');
            handles=guidata(hObject);
            handles.RecentFilePath=filepath;
            guidata(hObject,handles);
        end
    end
else
    if isfield(handles,'RecentFilePath')
        [filename,filepath]=uigetfile(['' handles.RecentFilePath '*.mat'],'Choose mat-file containing TD.');
        handles=guidata(hObject);
        handles.RecentFilePath=filepath;
        guidata(hObject,handles);
    else
        [filename,filepath]=uigetfile(['' pwd filesep '*.mat'],'Choose mat-file containing TD.');
        handles=guidata(hObject);
        handles.RecentFilePath=filepath;
        guidata(hObject,handles);
    end
end

set(handles.pushbutton1,'Enable','Off');
%set(handles.loadfile,'Enable','Off');
set(handles.plottcs,'Enable','Off');
drawnow;
if filename == 0
    set(handles.pushbutton1,'Enable','On');
    set(handles.loadfile,'Enable','On');
    set(handles.plottcs,'Enable','On');
else
    
    matfile=strcat(filepath,filename);
    set(handles.file,'string',['' matfile '']);
    drawnow;
    m=load(['' matfile '']);
    
    set(handles.niftifile,'string',['' m.BOLDfile '']);
    drawnow;
    
    if exist(m.BOLDfile,'file')~=2
        choice=questdlg(['The nifti-file ' m.BOLDfile ' from which TD in ' matfile ' was created cannot be found. It was probably moved to a new loaction. You can update that information by manually selecting the corresponding nifti-file.'],'Nifti file not found','Select File','Cancel','Select File');
        if strcmp(choice,'Select File')
           [BOLDfilename,BOLDfilepath]=uigetfile(['' pwd filesep '*.nii'],'MultiSelect','Off','Select nifti-file'); 
           BOLDfile=strcat(BOLDfilepath,BOLDfilename);
           save(matfile,'BOLDfile','-append');
           disp('BOLD-file directory successfully updated.');
           m=load(['' matfile '']);
           set(handles.niftifile,'string',['' m.BOLDfile '']);
           drawnow;
        end
    end
    
    [data,pixdim,rotate,dtype,slice_code] = readnifti(m.BOLDfile);
    
	handles=guidata(hObject);
    handles.data=data;
	handles.filename=filename;
	handles.filepath=filepath;
    handles.niifile=m.BOLDfile;
    if ~isfield(m,'mask2')
        handles.mask=m.mask;
        set(handles.listbox_rois2,'Visible','Off')
        set(handles.text32,'Visible','Off')
    elseif isfield(m,'mask2')
        handles.mask1=m.mask1;
        handles.mask2=m.mask2;
        set(handles.listbox_rois2,'Visible','On')
        set(handles.text32,'Visible','On')
    end
	guidata(hObject,handles);
    
    plot(handles.axes1,linspace(-5,-5,11),[0,0,0,0,0,0,0,0,0,0,0]); 
    plot(handles.axes2,linspace(-5,-5,11),[0,0,0,0,0,0,0,0,0,0,0]); 
    drawnow;
   
    
    set(handles.nanpairs,'Value',1.0);
    
    if ~isfield(m,'mask2')
        rois=m.ROIS_TIMECOURSES(1:m.noofROIS);
        percentage=m.noofNaN/(m.noofROIS*m.noofROIS);
        percentage=percentage*100;
        s=sprintf(['Number of Rois: ' num2str(m.noofROIS) '.\nNumber of NaN-entries: ' num2str(m.noofNaN) ' (~' num2str(percentage) '%%).\nMethod: ' m.Method '']);
        set(handles.listbox_rois1,'String',num2str(rois'));
    elseif isfield(m,'mask2')
        rois1=m.ROIS_TIMECOURSES1(1:m.noofROIS_mask1);
        rois2=m.ROIS_TIMECOURSES2(1:m.noofROIS_mask2);
        percentage=m.noofNaN/(m.noofROIS_mask1*m.noofROIS_mask2);
        percentage=percentage*100;
        s=sprintf(['Number of Rois1: ' num2str(m.noofROIS_mask1) '. Number of Rois2: ' num2str(m.noofROIS_mask2) '.\nNumber of NaN-entries: ' num2str(m.noofNaN) ' (~' num2str(percentage) '%%).\nMethod: ' m.Method '']);
        set(handles.listbox_rois1,'String',num2str(rois1'));
        set(handles.listbox_rois2,'String',num2str(rois2'));
    end
    set(handles.maxlag,'String',['' num2str(max(m.lags)) '']);
    maxlag_val=str2double(get(handles.maxlag,'String'));
    set(handles.maxlag_sec,'string',['' num2str(maxlag_val*2) ' sec']);
    drawnow;
    %set(handles.fileinfo,'string',['Number of Rois: ' num2str(m.noofROIS) '.  Number of NaN-entries: ' num2str(m.noofNaN) '']);
    set(handles.fileinfo,'string',s);
    drawnow;
    TD=m.TD;
    [rows, columns]=size(TD);
    nanpairs=[];
    set(handles.r_min,'string',['r(t=0)_min=' num2str(m.r_min) '']);
    set(handles.m_min,'string',['m_min=' num2str(m.m_min) '']);
    set(handles.maxlag_info,'string',['maxlag=' num2str(max(m.lags)) '']);
    drawnow;
    
    disp('Investigating NaN-pairs.');
    set(handles.nanpairs,'string','Loading. This may take a while.');
    drawnow;
    tic
    %Load only first 10.000
        for row=1:rows
            for col=1:columns
                if isnan(TD(row,col)) && length(nanpairs)<=10000
                    if isfield(m,'mask2')
                        nanpairs=[nanpairs; rois1(row), rois2(col)];
                    else
                        nanpairs=[nanpairs; rois(row), rois(col)];
                    end
                elseif length(nanpairs)>10000
                    row=rows;
                    break
                end
            end
        end
    toc
    
    
    set(handles.nanpairs,'string',num2str(nanpairs));
    drawnow;
    
    set(handles.pushbutton1,'Enable','On');
    set(handles.loadfile,'Enable','On');
    set(handles.loadstruct,'Enable','On');
    set(handles.plottcs,'Enable','On');
    drawnow;
   
end


% --- Executes during object creation, after setting all properties.
function roi1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function roi2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function roi1_Callback(hObject, eventdata, handles)
set(handles.roi1,'Value',1)
set(handles.roi2,'Value',0.0)


function maxlag_Callback(hObject, eventdata, handles)
% hObject    handle to maxlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxlag as text
%        str2double(get(hObject,'String')) returns contents of maxlag as a double
maxlag_val=str2double(get(handles.maxlag,'String'));
set(handles.maxlag_sec,'string',['' num2str(maxlag_val*2) ' sec']);
set(handles.maxlag,'string',maxlag_val);


% --- Executes during object creation, after setting all properties.
function maxlag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function maxlag_sec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxlag_sec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function roi2_Callback(hObject, eventdata, handles)
% hObject    handle to roi2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of roi2 as text
%        str2double(get(hObject,'String')) returns contents of roi2 as a double
set(handles.roi1,'Value',0.0);
set(handles.roi2,'Value',1);


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on selection change in listbox_rois.
function listbox_rois1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_rois contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_rois
value1=get(handles.roi1,'Value');
value2=get(handles.roi2,'Value');

contents = cellstr(get(handles.listbox_rois1,'String'));
roi=contents{get(handles.listbox_rois1,'Value')};

if value1
    set(handles.roi1,'string',roi);
elseif value2
    set(handles.roi2,'string',roi);
end
drawnow;

% --- Executes during object creation, after setting all properties.
function listbox_rois_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_rois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on selection change in nanpairs.
function nanpairs_Callback(hObject, eventdata, handles)
% hObject    handle to nanpairs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns nanpairs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from nanpairs
contents = cellstr(get(handles.nanpairs,'String'));
rois=contents{get(hObject,'Value')};
rois=strsplit(rois,' ');
if isempty(rois{1})
    roi1=rois{2};
    roi2=rois{3};
else
    roi1=rois{1};
    roi2=rois{2};
end
set(handles.roi1,'string',roi1);
set(handles.roi2,'string',roi2);
drawnow;

% --- Executes during object creation, after setting all properties.
function nanpairs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nanpairs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to maxlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxlag as text
%        str2double(get(hObject,'String')) returns contents of maxlag as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showTD.
function showTD_Callback(hObject, eventdata, handles)
% hObject    handle to showTD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showTD
value=get(handles.showTD,'Value');
file=get(handles.file,'string');
m=load(['' file '']);
if value
    set(handles.axes2,'Visible','off');
    set(handles.TDtable,'Visible','on');
else
    set(handles.axes2,'Visible','on');
    set(handles.TDtable,'Visible','off');
end
set(handles.TDtable,'Data',m.TD);
set(handles.TDtable,'ColumnName',m.ROIS_TIMECOURSES(1:m.noofROIS));
set(handles.TDtable,'RowName',m.ROIS_TIMECOURSES(1:m.noofROIS));


% --- Executes on button press in loadstruct.
function loadstruct_Callback(hObject, eventdata, handles)
% hObject    handle to loadstruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=handles.filename;
filepath=handles.filepath;


if iscell(filename)
    for iter=1:length(filename)
        matfile=strcat(filepath,filename{iter});
        file=load(['' matfile ''])
    end
else
    matfile=strcat(filepath,filename);
    file=load(['' matfile ''])
end

disp(['m=load(''' matfile ''')']);


% --- Executes on button press in plottcs.
function plottcs_Callback(hObject, eventdata, handles)
% hObject    handle to plottcs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi1=str2double(get(handles.roi1,'string'));
roi2=str2double(get(handles.roi2,'string'));
file=strcat(handles.filepath,handles.filename);

pushbutton1_Callback(hObject, eventdata, handles);

fig=get(handles.checkbox1,'Value');
m=load(['' file '']);

I1=find(m.ROIS_TIMECOURSES1(1:m.noofROIS_mask1)==roi1);
I2=find(m.ROIS_TIMECOURSES2(1:m.noofROIS_mask2)==roi2);

set(handles.TD,'string',['TD(' num2str(roi1) ',' num2str(roi2) ')=' num2str(m.TD(I1,I2)) ' TD(' num2str(roi2) ',' num2str(roi1) ')=' num2str(m.TD(I2,I1)) '']);

%Close existing Figure windows if desired
if get(handles.closefig,'value')
	fh=findobj(0,'type','figure');
	nfh=length(fh);
	for i=1:nfh
        if ~isempty(fh(i).Number)
            num=fh(i).Number;
            close(num);
        end
	end
end

x=linspace(1,length(m.ROIS_TIMECOURSES1(I1,2:end)),length(m.ROIS_TIMECOURSES1(I1,2:end)));
figure;
plot(x,m.ROIS_TIMECOURSES1(I1,2:end)); 
title(['ROI ' num2str(roi1) '']);
%set('XTick',xtick);
%set('XGrid','On');
%set('YGrid','On');
%annotation('textbox',[0.5 0.5 0.1 0.11],'string','blahblah');
movegui(gcf,'northwest');
figure;
plot(x,m.ROIS_TIMECOURSES2(I2,2:end));
title(['ROI ' num2str(roi2) '']);
%set('XTick',xtick);
%set('XGrid','On');
%set('YGrid','On');

disp(['roi_' num2str(roi1) '']);
disp(['roi_' num2str(roi2) '']);

assignin('base',['roi_' num2str(roi1) ''],m.ROIS_TIMECOURSES1(I1,2:end));
assignin('base',['roi_' num2str(roi2) ''],m.ROIS_TIMECOURSES2(I2,2:end));


% --- Executes on button press in closefig.
function closefig_Callback(hObject, eventdata, handles)
% hObject    handle to closefig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of closefig



function i_Callback(hObject, eventdata, handles)
% hObject    handle to i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of i as text
%        str2double(get(hObject,'String')) returns contents of i as a double


% --- Executes during object creation, after setting all properties.
function i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function j_Callback(hObject, eventdata, handles)
% hObject    handle to j (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of j as text
%        str2double(get(hObject,'String')) returns contents of j as a double


% --- Executes during object creation, after setting all properties.
function j_CreateFcn(hObject, eventdata, handles)
% hObject    handle to j (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function k_Callback(hObject, eventdata, handles)
% hObject    handle to k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of k as text
%        str2double(get(hObject,'String')) returns contents of k as a double


% --- Executes during object creation, after setting all properties.
function k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in loadnii.
function loadnii_Callback(hObject, eventdata, handles)
% hObject    handle to loadnii (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,filepath]=uigetfile(['' pwd filesep '*.nii']);

if filename == 0
    
else
    niifile=strcat(filepath,filename);
    set(handles.niftifile,'string',['' niifile '']);
    drawnow;

    %[niidata,pixdim,rotate,dtype,slice_code] = readnifti(niftifile);

    handles=guidata(hObject);
    %handles.niidata=niidata;
    handles.niifile=niifile;
    %handles.niifilepath=niifilepath;
    guidata(hObject,handles);
end



function nii_i_Callback(hObject, eventdata, handles)
% hObject    handle to nii_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nii_i as text
%        str2double(get(hObject,'String')) returns contents of nii_i as a double


% --- Executes during object creation, after setting all properties.
function nii_i_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nii_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nii_j_Callback(hObject, eventdata, handles)
% hObject    handle to nii_j (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nii_j as text
%        str2double(get(hObject,'String')) returns contents of nii_j as a double


% --- Executes during object creation, after setting all properties.
function nii_j_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nii_j (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nii_k_Callback(hObject, eventdata, handles)
% hObject    handle to nii_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nii_k as text
%        str2double(get(hObject,'String')) returns contents of nii_k as a double


% --- Executes during object creation, after setting all properties.
function nii_k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nii_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in plotniivox.
function plotniivox_Callback(hObject, eventdata, handles)
% hObject    handle to plotniivox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[niidata,pixdim,rotate,dtype,slice_code] = readnifti(handles.niifile);

data=niidata;
i=str2double(get(handles.nii_i,'string'));
j=str2double(get(handles.nii_j,'string'));
k=str2double(get(handles.nii_k,'string'));

%fsl vs matlab
i=i+1;
j=j+1;
k=k+1;

[rows,columns,pages,timesteps]=size(data);

for t=1:timesteps
    timecourse(t)=data(i,j,k,t);
end

%Close existing Figure windows if desired
if get(handles.closefig,'value')
	fh=findobj(0,'type','figure');
	nfh=length(fh);
	for iter=1:nfh
        if ~isempty(fh(iter).Number)
            num=fh(iter).Number;
            close(num);
        end
	end
end

%Determine which ROI the voxel belongs to
if isfield(handles,'mask')
    mask=handles.mask;
    [data_mask,pixdim_mask,orient_mask,dtype_mask,slice_code_mask] = readnifti(['' mask '']);
    ROI=data_mask(i,j,k);
    set(handles.roitxt,'string',['' num2str(ROI) '']);
else
    set(handles.roitxt,'string','no mask specified');
end

figure;
plot(linspace(1,timesteps,timesteps),timecourse); 
if isfield(handles,'mask')
    title(['i=' num2str(i-1) ', j=' num2str(j-1) ', k=' num2str(k-1) ' (ROI ' num2str(ROI) ')']);
else
    title(['i=' num2str(i-1) ', j=' num2str(j-1) ', k=' num2str(k-1) '']);
end

disp(['niitimecourse_' num2str(i-1) '_' num2str(j-1) '_' num2str(k-1) '']);

assignin('base',['niitimecourse_' num2str(i-1) '_' num2str(j-1) '_' num2str(k-1) ''],timecourse);

drawnow;


% --- Executes during object creation, after setting all properties.
function listbox_rois1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_rois1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_rois2.
function listbox_rois2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_rois2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_rois2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_rois2


% --- Executes during object creation, after setting all properties.
function listbox_rois2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_rois2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
