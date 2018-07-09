function varargout = rslaanalyzer(varargin)
% RSLAANALYZER MATLAB code for rslaanalyzer.fig
%      RSLAANALYZER, by itself, creates a new RSLAANALYZER or raises the existing
%      singleton*.
%
%      H = RSLAANALYZER returns the handle to a new RSLAANALYZER or the handle to
%      the existing singleton*.
%
%      RSLAANALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RSLAANALYZER.M with the given input arguments.
%
%      RSLAANALYZER('Property','Value',...) creates a new RSLAANALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before rslaanalyzer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to rslaanalyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help rslaanalyzer

% Last Modified by GUIDE v2.5 24-Apr-2018 09:25:59

% AUTHOR: Sascha Frölich, sascha.froelich@gmail.com

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @rslaanalyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @rslaanalyzer_OutputFcn, ...
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


% --- Executes just before rslaanalyzer is made visible.
function rslaanalyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to rslaanalyzer (see VARARGIN)

% Choose default command line output for rslaanalyzer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes rslaanalyzer wait for user response (see UIRESUME)
% uiwait(handles.lgmpsnl);


% --- Outputs from this function are returned to the command line.
function varargout = rslaanalyzer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in loadfile.
function loadfile_Callback(hObject, eventdata, handles)
% hObject    handle to loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'Recentfilepath')
	[filename,filepath]=uigetfile(['' handles.Recentfilepath '*.mat'],'MultiSelect','Off','Select groupfile (.mat)');
else
	[filename,filepath]=uigetfile(['' pwd filesep '*.mat'],'MultiSelect','Off');
end

if ~isequal(filepath,0)
	handles=guidata(hObject);
	handles.Recentfilepath=filepath;
	guidata(hObject,handles);
else
    return
end

set(handles.file,'string','Loading. Please wait...');
drawnow;

%Check if several files were selected or not
if iscell(filename) % several files selected (This should never happen)
    set(handles.performpca,'Enable','Off');
    handles=guidata(hObject);
    handles.filename=filename;
    handles.filepath=filepath;
    guidata(hObject,handles);
    set(handles.loadstruct,'Enable','On');
    set(handles.correlationviewer,'Enable','On');
    set(handles.file,'string',['' num2str(length(filename)) ' files selected.']);
elseif ~iscell(filename) %one file selected
    if filename==0 %not a single file selected
        %Leave descritptions as they were
        matfile=strcat(handles.filepath,handles.filename);
        set(handles.file,'string',['' matfile '']);
    else %Update GUI display
        set(handles.performpca,'Enable','Off');
        set(handles.file,'string','Loading, please wait.');
        set(handles.loadstruct,'Visible','On');
        set(handles.loadstruct,'Enable','On');
        set(handles.plotlagproj,'Enable','On');
        set(handles.plotlagmap,'Enable','On');
        set(handles.plottd,'Enable','On');
        set(handles.checktrans,'Enable','On');
        set(handles.checktrans,'Enable','On');
        set(handles.sumofthreads,'Enable','On');
        set(handles.nthreads,'Enable','On');
        matfile=strcat(filepath,filename);
        m=load(matfile);

        if strcmp(m.mask1,m.mask2)
            set(handles.choose_mask_buttongroup,'Visible','Off');
            set(handles.first_mask_lagproj,'Visible','Off');
            set(handles.second_mask_lagproj,'Visible','Off');
        else
            set(handles.choose_mask_buttongroup,'Visible','On');
            set(handles.first_mask_lagproj,'Visible','On');
            set(handles.second_mask_lagproj,'Visible','On');
        end

        if isfield(m,'ROIS_TIMECOURSES1')
            set(handles.correlationviewer,'Visible','On');
            set(handles.correlationviewer,'Enable','On'); 
        else
            set(handles.correlationviewer,'Visible','Off');
            set(handles.correlationviewer,'Enable','Off'); 
        end

        set(handles.mask2_panel,'Visible','On');
        set(handles.mask2_info,'Visible','On');
        set(handles.mask2_info,'string',m.mask2);
        if isfield(m,'BOLDfile')
            set(handles.NiftiFileInfofield,'string',m.BOLDfile);
        else
            set(handles.NiftiFileInfofield,'string','');
        end
        
        if isfield(m,'filenames') % Group level matlab file
            set(handles.stats,'enable','on');
            set(handles.savenewfile,'enable','on');
        else
            set(handles.stats,'enable','off');
            set(handles.savenewfile,'enable','off');
        end

        handles=guidata(hObject);
        handles.filename=filename;
        handles.filepath=filepath;
        handles.filestruct=m;
        % Set ploting=3 by default. This assingns a value to
        % "handles.plotting" and can be used when actually plotting
        handles.plotting=3; %0=lagthread 1=lagmaps 2=lag projection
        guidata(hObject,handles);


        %Fileinfos------------------------------------------------%
        percentage=m.noofNaN/(m.noofROIS_mask1*m.noofROIS_mask2);
        percentage=percentage*100;
        set(handles.uipanel4,'title','TD info'); 
        set(handles.noofrois,'string',['# ROIs (1st mask): ' num2str(m.noofROIS_mask1) '    # ROIs (2nd mask): ' num2str(m.noofROIS_mask2) '']);
        set(handles.noofnan,'string',['Number of NaN-entries: ' num2str(m.noofNaN) '(~' num2str(percentage) '%)']);
        set(handles.r_min,'string',['r_min=' num2str(m.r_min) '']);
        set(handles.m_min,'string',['m_min=' num2str(m.m_min) '']);
        set(handles.maxlag,'string',['maxlag=' num2str(max(m.lags)) '']);
        set(handles.file,'string',matfile);
        set(handles.mask,'string',m.mask1);

        if isfield(m,'ColumnMeans') %TDz has already been computed.
            set(handles.plotlagmap,'Enable','On');
        else %TDz has not yet been computed
            set(handles.performpca,'Enable','On');
        end

        if isfield(m,'COMPONENTS') %PCA has already been performed
            %set(handles.performpca,'Enable','Off');
            set(handles.performpca,'Enable','On');
            set(handles.evs,'Enable','On');
            set(handles.plotlagthread,'Enable','On');
        else
            set(handles.performpca,'Enable','On');
            set(handles.evs,'Enable','Off');
            set(handles.plotlagthread,'Enable','Off');
        end

        if issymmetric(m.TD,'skew')
            set(handles.fileinfo2,'string','TD is skew-symmetric');
        else
            set(handles.fileinfo2,'string','TD is not skew-symmetric, or contains NaN-entries.');
        end
        %/Fileinfos------------------------------------------------%
        
    end
end


% --- Executes on button press in performpca.
function performpca_Callback(hObject, eventdata, handles)
% hObject    handle to performpca (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%[comp,loadings,eigenvalues]=pca(TDz) returns the principal components of TDz as
%column vectors in matrix comp, with the leftmost column being the
%eigenvector corresponding to the greatest eigenvalue. The Column vector
%eigenvalues contains the non-zero eigenvalues.
filename=handles.filename;
filepath=handles.filepath;
if ~iscell(filename) %One file selected
    matfile=strcat(filepath,filename);
    
    storeredund=get(handles.storeredund,'Value');

	set(handles.performpca,'Enable','Off');
	drawnow;
    performpca(matfile,storeredund,0,1);
    
    m=load(matfile);
    handles=guidata(hObject);
    handles.filestruct=m;
    guidata(hObject,handles);
        
    set(handles.performpca,'string','Perform PCA');
    set(handles.performpca,'Enable','On');

    set(handles.plotlagmap,'Enable','On'); 
    set(handles.plotlagproj,'Enable','On');
    set(handles.plotlagthread,'Enable','On');
    set(handles.evs,'Enable','On');
    drawnow;

end


% --- Executes on button press in loadstruct.
function loadstruct_Callback(~, ~, handles)
% hObject    handle to loadstruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=handles.filename;
filepath=handles.filepath;
commandwindow
fprintf('Loading File. Please wait...\n');

if iscell(filename)
    for iter=1:length(filename)
        matfile=strcat(filepath,filename{iter});
        file=load(matfile)
    end
else
    matfile=strcat(filepath,filename);
    file=load(matfile)
end

disp(['m=load(''' matfile ''')']);
assignin('base','m',file);


% --- Executes on button press in storeredund.
%function storeredund_Callback(hObject, eventdata, handles)
% hObject    handle to storeredund (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of storeredund

function i_Callback(hObject, eventdata, handles)
% hObject    handle to i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'string') returns contents of i as text
%        str2double(get(hObject,'string')) returns contents of i as a double


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

% Hints: get(hObject,'string') returns contents of j as text
%        str2double(get(hObject,'string')) returns contents of j as a double


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

% Hints: get(hObject,'string') returns contents of k as text
%        str2double(get(hObject,'string')) returns contents of k as a double


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


% --- Executes on button press in plotlagmap.
function plotlagmap_Callback(hObject, eventdata, handles)
% hObject    handle to plotlagmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

filename=handles.filename;
filepath=handles.filepath;

file=strcat(handles.filepath,handles.filename);
m=handles.filestruct;

% If plotting==3, then this callback function was called directly from the
% button in the GUI (and not some other function in this script), i.e. plot
% the LAGMAP -> set plotting=1
plotting=handles.plotting; %0=lagthread 1=lagmaps 2=lag projection

if plotting==3
   plotting=1;
end

if plotting==0 % lagthread (i.e. Thread Cascade)
    if isnan(str2double(get(handles.lagthread,'string')))
        msgbox('Enter Lag-Thread No.','Error','error');
    else

        fprintf('\n!!!!!!! \n NB: Plotting column of m.L(). Make sure weighted sum of lag-threads allows for this (Alert 1).\n!!!!!!!\n');

        rois=m.TDrois_rows;

        seedthread=str2double(get(handles.lagthread,'string'));
        lagplot=m.L(:,seedthread);

        set(handles.variance,'string',['' num2str(round(m.explained(str2double(get(handles.lagthread,'string'))),1)) '% Variance explained.']);

    end
elseif plotting==1 %lagmap (i.e. seeded cascade). A Lag map is a column vector either of TD or TDz
    if isnan(str2double(get(handles.lagmapseed,'string')))
        msgbox('Enter Lag-Map seed (Error 1)','Error','error');
    else
        % The SEED is from a different mask than the ROIs that are
        % displayed!
        seed_intens=str2double(get(handles.lagmapseed,'string'));

        if get(handles.second_mask_lagproj,'Value')
            % rois is the list of ROIs on which to be plotted.
            rois=m.TDrois_rows;
            seedthread=find(seed_intens==m.TDrois_cols);
            disp(['Plotting TD(' num2str(seedthread) ',:)'])
            lagplot=m.TD(:,seedthread);
        elseif get(handles.first_mask_lagproj,'Value') % This is never the case
            rois=m.TDrois_rows;
            seedthread=find(seed_intens==m.TDrois_cols);
            disp(['Plotting TD(:,' num2str(seedthread) ')'])
            lagplot=m.TD(:,seedthread);
        end
    end
elseif plotting==2 % lag projection (i.e. Global Cascade) of TD
        % We need vectors "rois" and "lagplot".
        % The ROIs are the row-ROIs of TD.

        rois=m.TDrois_cols;
        if m.level==2
            lagplot=m.GroupLagProjection;
        elseif m.level==1
            if isfield(m,'ColumnMeans')
                lagplot=m.ColumnMeans;
            else
                lagplot=nanmean(m.TD);
            end
        end
        
        seedthread=0;
end

islice=str2num(get(handles.i,'string'));
jslice=str2num(get(handles.j,'string'));
kslice=str2num(get(handles.k,'string'));

if ~strcmp(m.mask1,m.mask2)
    if length(rois)==m.noofROIS_mask1
        maskfile=m.mask1;
    elseif length(rois)==m.noofROIS_mask2
        maskfile=m.mask2;
    end
else
    maskfile=m.mask1;
end

ifind=strfind(get(handles.i,'string'),':');
jfind=strfind(get(handles.j,'string'),':');
kfind=strfind(get(handles.k,'string'),':');

% Only if one view (sagittal, coronal or transveral) is plotted!
% Otherwise the option to plot several slices at once is not possible!
if ~isempty(ifind)
    itemp=get(handles.i,'string');
    if length(ifind)==1
        islice=str2double(itemp(1:ifind(1)-1)):str2double(itemp(ifind(1)+1:end));
    elseif length(ifind)==2
        islice=str2double(itemp(1:ifind(1)-1)):str2double(itemp(ifind(1)+1:ifind(2)-1)):str2double(itemp(ifind(2)+1:end));
    end
    jslice=str2double(get(handles.j,'string'));
    kslice=str2double(get(handles.k,'string'));
elseif ~isempty(jfind)
    jtemp=get(handles.j,'string');
    if length(jfind)==1
        jslice=str2double(itemp(1:jfind(1)-1)):str2double(jtemp(jfind(1)+1:end));
    elseif length(jfind)==2
        jslice=str2double(jtemp(1:jfind(1)-1)):str2double(jtemp(jfind(1)+1:jfind(2)-1)):str2double(jtemp(jfind(2)+1:end));
    end
    islice=str2double(get(handles.i,'string'));
    kslice=str2double(get(handles.k,'string'));
elseif ~isempty(kfind)
    ktemp=get(handles.k,'string');
    if length(kfind)==1
        kslice=str2double(ktemp(1:kfind(1)-1)):str2double(ktemp(kfind(1)+1:end));
    elseif length(kfind)==2
        kslice=str2double(ktemp(1:kfind(1)-1)):str2double(ktemp(kfind(1)+1:kfind(2)-1)):str2double(ktemp(kfind(2)+1:end));
    end
    islice=str2double(get(handles.i,'string'));
    jslice=str2double(get(handles.j,'string'));
else
    islice=str2double(get(handles.i,'string'));
    jslice=str2double(get(handles.j,'string'));
    kslice=str2double(get(handles.k,'string'));
end

ploti=get(handles.plot_i,'Value');
plotj=get(handles.plot_j,'Value');
plotk=get(handles.plot_k,'Value');
slicepos=get(handles.slicepos,'Value');

% If a slice is negative (e.g. islice=-5), then it shall not be plotted, but
% its position (abs value) can still be used to position the crosshair
if ploti==0
   islice=islice(1)*(-1); 
end

if plotj==0
   jslice=jslice(1)*(-1);  
end

if plotk==0
   kslice=kslice(1)*(-1); 
end

%set absolute color range for every plot?
if get(handles.revert_plot,'Value')
    lagplot=lagplot*(-1);
end

for islice_iter=islice
    for jslice_iter=jslice
        for kslice_iter=kslice

            if get(handles.set_c_range,'Value');
                createplots('none','none',maskfile,islice_iter,jslice_iter,kslice_iter,rois,lagplot,plotting,seedthread,slicepos,{get(handles.set_c_range,'Value'),str2double(get(handles.c_range_from,'String')),str2double(get(handles.c_range_to,'String'))},strcat(handles.filepath,handles.filename),{get(handles.SavePlotsPNG,'Value'),get(handles.PNGDIRTXT,'string')});
            else
                createplots('none','none',maskfile,islice_iter,jslice_iter,kslice_iter,rois,lagplot,plotting,seedthread,slicepos,{get(handles.set_c_range,'Value'),0,0},strcat(handles.filepath,handles.filename),{get(handles.SavePlotsPNG,'Value'),get(handles.PNGDIRTXT,'string')});
            end

        end
    end
end

if get(handles.linear_vector_plot,'Value')
    figure;imagesc(lagplot); colormap('jet');set(gca,'Position',[.10 .11 .685 .815]); cb2 = colorbar('Position',[.81 .11 .0675 .815]);
end
assignin('base','lagplot',lagplot)
    
handles=guidata(hObject);
handles.plotting=3; % 0=lagthread 1=lagmaps 2=lag projection
guidata(hObject,handles);

% --- Executes on button press in plotlagthread.
function plotlagthread_Callback(hObject, eventdata, handles)
% hObject    handle to plotlagthread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guidata(hObject);
handles.plotting=0; %0=lagthread 1=lagmaps 2=lag projection
guidata(hObject,handles);

plotlagmap_Callback(hObject,eventdata,handles);


function lagmapseed_Callback(hObject, eventdata, handles)
% hObject    handle to lagmapseed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'string') returns contents of lagmapseed as text
%        str2double(get(hObject,'string')) returns contents of lagmapseed as a double


% --- Executes during object creation, after setting all properties.
function lagmapseed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lagmapseed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in loadfile2.
function loadfile2_Callback(hObject, eventdata, handles)
% hObject    handle to loadfile2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filenames,filepath2]=uigetfile(['' pwd filesep '*.mat'],'MultiSelect','On','Select matlab-files containing TD.');

%Check if several files were selected or not
if iscell(filenames) %several files selected
    set(handles.file2,'string','Loading, please wait.');
    set(handles.loadfile2,'enable','off');
    set(handles.loadstruct2,'Enable','off');
    set(handles.groupleveltd,'Enable','Off');
    drawnow;
    handles=guidata(hObject);
    handles.filenames=filenames;
    handles.filepath2=filepath2;
    guidata(hObject,handles);
    set(handles.loadstruct2,'Visible','On');
    drawnow;
    
    generateTDpart_tf=zeros(1,length(filenames)); %1 for those files that needs to have a TDpart created
    generateTDz_tf=zeros(1,length(filenames)); %1 for those files that needs to have a TDz created
    nocomp=0;
    for iter=1:length(filenames)
       file = strcat(filepath2,filenames{iter});
       m=load(file);
       if m.noofNaN && ~isfield(m,'nodes')
           generateTDpart_tf(iter)=1;
       end
       if ~isfield(m,'ColumnMeans')
           generateTDz_tf(iter)=1;
       end
       
       %Check whether all files have same r_min and m:min
       if iter==1
          m1=m.m_min;
          r1=m.r_min;
          m2=m1;
          r2=r1;
          mask1=m.mask1;
          mask2=m.mask2;
          mask1COMP=m.mask1;
          mask2COMP=m.mask2;
       else
           m2=m.m_min;
           r2=m.r_min;
           mask1COMP=m.mask1;
           mask2COMP=m.mask2;
       end
       
       if (m2~=m1 || r2~=r1) && nocomp==0
           msgbox('r_min and m_min are not identical over all selected files.','Error','error');
           %nocomp==2 -> r/m
           nocomp=2;
       end
       
       if ~get(handles.allowmasks,'Value')
           if  (~strcmp(mask1,mask1COMP) || ~strcmp(mask2,mask2COMP)) && nocomp==0
               msgbox('Brainmasks are not identical over all selected files.','Error','error');
               %nocomp==1 -> masks
               nocomp=1;
           end
       end
       
       
    end
    
	if nocomp==2 || (nocomp==1 && ~get(handles.allowmasks,'Value'))
        set(handles.groupleveltd,'Enable','Off');
        drawnow;
    else
        set(handles.groupleveltd,'Enable','On');
        drawnow;
    end
    
    handles=guidata(hObject);
    handles.generateTDpart_tf=generateTDpart_tf;
    handles.generateTDz_tf=generateTDz_tf;
    guidata(hObject,handles);
    
    set(handles.file2,'string',['' num2str(length(filenames)) ' files selected.']);
    set(handles.loadfile2,'enable','on');
    set(handles.loadstruct2,'Enable','On');

elseif ~iscell(filenames) %one file selected
    if filenames == 0
        
    else
        set(handles.file2,'string','1 file selected.');
    end
end


% --- Executes on button press in loadstruct2.
function loadstruct2_Callback(~, ~, handles)
% hObject    handle to loadstruct2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filenames=handles.filenames;
filepath2=handles.filepath2;
if iscell(filenames)
    for iter=1:length(filenames)
        file = strcat(filepath2,filenames{iter});
        matfile=load(file);
        commandwindow
        disp(['m=load(''' file ''')']);
    end
else
    file = strcat(filepath2,filenames);
    matfile=load(file)
    commandwindow
    disp(['m=load(''' file ''')']);
end


% --- Executes on button press in groupleveltd.
function groupleveltd_Callback(hObject, eventdata, handles)
% hObject    handle to groupleveltd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filenames=handles.filenames;
filepath2=handles.filepath2;

disp('Preparing for groupTD computation.');

commandwindow;
disp('Starting Computation.');

set(handles.loadfile2,'Enable','Off');
set(handles.loadstruct2,'Enable','off');
drawnow;

handles=guidata(hObject);
groupleveltd(filenames,filepath2,'GROUPFILE',filepath2,0);
grouplagproj('GROUPFILE.mat',filepath2,0);

set(handles.groupleveltd,'Enable','On');
set(handles.loadfile2,'Enable','On');
set(handles.loadstruct2,'Enable','On');
msgbox('Group-Level TD successfully computed.');


% --- Executes on button press in storeredund2.
%function storeredund2_Callback(hObject, eventdata, handles)
% hObject    handle to storeredund2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of storeredund2



%function choosefilename_Callback(hObject, eventdata, handles)
% hObject    handle to choosefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of choosefilename as text
%        str2double(get(hObject,'String')) returns contents of choosefilename as a double


% --- Executes during object creation, after setting all properties.
function choosefilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to choosefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in closefig.
function closefig_Callback(hObject, eventdata, handles)
% hObject    handle to closefig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of closefig



%function lagthread_Callback(hObject, eventdata, handles)
% hObject    handle to lagthread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lagthread as text
%        str2double(get(hObject,'String')) returns contents of lagthread as a double


% --- Executes during object creation, after setting all properties.
function lagthread_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lagthread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in evs.
function evs_Callback(~, ~, handles)
% hObject    handle to evs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=handles.filename;
filepath=handles.filepath;

if iscell(filename) %several files selected
        %This will never be the case as pltolagmap will only be enabled if a
        %single file is selected
elseif ~iscell(filename) %one file selected
    file=strcat(handles.filepath,handles.filename);
    m=load(file);
    
    if length(m.explained)<20
    
        figure;
        scatter(1:length(m.explained),m.explained,'filled');
        title('Variance explained.');
        ylabel('% explained');
        xlabel('Eigenvalue no.');
        
    else
    
        figure;
        scatter(1:20,m.explained(1:20),'filled');
        title('Variance explained.');
        ylabel('% explained');
        xlabel('Eigenvalue no.');
        
    end

end


% --- Executes on button press in plotlagproj.
function plotlagproj_Callback(hObject, eventdata, handles)
% hObject    handle to plotlagproj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles=guidata(hObject);
handles.plotting=2; %0=lagthread 1=lagmaps 2=lag projection
guidata(hObject,handles);

plotlagmap_Callback(hObject,eventdata,handles);



%function transrows1_Callback(hObject, eventdata, handles)
% hObject    handle to transrows1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transrows1 as text
%        str2double(get(hObject,'String')) returns contents of transrows1 as a double


% --- Executes during object creation, after setting all properties.
function transrows1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transrows1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%function transcols1_Callback(hObject, eventdata, handles)
% hObject    handle to transcols1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transcols1 as text
%        str2double(get(hObject,'String')) returns contents of transcols1 as a double


% --- Executes during object creation, after setting all properties.
function transcols1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transcols1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checktrans.
function checktrans_Callback(hObject, eventdata, handles)
% hObject    handle to checktrans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=handles.filename;
filepath=handles.filepath;

if ~iscell(filename)
    file=strcat(handles.filepath,handles.filename);
    m=handles.filestruct;

    [ntrans,trans]=transitive(m.TD);
    
    msgbox(['Transitive triples: ' num2str(trans) '. Non-transitive triples:' num2str(ntrans) '']);
    disp(['Transitive triples: ' num2str(trans) '. Non-transitive triples:' num2str(ntrans) '']);    
    
    if get(handles.savenewfile,'Value')
        
        save(file,'ntrans','-append');
        save(file,'trans','-append');
        
    end
    
end



%function transrows2_Callback(hObject, eventdata, handles)
% hObject    handle to transrows2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transrows2 as text
%        str2double(get(hObject,'String')) returns contents of transrows2 as a double


% --- Executes during object creation, after setting all properties.
function transrows2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transrows2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%function transcols2_Callback(hObject, eventdata, handles)
% hObject    handle to transcols2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of transcols2 as text
%        str2double(get(hObject,'String')) returns contents of transcols2 as a double


% --- Executes during object creation, after setting all properties.
function transcols2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to transcols2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in allowmasks.
%function allowmasks_Callback(hObject, eventdata, handles)
% hObject    handle to allowmasks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allowmasks


% --- Executes on slider movement.
function islider_Callback(hObject, eventdata, handles)
% hObject    handle to islider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
value=round(get(handles.islider,'value'));
set(handles.islider,'value',value);
set(handles.i,'string',num2str(value));


% --- Executes on button press in slicepos.
function slicepos_Callback(hObject, eventdata, handles)
% hObject    handle to slicepos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slicepos


% --- Executes on button press in grouplagproj.
function grouplagproj_Callback(~, ~, handles)
% hObject    handle to grouplagproj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tic
filenames=handles.filenames;
filepath2=handles.filepath2;

%The following code has been moved to grouplagproj.m

grouplagproj(filenames,filepath2);
toc     


% --- Executes on button press in plot_i.
function plot_i_Callback(hObject, eventdata, handles)
% hObject    handle to plot_i (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_i


% --- Executes on button press in plot_j.
function plot_j_Callback(hObject, eventdata, handles)
% hObject    handle to plot_j (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_j


% --- Executes on button press in plot_k.
function plot_k_Callback(hObject, eventdata, handles)
% hObject    handle to plot_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plot_k


% --- Executes on button press in correlationviewer.
function correlationviewer_Callback(~, ~, handles)
% hObject    handle to correlationviewer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CorrViewer(handles.filepath,handles.filename);


% --- Executes on button press in set_c_range.
function set_c_range_Callback(~, ~, handles)
% hObject    handle to set_c_range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of set_c_range
if get(handles.set_c_range,'Value')
    set(handles.c_range_to,'Enable','On');
    set(handles.c_range_from,'Enable','On');
else
    set(handles.c_range_to,'Enable','Off');
    set(handles.c_range_from,'Enable','Off');
end




function c_range_from_Callback(hObject, eventdata, handles)
% hObject    handle to c_range_from (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c_range_from as text
%        str2double(get(hObject,'String')) returns contents of c_range_from as a double


% --- Executes during object creation, after setting all properties.
function c_range_from_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c_range_from (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c_range_to_Callback(hObject, eventdata, handles)
% hObject    handle to c_range_to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c_range_to as text
%        str2double(get(hObject,'String')) returns contents of c_range_to as a double


% --- Executes during object creation, after setting all properties.
function c_range_to_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c_range_to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in plottd.
function plottd_Callback(~, ~, handles)
% hObject    handle to plottd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

file=strcat(handles.filepath,handles.filename);
m=handles.filestruct;

figure;

if get(handles.reversetd,'Value')
    imagesc(-m.TD);
else
    imagesc(m.TD);
end

colormap('jet')
set(gca,'Position',[.10 .11 .685 .815]);
colorbar('Position',[.81 .11 .0675 .815]);
title(handles.filename);

if get(handles.clrrange,'value')
   caxis([-1000,1000]); 
end


%function no_iterations_Callback(hObject, eventdata, handles)
% hObject    handle to no_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of no_iterations as text
%        str2double(get(hObject,'String')) returns contents of no_iterations as a double


% --- Executes during object creation, after setting all properties.
function no_iterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to no_iterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stats.
function stats_Callback(hObject, eventdata, handles)
% hObject    handle to stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Close existing Figure windows if desired
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

if get(handles.loadstats,'Value')
    
    m=handles.filestruct;
    
    if ~isfield(m,'TDConvCorr')
        commandwindow
       error('No TD convergence statistics in file.'); 
    end
    
    figure;
    if isfield(m,'crr_rand_method1')
        plot(1:m.max_subset,m.crr,'green',1:m.max_subset,m.crr_rand_method1,'red',1:m.max_subset,m.crr_rand_method2,'blue',1:m.max_subset,m.crr_rand_method3,'cyan');
    else
        plot(1:m.TDConvMaxSubset,m.TDConvCorr,'green'); 
    end
    
    msgbox(['non-transitive triples: ' num2str(m.ntrans) '']);
    
else
    
    no_its=str2double(get(handles.no_iterations,'string'));
    TD_corr(strcat(handles.filepath,handles.filename),no_its,get(handles.savenewfile,'value')); 
    checktrans_Callback(hObject,eventdata,handles)
    set(handles.stats,'Enable','On');
    
end


% --- Executes on button press in savenewfile.
%function savenewfile_Callback(hObject, eventdata, handles)
% hObject    handle to savenewfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of savenewfile


% --- Executes on button press in reversetd.
%function reversetd_Callback(hObject, eventdata, handles)
% hObject    handle to reversetd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reversetd


% --- Executes on button press in clrrange.
%function clrrange_Callback(hObject, eventdata, handles)
% hObject    handle to clrrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of clrrange


% --- Executes on button press in loadstats.
function loadstats_Callback(hObject, eventdata, handles)
% hObject    handle to loadstats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of loadstats


function nthreads_Callback(hObject, eventdata, handles)
% hObject    handle to nthreads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nthreads as text
%        str2double(get(hObject,'String')) returns contents of nthreads as a double


% --- Executes during object creation, after setting all properties.
function nthreads_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nthreads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sumofthreads.
function sumofthreads_Callback(hObject, eventdata, handles)
% hObject    handle to sumofthreads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename=handles.filename;
filepath=handles.filepath;

file=strcat(handles.filepath,handles.filename);
m=load(file);

if ~isfield(m,'EIGENVALUES');
   error('You need to perform PCA on TD before this step.'); 
end

no_threads=str2double(get(handles.nthreads,'string'));

% This will store the sign of the Lag Threads, and the associated
% correlation value in the last column.
CorrPerms=zeros(no_threads+1,no_threads+1);

% The vector "signs" stores the signs for the sum (+ or - 1).
sign=ones(1,no_threads);

for no_minus=0:no_threads
    
    % Yes, that's intentional! (first loop: "for minus=1:0", i.e. loop is not entered)
    for minus=1:no_minus
        sign(minus)=-1;
    end
    
        MinusPerms=perms(sign);
        MinusPerms=unique(MinusPerms,'rows');
    
        B=bsxfun(@times,m.EIGENVALUES(1:no_threads)',MinusPerms);
       
        WeightedSum=m.L(:,[1:no_threads])*B';     % Columns of WeightedSums contain the weighted sum for all possible combination of minuses
        
        CorrVal=zeros(1,size(MinusPerms,1));
        
        for ii=1:size(MinusPerms,1)
           CorrVal(ii)=corr2(WeightedSum(:,ii),m.GroupLagProjection');
        end
        
        % Select the maximum correlation value for all possible
        % combinations of no_minus minuses
        [~,I]=max(abs(CorrVal));
        CorrPerms(no_minus+1,:)=[MinusPerms(I,:) CorrVal(I)];    % Store the order of minuses along with the associated maximum correlation value
	
end

LagThread=CorrPerms(:,1:end-1);
CorrelationValue=CorrPerms(:,end);

VariableNames=cell(1,no_threads+1);
VariableNames{end}='CorrelationValue';

for ii=1:no_threads
   VariableNames{ii}=['LagThread' num2str(ii) ''];
end
commandwindow
array2table([LagThread CorrelationValue],'VariableNames',VariableNames)

% --- Executes on button press in linear_vector_plot.
%function linear_vector_plot_Callback(hObject, eventdata, handles)
% hObject    handle to linear_vector_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of linear_vector_plot


% --- Executes on button press in difference_plot.
%function difference_plot_Callback(hObject, eventdata, handles)
% hObject    handle to difference_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of difference_plot


% --- Executes on button press in pushbutton40.
%function pushbutton40_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton40 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in psychometrics.
function psychometrics_Callback(hObject, eventdata, handles)
% hObject    handle to psychometrics (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Close existing Figure windows if desired
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

% choice=questdlg('Use groupfile (.mat) or select number of .nii-files?','Where to extract latency values from?','Groupfile (slow)','Nii (fast)','Nii (fast)');
choice='Nii (fast)';
switch choice
                case 'Groupfile (slow)'
                    disp(['Continuing with Groupfile. ']);
                    if isfield(handles,'RecentgroupfilePath')
                        [groupfileName,groupfilePath]=uigetfile(['' handles.RecentgroupfilePath '*.mat'],'MultiSelect','Off','Select groupfile (.mat)');
                        handles=guidata(hObject);
                        handles.RecentgroupfilePath=groupfilePath;
                        guidata(hObject,handles);
                    else
                        [groupfileName,groupfilePath]=uigetfile(['' pwd filesep '*.mat'],'MultiSelect','Off','Select groupfile (.mat)');
                        handles=guidata(hObject);
                        handles.RecentgroupfilePath=groupfilePath;
                        guidata(hObject,handles);
                    end
                    groupfile=strcat(groupfilePath,groupfileName);
                    
                    if isfield(handles,'RecentcsvfilePath')
                        [csvfileName,csvfilePath]=uigetfile(['' handles.RecentcsvfilePath '.csv'],'MultiSelect','off','Select .csv-file containing behavioral scores.');
                        handles=guidata(hObject);
                        handles.RecentcsvfilePath=csvfilePath;
                        guidata(hObject,handles);
                    else
                        [csvfileName,csvfilePath]=uigetfile(['' pwd filesep '*.csv'],'MultiSelect','Off','Select .csv-file containing behavioral scores.');
                        handles=guidata(hObject);
                        handles.RecentcsvfilePath=csvfilePath;
                        guidata(hObject,handles);
                    end
                    csvfile=strcat(csvfilePath,csvfileName)
                    
                    if get(handles.SaveCorrFig,'Value')
                        if isfield(handles,'RecentDIRECTORYNAME')
                            DIRECTORYNAME=uigetdir(['' handles.RecentDIRECTORYNAME ''],'Where to store correlation figures?');
                            handles=guidata(hObject);
                            handles.RecentDIRECTORYNAME=DIRECTORYNAME;
                            guidata(hObject,handles);
                        else
                            DIRECTORYNAME=uigetdir(['' pwd ''],'Where to store correlation figures?');
                            handles=guidata(hObject);
                            handles.RecentDIRECTORYNAME=DIRECTORYNAME;
                            guidata(hObject,handles);
                        end
                    end
                    
                    
                    fid = fopen(csvfile); % 'Data' worksheet

                    clear contents
                    contents=cell(1,1);
                    contentsTemp=cell(1,1);
                    EscLoop=0;
                    line=0;
                    while EscLoop~=1
                        line=line+1;
                        if line==1
                            contentsTemp{line,1}=fgetl(fid);
                        else
                            contentsTemp=[contentsTemp; fgetl(fid)];
                        end
                        if contentsTemp{line,1}==-1
                           EscLoop=1;
                           contents=contentsTemp(1:end-1,:);
                           clear contentsTemp
                        end
                    end

                    acols = length(find(contents{1,1}==','))+1; % number of columns
                    aformat = repmat('%s ', 1, acols); % create format
                    fclose(fid);
                    
                    fid = fopen(csvfile);
                    temp = textscan(fid,aformat,'Delimiter',',');
                    fclose(fid);
                    
                    ColumnNames=cell(1,size(temp,2));
                    no_lines=size(contents,1);
                    
                    disp(['No. of subjects: ' num2str(no_lines) '']);
                    
                    for clm=1:size(temp,2)
                       Column=temp{clm};
                       ColumnNames(clm)=Column(1);
                    end
                    [selection_temp,ok]=listdlg('PromptString','Choose one or more Psychometric(s)','SelectionMode','Single','ListString',ColumnNames,'SelectionMode','multiple');
                    
                    if ok
                        [nii_maskfilename,nii_maskpath]=uigetfile(['' pwd filesep '*.nii'],'MultiSelect','Off','Select mask-file (need not be binary)');
                        NiiMaskfile=strcat(nii_maskpath,nii_maskfilename);
                        [MaskFileData,~,~,~,~] = readnifti(NiiMaskfile);
                        m=load(groupfile);

                        %%% The following header parts should be identical for
                        %%% MaskFileData and the nii-files:
                        % pixdim
                        % qoffset_x
                        % srwo_x
                        % If that is not the case, they are likely to have
                        % flipped orientation
                        MaskHeader=readniftifileheader(NiiMaskfile);

                        NiiHeader=readniftifileheader(m.mask2);
                        [niidata,~,~,~,~] = readnifti(m.mask2);
                        if ~isequal(NiiHeader.srow_x,MaskHeader.srow_x) || ~isequal(NiiHeader.pixdim,MaskHeader.pixdim) || ~isequal(NiiHeader.qoffset_x,MaskHeader.qoffset_x)
                            msgbox('Header of ROI-Maskfile different from header of nifti-files. Check if orientations are compatible.','Error','error');
                            error('Header of ROI-Maskfile different from header of nifti-files. Check if orientations are compatible.');
                        elseif size(niidata)~=size(MaskFileData)
                            msgbox('Nifti-files and ROI-Maskfile do not have the same matrix siz.','Error','error');
                            error('Nifti-files and ROI-Maskfile do not have the same matrix size.');
                        end

                        MaskFileData(find(MaskFileData>0))=1;
                        RoisOvrlp=MaskFileData.*niidata;
                        roisUn=unique(RoisOvrlp);
                        if roisUn(1)==0
                            roisUn=roisUn(2:end);
                        end
                        rois_temp=cell(1,1);
                        rois_temp{1}=strcat(num2str(roisUn(1)));

                        if length(roisUn)>1
                            for idxi=2:length(roisUn)
                               rois_temp{1}=strcat(rois_temp{1},',',num2str(roisUn(idxi)));
                            end
                        end
                        for idxi=1:length(roisUn)
                            weighting(idxi)=length(find(RoisOvrlp==roisUn(idxi)));
                        end
                        
                        for selection=selection_temp
                            % roilistcontents=cellstr(get(handles.tdroiscols,'String'));
                            % roi=str2double(roilistcontents{get(handles.tdroiscols,'Value')});
                            clear semicolons
                            semicolons=strfind(rois_temp,';');

                            if isequal(semicolons,{[]})
                               semicolons=[]; 
                            end

                            for k=1:length(semicolons)+1
                                if isempty(semicolons)
                                    RoisTMP=strsplit(rois_temp{1},',');
                                    rois=[];
                                    for entry=1:length(RoisTMP)
                                       if ~isempty(RoisTMP(entry))
                                           rois=[rois, str2double(RoisTMP(entry))];
                                       end
                                    end

                                else
                                    rois_k=strsplit(rois_temp{1},';');
                                    rois_k_temp=strsplit(rois_k{k},',');
                                    rois=[];
                                    for entry=1:length(rois_k_temp)
                                       if  ~isempty(rois_k_temp(entry))
                                           rois=[rois, str2double(rois_k_temp(entry))];
                                       end
                                    end
                                end
                            
                                disp(['For psychometric measure ' ColumnNames{selection} ' using the following roi list: ' num2str(rois) '']);
                                
                                scores=cell(no_lines-1,2);
                                
                                for subj=2:no_lines
                                    ContentsLine=strsplit(contents{subj,1},',','CollapseDelimiters',false);
                                    scores{subj-1,1}=ContentsLine{1};
                                    scores{subj-1,2}=str2double(ContentsLine{selection});
                                end

                                disp('The subjects'' names in the csv-file have to be equal to the filenames in the selected groupfile');
                                latency=[];
                                psychometric_score=[];
                                AllSbjcts={};
                                for file=1:length(m.filenames)

                                    temp_filename=m.filenames{file};
                                    temp_filename=temp_filename(1:end-4);

                                    if ~isnan(scores{find(strcmp(temp_filename,scores)),2})

                                       SubjectTDfile=load(strcat(m.pathtofiles,m.filenames{file})); 

                                       % disp(['Checking file ' strcat(m.pathtofiles,m.filenames{file}) ''])

                                       I=[];
                                       for roi=1:length(rois)
                                            I=[I find(rois(roi)==SubjectTDfile.TDrois_cols)];
                                       end
                                       ColumnMeans=nanmean(SubjectTDfile.TD);

                                       lag_values=[];
                                       
                                        switch ROI_choice
                                            case 'Enter ROIs manually'
                                               for i=1:length(I)
                                                    lag_values=[lag_values ColumnMeans(I(i))];
                                               end
                                                latency=[latency nanmean(lag_values,2)];
                                               
                                            case 'Load Mask'
                                               for i=1:length(I)
                                                   % This works as vector
                                                   % rois() is sorted like
                                                   % vector weighting()
                                                    lag_values=[lag_values ColumnMeans(I(i))*weighting(i)];
                                               end
                                               
                                                not_nans=zeros(size(lag_values,1),1); % stores the no. of NOT-NaNs in every row of lag_values

                                                for RW=1:size(lag_values,1)
                                                    not_nans(RW)=length(find(~isnan(lag_values(RW,:))));
                                                end
                                                
                                                latency=[latency nanmean(lag_values,2).*not_nans/sum(weighting)];
                                                
                                        end
                                        psychometric_score=[psychometric_score scores{find(strcmp(temp_filename,scores)),2}];
                                        AllSbjcts=[AllSbjcts temp_filename];
                                    end

                                end

                                [r,p]=corr(psychometric_score',latency');

                                figure;set(gcf,'name',['' groupfile '']); scatter(psychometric_score,latency,'filled'); title(['' num2str(length(rois)) ' ROIs (lag) vs ' ColumnNames{selection} ' (r=' num2str(r) ', p=' num2str(p) ')']);
                                xlabel(['' ColumnNames{selection} ''])
                                ylabel('lag');

                                print(gcf, '-dpng', ['' DIRECTORYNAME '' filesep 'sel' num2str(selection) 'roiset' num2str(k) '']);
                                close(gcf);
                                
                                assignin('base','latency',latency)
                                assignin('base','psychometric_score',psychometric_score)
                                assignin('base','AllSbjcts',AllSbjcts)
                            end
                        end
                    end

                        set(handles.psychometrics,'Enable','On');
                        disp('Done!');
                        
        case 'Nii (fast)' % Was tested - works
            % set(handles.psychometrics,'Enable','Off');
                    if isfield(handles,'Recentnii_filepath')
                        [nii_filenames,nii_filepath]=uigetfile(['' handles.Recentnii_filepath '*.nii'],'MultiSelect','On','Select nifti-files');
                    else
                        [nii_filenames,nii_filepath]=uigetfile(['' pwd filesep '*.nii'],'MultiSelect','On','Select nifti-files');
                    end
                    
                    if ~isequal(nii_filepath,0)
                        handles=guidata(hObject);
                        handles.Recentnii_filepath=nii_filepath;
                        guidata(hObject,handles);
                    else
                       return 
                    end
                    
                    if isfield(handles,'RecentcsvfilePath')
                        [csvfileName,csvfilePath]=uigetfile(['' handles.RecentcsvfilePath '*.csv'],'MultiSelect','On','Select .csv-files containing behavioral scores.');
                    else
                        [csvfileName,csvfilePath]=uigetfile(['' pwd filesep '*.csv'],'MultiSelect','On','Select .csv-files containing behavioral scores.');
                    end
                    
                    if ~isequal(csvfilePath,0)
                        handles=guidata(hObject);
                        handles.RecentcsvfilePath=csvfilePath;
                        guidata(hObject,handles);
                    else
                       return 
                    end
                    
                    csvfile=strcat(csvfilePath,csvfileName)
                    
                    fid = fopen(csvfile); % 'Data' worksheet

                    clear contents
                    contents=cell(1,1);
                    contentsTemp=cell(1,1);
                    EscLoop=0;
                    line=0;
                    while EscLoop~=1
                        line=line+1;
                        if line==1
                            contentsTemp{line,1}=fgetl(fid);
                        else
                            contentsTemp=[contentsTemp; fgetl(fid)];
                        end
                        if contentsTemp{line,1}==-1
                           EscLoop=1;
                           contents=contentsTemp(1:end-1,:);
                           clear contentsTemp
                        end
                    end

                    acols = length(find(contents{1,1}==','))+1; % number of columns
                    aformat = repmat('%s ', 1, acols); % create format
                    fclose(fid);
                    
                    fid = fopen(csvfile);
                    temp = textscan(fid,aformat,'Delimiter',',');
                    fclose(fid);
                    
                    ColumnNames=cell(1,size(temp,2));
                    no_lines=size(contents,1);
                    
                    disp(['No. of subjects: ' num2str(no_lines) '']);
                    
                    for clm=1:size(temp,2)
                       Column=temp{clm};
                       ColumnNames(clm)=Column(1);
                    end
                    %%% Choice = Nii (fast)
                    [selection_temp,ok]=listdlg('PromptString','Choose one or more Psychometric(s)','SelectionMode','Single','ListString',ColumnNames,'SelectionMode','multiple');
                    
                    if ok
                            ROI_choice='Load Mask';
                            
                            if isfield(handles,'Recentnii_maskpath')
                                [nii_maskfilename,nii_maskpath]=uigetfile(['' handles.Recentnii_maskpath filesep '*.nii'],'MultiSelect','Off','Select ROI mask-file (need not be binary)');
                            else
                                [nii_maskfilename,nii_maskpath]=uigetfile(['' pwd filesep '*.nii'],'MultiSelect','Off','Select ROI mask-file (need not be binary)');
                            end
                            if ~isequal(nii_maskpath,0)
                                handles=guidata(hObject);
                                handles.Recentnii_maskpath=nii_maskpath;
                                guidata(hObject,handles);
                            else
                               return 
                            end

                            NiiMaskfile=strcat(nii_maskpath,nii_maskfilename);
                            [MaskFileData,~,~,~,~] = readnifti(NiiMaskfile);

                            % set(handles.psychometrics,'Enable','Off');

                            %%% The following header parts should be identical for
                            %%% MaskFileData and the nii-files:
                            % pixdim
                            % qoffset_x
                            % srwo_x
                            % If that is not the case, they are likely to have
                            % flipped orientation
                            MaskHeader=readniftifileheader(NiiMaskfile);

                            NiiHeader=readniftifileheader(strcat(nii_filepath,nii_filenames{1}));
                            [niidata,~,~,~,~] = readnifti(strcat(nii_filepath,nii_filenames{1}));
                            if ~isequal(NiiHeader.srow_x,MaskHeader.srow_x) || ~isequal(NiiHeader.pixdim,MaskHeader.pixdim) || ~isequal(NiiHeader.qoffset_x,MaskHeader.qoffset_x)
                                msgbox('Header of ROI-Maskfile different from header of nifti-files. Check if orientations are compatible.','Error','error');
                                error('Header of ROI-Maskfile different from header of nifti-files. Check if orientations are compatible.');
                            elseif size(niidata)~=size(MaskFileData)
                                msgbox('Nifti-files and ROI-Maskfile do not have the same matrix siz.','Error','error');
                                error('Nifti-files and ROI-Maskfile do not have the same matrix size.');
                            end

                            MaskFileData(find(MaskFileData>0))=1;
                            clear niidata
                        
                            CorrMode=questdlg('Analysis method?','What analysis do you want to perform?','Two-Sample Pearson Correlation','Partial Correlation','Two-Sample Pearson Correlation');
                            
                            if get(handles.SaveCorrFig,'Value')
                                if isfield(handles,'RecentDIRECTORYNAME')
                                    DIRECTORYNAME=uigetdir(['' handles.RecentDIRECTORYNAME ''],'Where to store correlation figures?');
                                    handles=guidata(hObject);
                                    handles.RecentDIRECTORYNAME=DIRECTORYNAME;
                                    guidata(hObject,handles);
                                else
                                    DIRECTORYNAME=uigetdir(['' pwd ''],'Where to store correlation figures?');
                                    handles=guidata(hObject);
                                    handles.RecentDIRECTORYNAME=DIRECTORYNAME;
                                    guidata(hObject,handles);
                                end
                            end
                        
                        switch CorrMode
                            case 'Partial Correlation'
                                
                                % Iterate through psychometric measures
                                scores=cell(no_lines-1,length(selection_temp));

                                for subj=2:no_lines
                                    ContentsLine=strsplit(contents{subj,1},',','CollapseDelimiters',false);
                                    scores{subj-1,1}=ContentsLine{1};
                                    for jj=1:length(selection_temp)
                                        scores{subj-1,jj+1}=str2double(ContentsLine{selection_temp(jj)});
                                    end
                                end
                                

                                disp('The subjects'' names as given in the csv-file must be present as string-pattern in the selected .nii-filenames (i.e. if a subject in the csv-file is "subjM006", then the corresponding nifty-file could be "swarsubjM006lag_projection.nii")');
                                latency=[];
                                psychometric_score=[];
                                AllSbjcts={};
                                for file=1:length(nii_filenames)

                                    temp_filename=nii_filenames{file};
                                    temp_filename=temp_filename(1:end-4);

                                    RowFind=0;
                                    for ROW=1:size(scores,1)
                                       if find(strfind(temp_filename,scores{ROW,1}))
                                           RowFind=ROW;
                                       end
                                    end

                                    % Make sure RowFind is not equal to
                                    % zero (this happens if a nifti-file is chosen which does not have a corresponding row in the csv file.)
                                    if RowFind
                                        if sum(isnan([scores{RowFind,2:end}]))==0

                                           [SubjectNiiData,~,~,~,~] = readnifti(strcat(nii_filepath,nii_filenames{file}));

                                           % disp(['Checking file ' strcat(nii_filepath,nii_filenames{file}) ''])

                                           % MaskFileData has already
                                           % been binarized
                                           MaskFileOvrlp=MaskFileData.*SubjectNiiData;
                                           MaskFileNotZero=find(MaskFileData>0);

                                            latency=[latency nanmean(MaskFileOvrlp(MaskFileNotZero))];
                                            temp_score=[scores{RowFind,2:end}]';
                                            psychometric_score=[psychometric_score temp_score];
                                            AllSbjcts=[AllSbjcts, temp_filename];
                                        end
                                    else
                                       disp(['NO ROW IN CSV FILE CORRESPONDING TO ' temp_filename '']); 
                                    end

                                end

                                X=[latency; psychometric_score]';
                                [rho,pval]=partialcorr(X);
                                commandwindow
                                fprintf('\n');
                                disp('------------------------------');
                                disp(['ROI maskfile: ' strcat(nii_maskpath,nii_maskfilename) '']);
                                disp(['Files-folder: ' nii_filepath '']);

                                str='';
                                for SEL=1:length(selection_temp)-1
                                    str=strcat(str,ColumnNames{selection_temp(SEL)},', ');
                                end
                                str=strcat(str,ColumnNames{selection_temp(end)});

                                disp(['Measures: ' str '']);
                                for kk=2:size(rho,1)
                                    fprintf('\n');
                                    disp(['Partial correlation of latency with measure ' ColumnNames{selection_temp(kk-1)} ' correcting for all other ' num2str(size(rho,1)-2) ' measures.']);
                                    disp(['r= ' num2str(rho(kk,1)) '']);
                                    disp(['p= ' num2str(pval(kk,1)) '']);

                                    if get(handles.SaveCorrFig,'Value')

                                        figure;set(gcf,'name','Partial Correlation from nii-files'); 
                                        scatter(psychometric_score(kk-1,:),latency,'filled'); 
                                        title(['Partial Corr of ' ColumnNames{selection_temp(kk-1)} ' vs ' nii_maskfilename ' (r=' num2str(rho(kk,1)) ', p=' num2str(pval(kk,1)) ')']);
                                        xlabel(['' ColumnNames{selection_temp(kk-1)} ''])
                                        ylabel('lag');

                                        print(gcf, '-dpng', ['' DIRECTORYNAME '' filesep 'PartCorr_' ColumnNames{selection_temp(kk-1)} 'VS' nii_maskfilename '']);
                                        close(gcf);

                                    end

                                end
                                fprintf('\n');
                                disp('------------------------------');
                                commandwindow

                                %disp('Spearman:');

                                %[rho,pval]=partialcorr(X,'type','Spearman')

                            case 'Two-Sample Pearson Correlation'

                                % Iterate through psychometric measures
                                for selection=selection_temp

                                    % scores will contain the
                                    % corresponding psychometric score
                                    % of ALL subjects (as listed in
                                    % csv-file). Those needed are then
                                    % filtered by use of the selected
                                    % .nii-files, where the name of the
                                    % subjects in the csv-file musst be
                                    % present as a string in the
                                    % filenames.
                                    scores=cell(no_lines-1,2);

                                    for subj=2:no_lines
                                            % If one cell was empty, this is presented by two adjacent commas
                                            % in the csv-file ('a,,b'). ')'CollapseDelimiters',false accounts for that.
                                            ContentsLine=strsplit(contents{subj,1},',','CollapseDelimiters',false);
                                            scores{subj-1,1}=ContentsLine{1};
                                            scores{subj-1,2}=str2double(ContentsLine{selection});

                                    end

                                    disp('The subjects'' names as given in the csv-file must be present as string-pattern in the selected .nii-filenames (i.e. if a subject in the csv-file is "subjM006", then the corresponding nifty-file could be "swarsubjM006lag_projection.nii")');
                                    latency=[];
                                    psychometric_score=[];
                                    AllSbjcts={};
                                    for file=1:length(nii_filenames)

                                        temp_filename=nii_filenames{file};
                                        temp_filename=temp_filename(1:end-4);

                                        RowFind=0;
                                        for ROW=1:size(scores,1)
                                           if find(strfind(temp_filename,scores{ROW,1}))
                                               RowFind=ROW;
                                           end
                                        end
                                        
                                        % Make sure RowFind is not equal to
                                        % zero (this happens if a nifti-file is chosen which does not have a corresponding row in the csv file.)
                                        if RowFind
                                            if ~isnan(scores{RowFind,2})

                                               [SubjectNiiData,~,~,~,~] = readnifti(strcat(nii_filepath,nii_filenames{file}));

                                               % disp(['Checking file ' strcat(nii_filepath,nii_filenames{file}) ''])

                                               % MaskFileData has already
                                               % been binarized
                                               MaskFileOvrlp=MaskFileData.*SubjectNiiData;
                                               MaskFileNotZero=find(MaskFileData>0);

                                              	latency=[latency nanmean(MaskFileOvrlp(MaskFileNotZero))];
                                                psychometric_score=[psychometric_score scores{RowFind,2}];
                                                AllSbjcts=[AllSbjcts, temp_filename];
                                            end
                                        else
                                            disp(['NO ROW IN CSV FILE CORRESPONDING TO ' temp_filename '']); 
                                        end

                                    end

                                    % use [r,p]=corr(a,b)
                                    [r,p]=corr(psychometric_score',latency');

                                    fprintf('\n');
                                    disp('------------------------------');
                                    disp(['Measure: ' ColumnNames{selection} '']);
                                    disp(['ROI-mask: ' nii_maskfilename '']);
                                    disp(['File-folder: ' nii_filepath '']);
                                    fprintf(['r=' num2str(r) ',\np=' num2str(p) '\n']);
                                    disp('------------------------------');
                                    fprintf('\n');
                                    % disp('Spearman:');

                                    % [r,p]=corr(psychometric_score',latency','type','Spearman')

                                    commandwindow

                                    if get(handles.SaveCorrFig,'Value')

                                        figure;set(gcf,'name','Correlation from nii-files'); 
                                        scatter(psychometric_score,latency,'filled'); 
                                        title(['' ColumnNames{selection} ' vs ' nii_maskfilename ' (r=' num2str(r) ', p=' num2str(p) ')']);
                                        xlabel(['' ColumnNames{selection} ''])
                                        ylabel('lag');

                                        print(gcf, '-dpng', ['' DIRECTORYNAME '' filesep '' ColumnNames{selection} 'VS' nii_maskfilename '']);
                                        close(gcf);

                                    end

                                    assignin('base','latency',latency)
                                    VariableName=ColumnNames{selection};
                                    % - and space cannot be part of a Variable Name
                                    VariableName(strfind(VariableName,'-'))='_';
                                    VariableName(strfind(VariableName,' '))='_';
                                    assignin('base',['' VariableName ''],psychometric_score)
                                    assignin('base','AllSbjcts',AllSbjcts)

                                end
                        end
                    end
                    
                    set(handles.psychometrics,'Enable','On');
            
end


% --- Executes on button press in revert_plot.
function revert_plot_Callback(hObject, eventdata, handles)
% hObject    handle to revert_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of revert_plot


% --- Executes on button press in smooth.
function smooth_Callback(hObject, eventdata, handles)
% hObject    handle to smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%%%%%%%%%
% NB: THIS FUNCTION CREATES AN UNSMOOTHED VERSION OF EVERY NIFTI-FILE PLUS A
% SMOOTHED VERSION
%%%%%%%%%%%%%%

clc;

%set(handles.smooth,'Enable','Off');
if isfield(handles,'Recentfilepath')
	[filename,filepath]=uigetfile(['' handles.Recentfilepath '*.mat'],'MultiSelect','On','Select mat-files.');
else
	[filename,filepath]=uigetfile(['' pwd filesep '*.mat'],'MultiSelect','On','Select mat-files.');
end

if ~isequal(filepath,0)
	handles=guidata(hObject);
	handles.Recentfilepath=filepath;
	guidata(hObject,handles);
else
    return
end


if get(handles.groupnif,'Value');
   files=cell(1,1); 
   files_smooth=cell(1,1); 
end

if ~get(handles.smooth_lagmap,'Value')
    if isfield(handles,'RecentDIRECTORYNAME')
        DIRECTORYNAME=uigetdir(['' handles.RecentDIRECTORYNAME ''],'Where do you want the generated nifti-files to be stored?')
    else
        DIRECTORYNAME=uigetdir(['' pwd ''],'Where do you want the generated nifti-files to be stored?')
    end

    if ~isequal(DIRECTORYNAME,0)
        handles=guidata(hObject);
        handles.RecentDIRECTORYNAME=DIRECTORYNAME;
        guidata(hObject,handles);
    else
        return
    end
end

if get(handles.smooth_lagproj,'Value') % Global cascade
    % uigetfile returns a cell if a number of files are selected,
    % while it returns a character array if one file is selected.
    % Make sure from this point we always have a cell array. The
    % length of the cell array corresponds to the number of selected
    % files
    if iscell(filename)
        filenames=filename;
    else
       filenames=cell(1,1);
       filenames{1,1}=filename;
    end
    
    commandwindow;
   for file=1:length(filenames)
        disp(['' num2str(length(filenames)-file+1) ' file(s) to go.']);
        disp(['Next file: ' strcat(filepath,filenames{file}) '']);
        m=load(strcat(filepath,filenames{file}));

        %Load corresponding brainmask in MNI-space:
        if ~strcmp(m.mask1,m.mask2)
            commandwindow
            error('Using two different brainmasks','error','error');
        else
            [niidata,~,~,~,~] = readnifti(m.mask1);

            %niidata_binary is a matrix with 1's where niidata has a
            %value different from 0, otherwise it's 0.
            niidata_binary=niidata;
            niidata_binary_ref=niidata;

            niidata_binary(find(niidata~=0))=1;
            niidata_binary(find(niidata==0))=0;

            niidata_binary_ref(find(niidata~=0))=1;
            niidata_binary_ref(find(niidata==0))=0;

        end

        % rois_idx contains the columns in TD which correspond to the
        % lag projection

        % TD_rois_cols contains the list of rois along the second
        % dimension of TD (as many ROIs as columns)
        rois_idx=m.TDrois_cols;


        % Walk through every entry of the lag projection. The ith
        % entry is the value which has to be assigned to all voxels in
        % that ROI, the ROI being m.TDrois_cols(i)
        niidata_temp=niidata;
        for roi=1:length(m.TDrois_cols)
            roi_intens=m.TDrois_cols(roi);

            % Remember: niidata contains the niftidata matrix from
            % m.mask1. Find all entries in the mask-matrix that
            % correspond to the ROI of the lagprojection we are
            % currently iterating over
            idx=find(roi_intens==niidata);

            if ~isempty(idx)
                niidata_temp(idx)=m.ColumnMeans(roi);
                niidata_binary(idx)=-1;
            else
                commandwindow
                error('idx empty','error','error');
            end

        end
        niidata=niidata_temp;
        clear niidata_temp

        %Now, if niidata_binary+niidata_binary_ref is zero everywhere, only the
        %entries where niidata was different from zero were changed.

        if length(find(niidata_binary+niidata_binary_ref~=0))
            commandwindow
            error('Difference larger than zero','error','error');
        end
        clear niidata_binary niidata_binary_ref

        current_filename=filenames{file};
        name=strcat(current_filename(1:end-4),'_LagProjection');

        [niidata_rows,niidata_columns,niidata_pages]=size(niidata);

        if niidata_rows~=91 || niidata_columns~=109 || niidata_pages~=91
            commandwindow
            error('niidata is of compromising matrix size. Only implemented for 2mmx2mmx2mm nifti files.');
        else
            write_nifti(m.mask1,strcat(DIRECTORYNAME,filesep,name,'.nii'),niidata,[3 91 109 91 0 0 0 0]);
        end

        disp(['Number of nan-entries: ' num2str(length(find(isnan(m.ColumnMeans)))) ' (of ' num2str(length(m.ColumnMeans)) ')']);

        sx=str2double(get(handles.sx,'string'))
        sy=str2double(get(handles.sy,'string'))
        sz=str2double(get(handles.sz,'string'))

        % Smooth the created nifti-file while creating a new,
        % smoothed file
        spm_smooth(strcat(DIRECTORYNAME,filesep,name,'.nii,1'),strcat(DIRECTORYNAME,filesep,['S' num2str(sx) num2str(sy) num2str(sz) ''],name,'.nii,1'),[str2double(get(handles.sx,'string')) str2double(get(handles.sy,'string')) str2double(get(handles.sz,'string'))],0);
        disp('Smoothing done.');

        if file==1
           files{1}=strcat(name,'.nii');
           files_smooth{1}=strcat(['S' num2str(sx) num2str(sy) num2str(sz) ''],name,'.nii');
        else
            files=[files strcat(name,'.nii')];
            files_smooth=[files_smooth, strcat(['S' num2str(sx) num2str(sy) num2str(sz) ''],name,'.nii')];
        end
   end
%Smooth lag map
elseif get(handles.smooth_lagmap,'Value')
    % Load a Mask-File. Every voxel which is not zero will be
    % considered
    if isfield(handles,'RecentMapSeedMaskPath')
        [MapSeedMasknames,MapSeedMaskPath]=uigetfile(['' handles.RecentMapSeedMaskPath '*.nii'],'MultiSelect','On','Load Seed-ROI-Masks for Lag-Maps (need not be binary).');
    else
        [MapSeedMasknames,MapSeedMaskPath]=uigetfile(['' pwd filesep '*.nii'],'MultiSelect','On','Load Seed-ROI-Masks for Lag-Maps (need not be binary).');
    end

    if ~isequal(MapSeedMaskPath,0)
        handles=guidata(hObject);
        handles.RecentMapSeedMaskPath=MapSeedMaskPath;
        guidata(hObject,handles);
    else
        return
    end

    % Make MapSeedMasknames into 1x1 cell if it is not a cell
    % to begin with
    if ~iscell(MapSeedMasknames)
        MapSeedMasknames={MapSeedMasknames};
    end

    % Choose storage directories for every single mask
    StrgDIRECTORIES=cell(length(MapSeedMasknames),1);
    Dir=0;
    for MapSeedMasknameITR=MapSeedMasknames
        MapSeedMaskname=MapSeedMasknameITR{1};
        Dir=Dir+1;
        if isfield(handles,'RecentStrgDIRECTORIES')
            DIRECTORYNAME=uigetdir(['' handles.RecentStrgDIRECTORIES ''],['Storage-Dir for ' MapSeedMaskname ' files?'])
        else
            DIRECTORYNAME=uigetdir(['' pwd ''],['Storage-Dir for ' MapSeedMaskname ' files?'])
        end

        if ~isequal(DIRECTORYNAME,0)
            handles=guidata(hObject);
            handles.RecentStrgDIRECTORIES=DIRECTORYNAME;
            guidata(hObject,handles);
            StrgDIRECTORIES{Dir,1}=DIRECTORYNAME;
        else
            return
        end
    end
    clear Dir

    MaskIter=0;
    for MapSeedMasknameITR=MapSeedMasknames
        MapSeedMaskname=MapSeedMasknameITR{1};
        disp(['Next Lag-Map Seed: ' MapSeedMaskname '']);
        MaskIter=MaskIter+1;

        MapSeed=strcat(MapSeedMaskPath,MapSeedMaskname);

        file=strcat(StrgDIRECTORIES{MaskIter,1},filesep,'maskfile.txt');
        fiD=fopen(file,'w');
        fprintf(fiD,['' MapSeed '']);
        fclose(fiD);

        % Check if chosen seedfile is compatible with mask
        SeedMaskHeader=readniftifileheader(MapSeed);

        if iscell(filename)
            ChckFilename=filename{1};
        else
            ChckFilename=filename;
        end


        ChckFile=load(strcat(filepath,ChckFilename));
        ChckHeader=readniftifileheader(ChckFile.mask2);

        if iscell(filename)
            for VerificationFileNo=2:length(filename)
                VerificationFile=load(strcat(filepath,filename{VerificationFileNo}));
                if ~isequal(VerificationFile.mask1,ChckFile.mask1) || ~isequal(VerificationFile.mask2,ChckFile.mask2)
                    commandwindow
                    error(['File ' filename{VerificationFileNo} ' does not have same mask as first file (' ChckFilename '). NB: This error might occur if the mask-files are indeed the same but stored as two files with diverging names or in two different directories. In this case, this error message can be commented out before trying again.'],'error','error');
                end
            end
        end


        [SeedMaskData,~,~,~,~] = readnifti(MapSeed);
        [ChckData,~,~,~,~] = readnifti(ChckFile.mask2);

        if ~isequal(ChckHeader.srow_x,SeedMaskHeader.srow_x) || ~isequal(ChckHeader.pixdim,SeedMaskHeader.pixdim) || ~isequal(ChckHeader.qoffset_x,SeedMaskHeader.qoffset_x)
            msgbox('Header of mask (i.e. nifti file of seed region) is different from header of nifti-files. Check if orientations are compatible.','Error','error');
            error('Header of mask is different from header of nifti-files. Check if orientations are compatible.');
        elseif size(SeedMaskData)~=size(ChckData)
            msgbox('Nifti-files and mask (i.e. nifti file of seed region) do not have the same matrix siz.','Error','error');
            error('Nifti-files and mask do not have the same matrix size.');
        end

        % Seed-ROIs are given in terms of mask2. The number of
        % voxels per ROI that are actually in the seed-file
        % determines the weighting of the specific ROI.

        % In OverlapData, only those entries are ~=0, which are ~=0
        % in the given seed mask. The value of those entries is the
        % corresponding ROI-entry from "mask2" (which is the mask
        % used to specify the columns of TD and was used for TD
        % computation).
        SeedMaskData(find(SeedMaskData~=0))=1;

        % ChckData is the nii-data belonging to the mask which
        % makes up the columns of the TD-matrix. 
        OverlapData=SeedMaskData.*ChckData;

        % Now code is compatible with later code
        % (which was originally coded for manual entry of seed regions)
        rois=unique(OverlapData)';

        if rois(1)==0
           rois=rois(2:end); 
        end

        % Now for every ROI in OverlapData, determine of how many
        % voxels it consists. Store the value in the weightind
        % vector WghtVec
        WghtVec=zeros(length(rois),1);
        for ii=1:length(rois)
            WghtVec(ii)=length(find(OverlapData==rois(ii)));
        end
        clear roi

        % Average the lag maps of a number of ROIs to obtain a single lag
        % map
        disp(['Using the following roi list: ' num2str(rois) '']);

        % uigetfile returns a cell if a number of files are selected,
        % while it returns a character array if one file is selected.
        % Make sure from this point we always have a cell array. The
        % length of the cell array corresponds to the number of selected
        % files
        if iscell(filename)
            filenames=filename;
        else
           filenames=cell(1,1);
           filenames{1,1}=filename;
        end

        m=load(strcat(filepath,filenames{1}));
        % TDrois_cols = rois_in_mask2
        if any(ismember(rois,m.TDrois_cols)==0)
                msgbox('The defined ROIs are not in mask2','Error','error');
        end

        for file=1:length(filenames)
            disp(['' num2str(length(filenames)-file+1) ' file(s) to go for this seed.']);
            disp(['Next file: ' strcat(filepath,filenames{file}) '']);
            m=load(strcat(filepath,filenames{file}));


            % Load corresponding brainmask in MNI-space:
            % LagMapMask = m.mask2;
            niidatamask= m.mask1;

            % niidata will be the data which will be written to the
            % newly generated nii-file. It is of the dimensions of the
            % mask other than the LagMapMask
            [niidata,~,~,~,~] = readnifti(niidatamask);

            % niidata_binary is a matrix with 1's where niidata has a
            % value different from 0, otherwise it's 0.
            niidata_binary=niidata;
            niidata_binary_ref=niidata;

            niidata_binary(find(niidata~=0))=1;
            niidata_binary(find(niidata==0))=0;

            niidata_binary_ref(find(niidata~=0))=1;
            niidata_binary_ref(find(niidata==0))=0;


            % rois_idx contains the columns (rows) in TD which correspond to the
            % lag maps
            rois_idx=[];

            for roi=1:length(rois)
                rois_idx=[rois_idx find(rois(roi)==m.TDrois_cols)]; 
            end

            % lagmaps stores all columns (rows) over which to average
            % next to each other as columns. lagmaps is thus a matrix
            % where the no of columns corresponds to the number of ROIs
            % over which to average, and the no of rows of lagmaps is
            % the no. of ROIs the resulting (averaged) lagmap will
            % contain.
            lagmaps=[];

            for roi=1:length(rois)
                lagmaps=[lagmaps, m.TD(:,rois_idx(roi))*WghtVec(roi)];
            end
            not_nans=zeros(size(lagmaps,1),1); % stores the no. of NOT-NaNs in every row of lagmaps

            for RW=1:size(lagmaps,1)
                not_nans(RW)=length(find(~isnan(lagmaps(RW,:))));
            end

            lagmap_avrg=nanmean(lagmaps,2).*not_nans/sum(WghtVec);

            % Now, lagmap_avrg is the averaged lagmap. 

            % Walk through every entry of the averaged lagmap. The ith
            % entry is the value which has to be assigned to all voxels in
            % that ROI, the ROI being m.TDrois_rows(i)
            niidata_temp=niidata;
            for roi=1:length(m.TDrois_rows)
                roi_intens=m.TDrois_rows(roi);

                % Find at which voxel-positions in the niidata-matrix of
                % the chosen seed-mask this ROI is located.
                % niidata is the data which will be written to the
                % newly generated nii-file
                idx=find(roi_intens==niidata);

                if ~isempty(idx)
                    niidata_temp(idx)=lagmap_avrg(roi);
                    % set the corresponding values to -1 in niidata_binary 
                    niidata_binary(idx)=-1;
                else
                    commandwindow
                    error('idx empty','error','error');
                end
            end
            niidata=niidata_temp;
            clear niidata_temp

            % Now, if niidata_binary+niidata_binary_ref is zero everywhere, only the
            % entries where niidata was different from zero were changed.

            if length(find(niidata_binary+niidata_binary_ref~=0))
                commandwindow
                error('difference larger than zero','error','error');
            end
            clear niidata_binary niidata_binary_ref

            % Create a nifti-file with (not so) random name in the temp-folder, then smooth using spm, delete nii-file, then output
            NewfileName=strcat(filenames{file},MapSeedMaskname(1:end-4),'-lagmap')

            [niidata_rows,niidata_columns,niidata_pages]=size(niidata);
            if niidata_rows~=91 || niidata_columns~=109 || niidata_pages~=91
                commandwindow
                error('niidata is of compromising matrix size. Only implemented for 2mmx2mmx2mm nifti files.');
            else
                write_nifti(m.mask1,strcat(StrgDIRECTORIES{MaskIter,1},filesep,NewfileName,'.nii'),niidata,[3 91 109 91 0 0 0 0]);
            end

            disp(['Number of nan-entries: ' num2str(length(find(isnan(lagmap_avrg)))) ' (of ' num2str(length(lagmap_avrg)) ')']);

            sx=str2double(get(handles.sx,'string'))
            sy=str2double(get(handles.sy,'string'))
            sz=str2double(get(handles.sz,'string'))

            spm_smooth(strcat(StrgDIRECTORIES{MaskIter,1},filesep,NewfileName,'.nii,1'),strcat(StrgDIRECTORIES{MaskIter,1},filesep,['S' num2str(sx) num2str(sy) num2str(sz) ''],NewfileName,'.nii,1'),[str2double(get(handles.sx,'string')) str2double(get(handles.sy,'string')) str2double(get(handles.sz,'string'))],0);
            disp('Smoothing done.');

            if file==1
               files{1}=strcat(NewfileName,'.nii');
               files_smooth{1}=strcat(['S' num2str(sx) num2str(sy) num2str(sz) ''],NewfileName,'.nii');
            else
                files=[files strcat(NewfileName,'.nii')];
                files_smooth=[files_smooth, strcat(['S' num2str(sx) num2str(sy) num2str(sz) ''],NewfileName,'.nii')];
            end

        end

        if get(handles.groupnif,'Value');
            groupnifti(files,strcat(StrgDIRECTORIES{MaskIter,1},filesep));
            groupnifti(files_smooth,strcat(StrgDIRECTORIES{MaskIter,1},filesep));
        end
        clear files files_smooth

    end
    clear MaskIter
% Smooth LagThread (Thread Cascade)
elseif get(handles.smooth_lagthread,'Value')            
    if iscell(filename)
        filenames=filename;
    else
       filenames=cell(1,1);
       filenames{1,1}=filename;
    end
    commandwindow
    for file=1:length(filenames)
        disp(['' num2str(length(filenames)-file+1) ' file(s) to go for this seed.']);
        disp(['Next file: ' filenames{file} '']);
        m=load(strcat(filepath,filenames{file}));

        if ~isfield(m,'L')
           commandwindow
           error('No field ''L'' in file. Perform PCA first to create Lag-Threads (i.e. Thread Cascades).');
        end

        %Load corresponding brainmask in MNI-space:
        if ~strcmp(m.mask1,m.mask2)
            commandwindow
            error('Using two different brainmasks','error','error');
        else
            [niidata,~,~,~,~] = readnifti(m.mask1);

            % niidata_binary is a matrix with 1's where niidata has a
            % value different from 0, otherwise it's 0.
            niidata_binary=niidata;
            niidata_binary_ref=niidata;

            niidata_binary(find(niidata~=0))=1;
            niidata_binary(find(niidata==0))=0;

            niidata_binary_ref(find(niidata~=0))=1;
            niidata_binary_ref(find(niidata==0))=0;

        end

        % rois_idx contains the columns in TD which correspond to the
        % lag thread

        % TD_rois_cols contains the list of rois along the second
        % dimension of TD (as many ROIs as columns)
        rois_idx=m.TDrois_cols;


        % Walk through every entry of the lag thread. The ith
        % entry is the value which has to be assigned to all voxels in
        % that ROI, the ROI being m.TDrois_cols(i)
        niidata_temp=niidata;
        for roi=1:length(m.TDrois_cols)
            roi_intens=m.TDrois_cols(roi);

            %Remember: niidata contains the niftidata matrix from
            %m.mask1. Find all entries in the mask-matrix that
            %correspond to the ROI of the lagprojection we are
            %currently iterating over
            idx=find(roi_intens==niidata);

            if ~isempty(idx)
                if ~get(handles.NegThread,'Value')
                    niidata_temp(idx)=m.L(roi,str2double(get(handles.LagThrNo,'string')));
                elseif get(handles.NegThread,'Value')
                    niidata_temp(idx)=-m.L(roi,str2double(get(handles.LagThrNo,'string')));
                end
                niidata_binary(idx)=-1;
            else
                commandwindow
                error('idx empty','error','error');
            end
        end
        niidata=niidata_temp;
        clear niidata_temp

        %Now, if niidata_binary+niidata_binary_ref is zero everywhere, only the
        %entries where niidata was different from zero were changed.

        if length(find(niidata_binary+niidata_binary_ref~=0))
            commandwindow
            error('Difference larger than zero','error','error');
        end
        clear niidata_binary niidata_binary_ref

        current_filename=filenames{file};
        name=strcat(current_filename(1:end-4),'-ThreadCascadeNo',get(handles.LagThrNo,'string'));

        [niidata_rows,niidata_columns,niidata_pages]=size(niidata);

        if niidata_rows~=91 || niidata_columns~=109 || niidata_pages~=91
            commandwindow
            error('niidata is of compromising matrix size. Only implemented for 2mmx2mmx2mm nifti files.');
        else
            write_nifti(m.mask1,strcat(DIRECTORYNAME,filesep,name,'.nii'),niidata,[3 91 109 91 0 0 0 0]);
        end

        disp(['Number of nan-entries: ' num2str(length(find(isnan(m.ColumnMeans)))) ' (of ' num2str(length(m.ColumnMeans)) ')']);

        sx=str2double(get(handles.sx,'string'))
        sy=str2double(get(handles.sy,'string'))
        sz=str2double(get(handles.sz,'string'))

        %Smooth the created nifti-file while creating a new,
        %smoothed file
        spm_smooth(strcat(DIRECTORYNAME,filesep,name,'.nii,1'),strcat(DIRECTORYNAME,filesep,['S' num2str(sx) num2str(sy) num2str(sz) ''],name,'.nii,1'),[str2double(get(handles.sx,'string')) str2double(get(handles.sy,'string')) str2double(get(handles.sz,'string'))],0);
        disp('Smoothing done.');

        if file==1
           files{1}=strcat(name,'.nii');
           files_smooth{1}=strcat(['S' num2str(sx) num2str(sy) num2str(sz) ''],name,'.nii');
        else
            files=[files strcat(name,'.nii')];
            files_smooth=[files_smooth, strcat(['S' num2str(sx) num2str(sy) num2str(sz) ''],name,'.nii')];
        end
    end
end
disp('nii-file generation terminated.');
set(handles.smooth,'Enable','On');

if ~get(handles.smooth_lagmap,'Value')
    if get(handles.groupnif,'Value');
        groupnifti(files,strcat(DIRECTORYNAME,filesep));
        groupnifti(files_smooth,strcat(DIRECTORYNAME,filesep));
    end
end

function sx_Callback(hObject, eventdata, handles)
% hObject    handle to sx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sx as text
%        str2double(get(hObject,'String')) returns contents of sx as a double

% --- Executes during object creation, after setting all properties.
function sx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sy_Callback(hObject, eventdata, handles)
% hObject    handle to sy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sy as text
%        str2double(get(hObject,'String')) returns contents of sy as a double


% --- Executes during object creation, after setting all properties.
function sy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function sz_Callback(hObject, eventdata, handles)
% hObject    handle to sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sz as text
%        str2double(get(hObject,'String')) returns contents of sz as a double


% --- Executes during object creation, after setting all properties.
function sz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in smooth_lagthread.
%function smooth_lagthread_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_lagthread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of smooth_lagthread


% --- Executes on button press in smooth_lagmap.
%function smooth_lagmap_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_lagmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of smooth_lagmap


% --- Executes on button press in smooth_projection.
%function smooth_projection_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_projection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of smooth_projection


% --- Executes on button press in plotnii.
function plotnii_Callback(hObject, eventdata, handles)
% hObject    handle to plotnii (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

if isfield(handles,'RecentNiiPlotPath')
    [niiname,niipath]=uigetfile(['' handles.RecentNiiPlotPath '*.nii'],'MultiSelect','On','Pick a nii-file to plot');
else
    [niiname,niipath]=uigetfile(['' pwd filesep '*.nii'],'MultiSelect','On','Pick a nii-file to plot');
end
if ~isequal(niipath,0)
    handles=guidata(hObject);
    handles.RecentNiiPlotPath=niipath;
    guidata(hObject,handles);
else
    return
end

    ThisFile=mfilename('fullpath');
    path=strcat(ThisFile(1:end-length(mfilename)),'masks',filesep);
    [maskname,maskpath]=uigetfile(['' path '*.nii'],'MultiSelect','On','Pick a mask-file');
    clear ThisFile path

    islice=str2double(get(handles.i,'string'));
    jslice=str2double(get(handles.j,'string'));
    kslice=str2double(get(handles.k,'string'));

    ifind=strfind(get(handles.i,'string'),':');
    jfind=strfind(get(handles.j,'string'),':');
    kfind=strfind(get(handles.k,'string'),':');
    
    % Only if one view (sagittal, coronal or transveral) is plotted!
    % Otherwise the option to plot several slices at once is not possible!
    if ~isempty(ifind)
        itemp=get(handles.i,'string');
        if length(ifind)==1
            islice=str2double(itemp(1:ifind(1)-1)):str2double(itemp(ifind(1)+1:end));
        elseif length(ifind)==2
            islice=str2double(itemp(1:ifind(1)-1)):str2double(itemp(ifind(1)+1:ifind(2)-1)):str2double(itemp(ifind(2)+1:end));
        end
    end
    if ~isempty(jfind)
        jtemp=get(handles.j,'string');
        if length(jfind)==1
            jslice=str2double(itemp(1:jfind(1)-1)):str2double(jtemp(jfind(1)+1:end));
        elseif length(jfind)==2
            jslice=str2double(jtemp(1:jfind(1)-1)):str2double(jtemp(jfind(1)+1:jfind(2)-1)):str2double(jtemp(jfind(2)+1:end));
        end
    end
    if ~isempty(kfind)
        ktemp=get(handles.k,'string');
        if length(kfind)==1
            kslice=str2double(ktemp(1:kfind(1)-1)):str2double(ktemp(kfind(1)+1:end));
        elseif length(kfind)==2
            kslice=str2double(ktemp(1:kfind(1)-1)):str2double(ktemp(kfind(1)+1:kfind(2)-1)):str2double(ktemp(kfind(2)+1:end));
        end
    end

    ThisFile=mfilename('fullpath');
    brainmask_backdrop=strcat(ThisFile(1:end-length(mfilename)),'masks',filesep,'fsl_MNI152_T1_2mm.nii');
    clear ThisFile
    
    if str2double(get(handles.c_range_from,'String'))==0 && str2double(get(handles.c_range_to,'String'))==0 
        error('Not a valid Color Range Option.')
    end
    
    for islice_iter=islice
        for jslice_iter=jslice
            for kslice_iter=kslice
                create_nii_plots(strcat(niipath,niiname),strcat(maskpath,maskname),{islice_iter,jslice_iter,kslice_iter},{get(handles.plot_i,'Value'),get(handles.plot_j,'Value'),get(handles.plot_k,'Value')},get(handles.slicepos,'Value'),brainmask_backdrop,{str2double(get(handles.c_range_from,'String')),str2double(get(handles.c_range_to,'String')),get(handles.set_c_range,'Value')},get(handles.revert_plot,'Value'),{get(handles.SavePlotsPNG,'Value'),get(handles.PNGDIRTXT,'string')});
            end
        end
    end


% --- Executes on button press in groupnifti.
function groupnifti_Callback(hObject, eventdata, handles)
% hObject    handle to groupnifti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc; 
set(handles.groupnifti,'Enable','Off');
drawnow;

[filename,filepath]=uigetfile(['' pwd filesep '*.nii'],'MultiSelect','On','Select files over which to average.');
if ~isequal(filepath,0)
    set(handles.groupnifti,'Enable','On');
    groupnifti(filename,filepath);
else
    set(handles.groupnifti,'Enable','On');
    return
end


% --- Executes on button press in smooth_lagproj.
function smooth_lagproj_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_lagproj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of smooth_lagproj
set(handles.groupnif,'Visible','On');
set(handles.NegThread,'Visible','Off');
set(handles.LagThrNo,'Visible','Off');
set(handles.text36,'Visible','Off');


% --- Executes during object creation, after setting all properties.
function groupnifti_CreateFcn(hObject, eventdata, handles)
% hObject    handle to groupnifti (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox26.
%function checkbox26_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox26


% --- Executes during object creation, after setting all properties.
function plotlagmap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotlagmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in groupnif.
%function groupnif_Callback(hObject, eventdata, handles)
% hObject    handle to groupnif (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of groupnif


% --- Executes on button press in SaveCorrFig.
function SaveCorrFig_Callback(hObject, eventdata, handles)
% hObject    handle to SaveCorrFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SaveCorrFig


% --- Executes on button press in SavePlotsPNG.
function SavePlotsPNG_Callback(hObject, eventdata, handles)
% hObject    handle to SavePlotsPNG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SavePlotsPNG

if get(hObject,'Value')
    set(handles.PNGDirectory,'Visible','On');
    if ~isfield(handles,'PNGDIRECTORY')
        DIRECTORYNAME=uigetdir(pwd,'Select Storage Directory for PNG files.');
            if ~isequal(DIRECTORYNAME,0) % Make sure user didn't hit 'cancel'
             	handles=guidata(hObject);
              	handles.PNGDIRECTORY=DIRECTORYNAME;
                guidata(hObject,handles);
                set(handles.PNGDIRTXT,'string',DIRECTORYNAME);
            end

        
    end
    
else
    set(handles.PNGDirectory,'Visible','Off');
end


% --- Executes on button press in ChngPNGDIR.
function ChngPNGDIR_Callback(hObject, eventdata, handles)
% hObject    handle to ChngPNGDIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'PNGDIRECTORY')
                        DIRECTORYNAME=uigetdir(['' handles.PNGDIRECTORY ''],'Select Storage Directoryf or PNG files.');
else
                        DIRECTORYNAME=uigetdir(pwd,'Select Storage Directory for PNG files.');
end

if ~isequal(DIRECTORYNAME,0) % Make sure user didn't hit 'cancel'
                            handles=guidata(hObject);
                            handles.PNGDIRECTORY=DIRECTORYNAME;
                            guidata(hObject,handles);
                            set(handles.PNGDIRTXT,'string',DIRECTORYNAME);
end


% --- Executes on button press in allowmasks.
function allowmasks_Callback(hObject, eventdata, handles)
% hObject    handle to allowmasks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allowmasks


% --- Executes on button press in groupleveltd.
function pushbutton51_Callback(hObject, eventdata, handles)
% hObject    handle to groupleveltd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in storeredund2.
function storeredund2_Callback(hObject, eventdata, handles)
% hObject    handle to storeredund2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of storeredund2


function choosefilename_Callback(hObject, eventdata, handles)
% hObject    handle to choosefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of choosefilename as text
%        str2double(get(hObject,'String')) returns contents of choosefilename as a double


% --- Executes on button press in smooth_lagthread.
function smooth_lagthread_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_lagthread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of smooth_lagthread

set(handles.NegThread,'Visible','On');
set(handles.LagThrNo,'Visible','On');
set(handles.text36,'Visible','On');


% --- Executes on button press in smooth_lagmap.
function smooth_lagmap_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_lagmap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of smooth_lagmap

set(handles.groupnif,'Visible','On');
set(handles.NegThread,'Visible','Off');
set(handles.LagThrNo,'Visible','Off');
set(handles.text36,'Visible','Off');


% --- Executes on button press in NegThread.
function NegThread_Callback(hObject, eventdata, handles)
% hObject    handle to NegThread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of NegThread


function LagThrNo_Callback(hObject, eventdata, handles)
% hObject    handle to LagThrNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LagThrNo as text
%        str2double(get(hObject,'String')) returns contents of LagThrNo as a double

% --- Executes during object creation, after setting all properties.
function LagThrNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LagThrNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LagThreadsPermutationTest.
function LagThreadsPermutationTest_Callback(hObject, eventdata, handles)
% hObject    handle to LagThreadsPermutationTest (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This function performs a permutation test on all combinations of Lag-Threads from group 1 (Threads 1:LTGRP1max) and group 2 (Threads 1:LTGRP2max) and
% checks if they differ significantly between groups.

if isfield(handles,'RecentPermTestFilespathGroup1')
    [group1FILENAME, group1PATHNAME] = uigetfile(['' handles.RecentPermTestFilespathGroup1 '*.mat'],'Choose .mat-groupfile of group 1.');
else
    [group1FILENAME, group1PATHNAME] = uigetfile(['' pwd filesep '*.mat'],'Choose .mat-groupfile of group 1.');
end

if ~isequal(group1PATHNAME,0)
    handles=guidata(hObject);
    handles.RecentPermTestFilespathGroup1=group1PATHNAME;
    guidata(hObject,handles);
else
    return
end

if isfield(handles,'RecentPermTestFilespathGroup2')
    [group2FILENAME, group2PATHNAME] = uigetfile(['' handles.RecentPermTestFilespathGroup2 '*.mat'],'Choose .mat-groupfile of group 2.');
else
    [group2FILENAME, group2PATHNAME] = uigetfile(['' pwd filesep '*.mat'],'Choose .mat-groupfile of group 2.');
end

if ~isequal(group2PATHNAME,0)
    handles=guidata(hObject);
    handles.RecentPermTestFilespathGroup2=group2PATHNAME;
    guidata(hObject,handles);
else
    return
end

GRPfile1=strcat(group1PATHNAME,group1FILENAME);
groupfile1=load(GRPfile1);

GRPfile2=strcat(group2PATHNAME,group2FILENAME);
groupfile2=load(GRPfile2);

if get(handles.PermutationTestThreadCascades,'Value')
    commandwindow;
    disp('Performing Permutation Test on Lag-Threads.');
    fprintf('\n...\n');
    
    % NoSubjects = number of subjects
    NoSubjects=size(groupfile1.filenames,1)+size(groupfile2.filenames,1);
    NoSbjctsGrp1=size(groupfile1.filenames,1);
    NoSbjctsGrp2=size(groupfile2.filenames,1);

    if ~isfield(groupfile1,'L') || ~isfield(groupfile1,'L')
        commandwindow;
        error('No Lag-Thread matrix L in one or both of the selected files. Perform PCA on group files first.');
    end

    NoThreads=size(groupfile1.L,2);
    NoRegions=size(groupfile1.L,1);

    
    %%% USER INPUT %%%%
    version=2;
    % version=1: Use PCs, version=2: use actual lag threads

    iterations=str2double(get(handles.PermutationTestIterations,'string'));

    LTGRP1max=str2double(get(handles.MaxThread,'string'));
    LTGRP2max=str2double(get(handles.MaxThread,'string'));;
    %%%%%%%%%%%%%%%%%%%
    
    TDs=zeros(NoThreads,NoThreads,NoSubjects);
    subjects=cell(1,NoSubjects);
    subjects(1:NoSbjctsGrp1)=groupfile1.filenames;
    subjects(NoSbjctsGrp1+1:NoSubjects)=groupfile2.filenames;

    group_1_actual=1:NoSbjctsGrp1;
    group_2_actual=NoSbjctsGrp1+1:NoSubjects;

    for grp1_subj=1:NoSbjctsGrp1
        subjectfile=groupfile1.filenames{grp1_subj};
        m=load(strcat(groupfile1.pathtofiles,subjectfile));

        TDs(:,:,grp1_subj)=m.TD;
    end

    for grp2_subj=1:NoSbjctsGrp2
        subjectfile=groupfile2.filenames{grp2_subj};
        m=load(strcat(groupfile2.pathtofiles,subjectfile));

        TDs(:,:,grp1_subj+grp2_subj)=m.TD;
    end

    groupTD_1_actual=nanmean(TDs(:,:,group_1_actual),3);
    groupTD_2_actual=nanmean(TDs(:,:,group_2_actual),3);

    groupTDz_1_actual=bsxfun(@minus,groupTD_1_actual,nanmean(groupTD_1_actual));
    groupTDz_2_actual=bsxfun(@minus,groupTD_2_actual,nanmean(groupTD_2_actual));

    [COMPONENTS_1_actual,~,~]=pcacov(cov(groupTDz_1_actual));
    [COMPONENTS_2_actual,~,~]=pcacov(cov(groupTDz_2_actual));

    L1_actual=(-1)*(groupTDz_1_actual*COMPONENTS_1_actual)/sqrt(size(groupTDz_1_actual,1));
    L2_actual=(-1)*(groupTDz_2_actual*COMPONENTS_2_actual)/sqrt(size(groupTDz_2_actual,1));
        
    p_coll=[];
    for grp1_comp=1:LTGRP1max
        for grp2_comp=1:LTGRP2max
        
        if version==1 %PCs
            actual_difference=COMPONENTS_1_actual(:,grp1_comp)-COMPONENTS_2_actual(:,grp2_comp);
        elseif version==2 %Lag Threads
            actual_difference=L1_actual(:,grp1_comp)-L2_actual(:,grp2_comp);
        end

        difference_rand=[];

        for iteration=1:iterations

            subject_order_rand=randperm(NoSubjects);
            group_1_rand=subject_order_rand(1:NoSbjctsGrp1);
            group_2_rand=subject_order_rand(NoSbjctsGrp1+1:NoSubjects);

            % Compute groupTD of both random groups
            groupTD_1=nanmean(TDs(:,:,group_1_rand),3);
            groupTD_2=nanmean(TDs(:,:,group_2_rand),3);

            groupTDz_1=bsxfun(@minus,groupTD_1,nanmean(groupTD_1));
            groupTDz_2=bsxfun(@minus,groupTD_2,nanmean(groupTD_2));

            [COMPONENTS_1,~,~]=pcacov(cov(groupTDz_1));
            [COMPONENTS_2,~,~]=pcacov(cov(groupTDz_2));

            L1_rand=(-1)*(groupTDz_1*COMPONENTS_1)/sqrt(size(groupTDz_1,1));
            L2_rand=(-1)*(groupTDz_2*COMPONENTS_2)/sqrt(size(groupTDz_2,1));

            % For every iteration, a new column is stored in
            % difference_rand where each column's row entry corresponds to the
            % group differences in the given region.
            if version==1
                difference_rand=[difference_rand, COMPONENTS_1(:,grp1_comp)-COMPONENTS_2(:,grp2_comp)];
            elseif version==2
                difference_rand=[difference_rand, L1_rand(:,grp1_comp)-L2_rand(:,grp2_comp)];
            end

        end

        % Test normality
        % kstest tests for a standard normal distribution
        % (h=0 -> normal distro)
        for i=1:NoRegions % The number of threads is assumed to be equal to the number of regions.
            kstest((difference_rand(i,:)-mean(difference_rand(i,:)))/std(difference_rand(i,:)));
        end

        % Now, vector p stores the p value for every region
        for entry=1:NoRegions
            p(entry)=length(find(abs(difference_rand(entry,:))>=abs(actual_difference(entry))))/iterations;
        end

        p_coll=[p_coll; grp1_comp, grp2_comp, p];
        end
    end

    assignin('base','p_coll',p_coll)
    p_coll;
    p_values=p_coll(:,3:end);
    
    ThreadGroup1=p_coll(:,1);
    ThreadGroup2=p_coll(:,2);
    VariableNames=cell(1,3);
    VariableNames{1}='ThreadGroup1';
    VariableNames{2}='ThreadGroup2';
    for reg=1:NoRegions
        ROI=p_coll(:,2+reg);
        VariableNames{3}=['ROI' num2str(groupfile1.TDrois_rows(reg)) ''];
        T=table(ThreadGroup1,ThreadGroup2,ROI,'VariableNames',VariableNames)
        fprintf(['p-values for ROI with value ' num2str(groupfile1.TDrois_rows(reg)) '.\n']);
        fprintf('------------------------------------------------\n')
    end

    StoreFile=questdlg('Store Results?','Store result of Analysis as mat-file??','Yes','No','Yes');

    switch StoreFile

        case 'Yes'

            DIRECTORYNAME=uigetdir('Where to store permutation results?');

            if version==1
                description=['used principal components, ' num2str(iterations) ' iterations. First column describes used PC of group 1, second column used PC of group 2, the following columns belong to the 17 mask regions.'];
                save(['' DIRECTORYNAME filesep 'permutation_test_PCs'],'p_coll','description');
            else
                description=['used lag threads, ' num2str(iterations) ' iterations. First column describes used lag thread of group 1, second column used lag thread of group 2, the following columns belong to the 17 mask regions.'];
                save(['' DIRECTORYNAME filesep 'permutation_test_lagthreads'],'ThreadGroup1','ThreadGroup1','p_values','description');
            end

    end
    
    clear p_coll
elseif get(handles.PermutationTestGlobalCascades,'Value')
    % Perform ttest on Global Casacdes. This is to be used when the
    % brainmask consists of discontinuous regions
    
    % NoSubjects = number of subjects
    NoSubjects=size(groupfile1.filenames,1)+size(groupfile2.filenames,1);
    NoSbjctsGrp1=size(groupfile1.filenames,1);
    NoSbjctsGrp2=size(groupfile2.filenames,1);

    NoColumns=size(groupfile1.TD,2);
    NoRows=size(groupfile1.TD,1);
    
    GlobalCascades=zeros(NoSubjects,NoColumns);
    subjects=cell(1,NoSubjects);
    subjects(1:NoSbjctsGrp1)=groupfile1.filenames;
    subjects(NoSbjctsGrp1+1:NoSubjects)=groupfile2.filenames;

    GroupAffiliation=zeros(NoSubjects,1);
    for grp1_subj=1:NoSbjctsGrp1
        subjectfile=groupfile1.filenames{grp1_subj};
        m=load(strcat(groupfile1.pathtofiles,subjectfile));
        GlobalCascades(grp1_subj,:)=m.ColumnMeans;
        GroupAffiliation(grp1_subj,1)=1;
        GlobalCascadesGrp1=GlobalCascades(grp1_subj,:);
    end

    for grp2_subj=1:NoSbjctsGrp2
        subjectfile=groupfile2.filenames{grp2_subj};
        m=load(strcat(groupfile2.pathtofiles,subjectfile));
        GlobalCascades(grp1_subj+grp2_subj,:)=m.ColumnMeans;
        GroupAffiliation(grp1_subj+grp2_subj,1)=2;
        GlobalCascadesGrp2=GlobalCascades(grp1_subj+grp2_subj,:);
    end
    
    pValue=zeros(NoColumns,1);
    normality=zeros(NoColumns,2);
    for Column=1:NoColumns
       % Test normality with kstest
       % kstest tests for a standard normal distribution
       % (h=0 -> normal distro)
       normality(Column,1)=~kstest((GlobalCascades(1:NoSbjctsGrp1,Column)-mean(GlobalCascades(1:NoSbjctsGrp1,Column)))/std(GlobalCascades(1:NoSbjctsGrp1,Column)));
       normality(Column,2)=~kstest((GlobalCascades(NoSbjctsGrp1+1:NoSubjects,Column)-mean(GlobalCascades(NoSbjctsGrp1+1:NoSubjects,Column)))/std(GlobalCascades(NoSbjctsGrp1+1:NoSubjects,Column)));
       if 0
           % Two-sample ttest
           [h,p]=ttest2(GlobalCascades(1:NoSbjctsGrp1,Column),GlobalCascades(NoSbjctsGrp1+1:NoSubjects,Column),'tail','both');
           PVALUES(Column)=p;
       end
       % ANOVA for covariate inclusion
       [p,~]=anovan(GlobalCascades(:,Column),{GroupAffiliation},'display','off');
       pValue(Column)=p;
    end
%     disp('p values determined by 2-sample ttest.');
    commandwindow
    Region=groupfile1.TDrois_cols;
%     table(Region,normality,pValue)
    table(Region,pValue)
    
elseif get(handles.PermutationTestSeededCascades,'Value')
    commandwindow;
    disp('Performing Statistical Significance Test on Lag Maps.');
    fprintf('\n...\n');
    
    
    % NoSubjects = number of subjects
    NoSubjects=size(groupfile1.filenames,1)+size(groupfile2.filenames,1);
    NoSbjctsGrp1=size(groupfile1.filenames,1);
    NoSbjctsGrp2=size(groupfile2.filenames,1);

    NoColumns=size(groupfile1.TD,2);
    NoRows=size(groupfile1.TD,1);
    % The Cascades will be stored as rows (onw row for one subjects)
    SeededCascades=zeros(NoSubjects,NoRows);
    subjects=cell(1,NoSubjects);
    subjects(1:NoSbjctsGrp1)=groupfile1.filenames;
    subjects(NoSbjctsGrp1+1:NoSubjects)=groupfile2.filenames;

    pValue=zeros(NoRows,1);

    GroupAffiliation=cell(NoSubjects,1);
    % Iterate through every Lag-Map. Lag-Maps are stored as ROWS. Each
    % subject one row
    COL=find(groupfile1.TDrois_cols==str2double(get(handles.PermutationTestIterations,'string')));
    for grp1_subj=1:NoSbjctsGrp1
        subjectfile=groupfile1.filenames{grp1_subj};
        m=load(strcat(groupfile1.pathtofiles,subjectfile));
        GroupAffiliation{grp1_subj,1}='p'; % premature
        SeededCascades(grp1_subj,:)=m.TD(:,COL)';
    end

    for grp2_subj=1:NoSbjctsGrp2
        subjectfile=groupfile2.filenames{grp2_subj};
        m=load(strcat(groupfile2.pathtofiles,subjectfile));
        GroupAffiliation{grp1_subj+grp2_subj,1}='m'; % mature
        SeededCascades(grp1_subj+grp2_subj,:)=m.TD(:,COL)';
    end

    GroupAffiliation=categorical(GroupAffiliation);
    % Now for the Lag-Map, look at every region
    for REGION=1:NoRows
        [p,~]=anovan(SeededCascades(:,REGION),{GroupAffiliation},'varnames',{'Group'},'display','off');
        pValue(REGION)=p; % Group
    end
    clear GroupAffiliation
    Region=groupfile1.TDrois_cols;
    table(Region,pValue)
end

disp('Done.');


% --- Executes on button press in PermutationTestLagThreads.
function PermutationTestThreadCascades_Callback(hObject, eventdata, handles)
% hObject    handle to PermutationTestLagThreads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PermutationTestLagThreads

if get(handles.PermutationTestThreadCascades,'Value')
    set(handles.MaxThread,'Visible','on');
    set(handles.text40,'Visible','on');
    set(handles.text39,'string','# iterations');
    set(handles.text39,'Visible','on');
    set(handles.PermutationTestIterations,'Visible','on');
    set(handles.PermutationTestIterations,'string','10000');
end

% --- Executes on button press in PermutationTestSeededCascades.
function PermutationTestSeededCascades_Callback(hObject, eventdata, handles)
% hObject    handle to PermutationTestSeededCascades (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PermutationTestSeededCascades

if get(handles.PermutationTestSeededCascades,'Value')
    set(handles.MaxThread,'Visible','off');
    set(handles.text40,'Visible','off');
    set(handles.text39,'string','Lag Map Seed');
    set(handles.text39,'Visible','on');
    set(handles.PermutationTestIterations,'Visible','on');
    set(handles.PermutationTestIterations,'string','Choose Seed');
end


% --- Executes on button press in PermutationTestGlobalCascades.
function PermutationTestGlobalCascades_Callback(hObject, eventdata, handles)
% hObject    handle to PermutationTestGlobalCascades (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PermutationTestGlobalCascades

if get(handles.PermutationTestGlobalCascades,'Value')
    set(handles.MaxThread,'Visible','off');
    set(handles.text40,'Visible','off');
    set(handles.text39,'Visible','off');
    set(handles.PermutationTestIterations,'Visible','off');
end


function MaxThread_Callback(hObject, eventdata, handles)
% hObject    handle to MaxThread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxThread as text
%        str2double(get(hObject,'String')) returns contents of MaxThread as a double


% --- Executes during object creation, after setting all properties.
function MaxThread_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxThread (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PermutationTestIterations_Callback(hObject, eventdata, handles)
% hObject    handle to PermutationTestIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PermutationTestIterations as text
%        str2double(get(hObject,'String')) returns contents of PermutationTestIterations as a double


% --- Executes during object creation, after setting all properties.
function PermutationTestIterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PermutationTestIterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
