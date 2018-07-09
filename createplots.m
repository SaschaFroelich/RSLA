function [] = createplots(varargin)
% "rois" are the corresponding intensity values of the ROIs, "lagmap" the
% corresponding lag-vector. these two vectors have to be of the same length
% plotting: 0=lagthread 1=lagmaps 2=lag projection

ThisFile=mfilename('fullpath');
path=strcat(ThisFile(1:end-length(mfilename)),'masks',filesep);
brainmask_backdrop=strcat(path,'fsl_MNI152_T1_2mm.nii');
clear ThisFile path

h0bject=varargin{1};
eventdata=varargin{2};
file=varargin{3};
filename=file(end-4:end);

rois=varargin{7};
lagmap=varargin{8};
plotting=varargin{9};
seedthread=varargin{10};
slicepos=varargin{11};

CRANGE=varargin{12};
SetRange=CRANGE{1};
RangeFrom=CRANGE{2};
RangeTo=CRANGE{3};

SAVEPNG=varargin{end};
SavePic=SAVEPNG{1};
SaveDir=SAVEPNG{2};

if isequal(h0bject,'none') % Fct was not called from slider callback
    islice=varargin{4};
    jslice=varargin{5};
    kslice=varargin{6};
    
    %Make "rois" into a row vector
    [rois1, rois2]=size(rois);
    if rois1~=1
       rois=rois'; 
    end

    % Make "lagmap" a row vector
    % THE ORDER IN LAGMAP DETERMINES THE COLOR MAP AND HAS TO BE IN THE SAME
    % ORDER AS VECTOR ROIS
    [lagmap1, lagmap2]=size(lagmap);
    if lagmap1~=1
       lagmap=lagmap'; 
    end
    
else
    if varargin{4}>=0
        islice=round(h0bject.Value);
        jslice=varargin{5};
        kslice=varargin{6};
    elseif varargin{5}>=0        
        islice=varargin{4};            %imgdata is the data of the maskfile
        jslice=round(h0bject.Value);
        kslice=varargin{6};
    elseif varargin{6}>=0
        islice=varargin{4};
        jslice=varargin{5};
        kslice=round(h0bject.Value);
    end
    
    % Now clear all the axes of the parent of h0bject, which is the figure
    % window on which the slider resides. Also clear all text, titles, and
    % remove the colorbar
    parent=h0bject.Parent;
    delete(findall(findall(gcf,'Type','axe'),'Type','text'));
    arrayfun(@cla,findall(parent,'Type','axes')); 
    arrayfun(@delete,findall(parent,'Type','Colorbar'));
    
end

% Determine how to name the figure windows
mtfile=varargin{13};
m=load(mtfile);
if plotting ~=2
    SeedRoiIntens=m.TDrois_rows(seedthread);
end
fig_names=strsplit(mtfile,filesep);
fig_names=strcat(fig_names{end-2},filesep,fig_names{end-1},filesep,fig_names{end});

tic
%R->L: x+ (91)
%P->A: y+ (109)
%I->S: z+ (91)

[data,~,~,~,~]=readnifti(file);
[data_backdrop,~,~,~,~]=readnifti(brainmask_backdrop);
header=readniftifileheader(file);

%There are several methods to display nifti-data visually. If sform_code>0,
%then "method 3" can be used (the one used here). See: https://brainder.org
%-> "The NIFTI-2 file format", or "The NIFTI file format"
if header.sform_code>0
    %Determine the spatial positions of the entries in "data"
    transf(1,:) = header.srow_x';
    transf(2,:) = header.srow_y';
    transf(3,:) = header.srow_z';
    
    transf=cat(1,header.srow_x',header.srow_y',header.srow_z');
    
    [rows,columns,pages]=size(data);
    
    %now, X,Y,Z symbolize the i,j,k for each voxel
    [X,Y,Z]=meshgrid(0:columns-1,0:rows-1,0:pages-1);
    
    
    for row=1:rows
        for col=1:columns
            for page=1:pages
                %Matlab starts at 1, but the voxel coord start at 0
                i=row-1;%=X(row,col,page)
                j=col-1;%=Y(row,col,page)
                k=page-1;%=Z(row,col,page)

                xyz=transf*[i;j;k;1];

                %R-L
                X(row,col,page)=xyz(1);
                %P-A
                Y(row,col,page)=xyz(2);
                %I-S
                Z(row,col,page)=xyz(3);

                %Now, X,Y,Z define a meshgrid of the actual positions in
                %space of the data points in data;
                
            end
        end
    end 
    
    %lagmapsanal is the first figure, so nfh will be 1, even if no actual
    %figure is opened.
    fh=findobj(0,'type','figure');
    nfh=length(fh);
    fig1=nfh+1;
    fig2=fig1+1;
    fig3=fig2+1;
    
    %fig 1 (riding along i = left - right)
    if islice >=0 
        ipage=islice+1;
    else
        ipage=islice-1;
    end
    %fig 2 (riding along j = anterior - posterior)
    if jslice>=0
        jpage=jslice+1;
    else
        jpage=jslice-1;
    end
    %fig 3 (riding along k = inferior - superior)
    if kslice>=0
        kpage=kslice+1;
    else
        kpage=kslice-1;
    end
    
    if ipage>0
        %-----FIGURE 1 (axes 1,4,7,8)----------------------%
        if ~isequal(h0bject,'none') % If createplots.m was called as callback fct from sliders
            f1=figure(parent.Number);
            ax1=axes('parent',f1);
            set(f1,'name',['' fig_names '']);
        else
            f1=figure(fig1);
            ax1=axes('parent',f1);
            movegui(f1,'northwest');
            set(f1,'name',['' fig_names ', sagittal view']);
        end
        
        % 109 in P->A Richtung
        for page=1:pages
            %ydata,zdata,imgdata etc are 2D meshgrids, while X,Y,Z, and
            %data re 3D matrices
            ydata(page,:)=Y(ipage,:,page);
            zdata(page,:)=Z(ipage,:,page);
            %imgdata is the data of the maskfile
            imgdata(page,:)=data(ipage,:,page);
            plot_imgdata(page,:)=data_backdrop(ipage,:,page);
        end
        %Plot gray background
        surf(ax1,ydata,zdata,plot_imgdata,'EdgeColor','None');
        colormap(ax1,'gray');
        xlabel(ax1,'I');
        ylabel(ax1,'P');
        %axis(ax1,[-126 90 -72 108]);
        %axis equal tight
        if plotting==1
            if slicepos
                title(ax1,['Lag map with seed ' num2str(SeedRoiIntens) '. (i=' num2str(abs(islice)) ',j=' num2str(abs(jslice)) ',k=' num2str(abs(kslice)) ')']);
            else
                title(ax1,['Lag map with seed ' num2str(SeedRoiIntens) '']);
            end
        elseif plotting==0
            if slicepos
                title(ax1,['Lag Thread No. ' num2str(seedthread) '. (i=' num2str(abs(islice)) ',j=' num2str(abs(jslice)) ',k=' num2str(abs(kslice)) ')']);
            else
                title(ax1,['Lag Thread No. ' num2str(seedthread) '']);
            end
        elseif plotting==2 %lag projection
            if slicepos
                title(ax1,['Lag projection (i=' num2str(abs(islice)) ',j=' num2str(abs(jslice)) ',k=' num2str(abs(kslice)) ')']);
            else
                title(ax1,'Lag Projection');
            end
        end
        view(ax1,2);

        %Find the points in "imgdata" (=mask) that correspond to the nodes (=rois)
        %I1 is a "roi-to-voxel-map": If the first ROI in rois has positions
        %113,114,200,201 in imgdata (the mask-file), then I1 starts off
        %with I1=[113,114,200,201,...]
        I1=[];
        %create colormap
        color1=[];
        for iter=1:length(rois)
            %Find all voxels in imgdata that belong to that specific ROI
            I=find(imgdata==rois(iter))';
            I1=[I1, I];
            color=ones(1,length(I))*lagmap(iter);
            %For every plotted point (=voxel) there has to be a color entry
            %in the vector "color1".
            color1=[color1, color];
        end

        [rows,columns]=size(imgdata);

        %The point with index I is in the following column in imgdata.
        %Set col=[] and row =[] in case I1 is empty, in order to delete the old
        %values;
        col=[];
        row=[];
        for iter=1:length(I1)
            col(iter)=ceil(I1(iter)/rows);
            row(iter)=I1(iter)-(col(iter)-1)*rows;
        end
        
        
        %Take any column in zdata, as all columns are the same
        %perhaps put +1 for better visual representation
        row_I1=zdata(row,1)+1; %+1
        %Take any row in ydata, as all rows are the same
        col_I1=ydata(1,col)+1; %+1

        ax4 = axes('parent',f1);
        %REMEMBER: surf works with meshgrid function, while scatter does
        %not (scatter only places single points)
        
        %scatter onto the SAME points as specified by the maskfile, only
        %with color as specified by color1, which is defined via the lagmap
        scatter(ax4,col_I1,row_I1,18,color1,'filled')
        %Link them together
        linkaxes([ax1,ax4])
        % Hide the top axes
        ax4.Visible = 'off';
        ax4.XTick = [];
        ax4.YTick = [];
        %Give each one its own colormap
        colormap(ax4,'jet')
        % Then add colorbars and get everything lined up
        %position of axes
        set([ax1,ax4],'Position',[.10 .11 .685 .815]);
        %position of colorbar
        cb2 = colorbar(ax4,'Position',[.81 .11 .0675 .815]);
        %colormap
        if ~SetRange
            caxis(ax4,[floor(min(lagmap)),ceil(max(lagmap))]);
        else
            caxis(ax4,[RangeFrom,RangeTo]);
        end
        axis([ax1,ax4],[-126 90 -72 108]);
        
        
        slider_i=uicontrol('Style','Slider','Min',0,'Max',90,'SliderStep',[1 1]/91,'Position',[270,6,250,16],'Value',islice);
        i_text=uicontrol('Style','Text','String',num2str(abs(islice)),'Position',[530,10,26,14]);
        set(slider_i,'Callback',{@createplots,file,slider_i.Value,-abs(jslice),-abs(kslice),rois,lagmap,plotting,seedthread,slicepos,CRANGE,mtfile,SAVEPNG});
        
        if slicepos
            %ax7 for positioning the transverse slice position
            ax7 = axes('parent',f1);
            %entries of Z only change along third index of Z
            %transv is position along I-S axis
            transv=Z(1,1,abs(kpage));
            scatter(ax7,linspace(min(min(ydata)),max(max(ydata)),2*length(ydata(1,:))),ones(1,218)*transv,1,'white','filled')
            %Link them together
            linkaxes([ax1,ax4,ax7])
            % Hide the top axes
            ax7.Visible = 'off';
            ax7.XTick = [];
            ax7.YTick = [];
            %position of axes
            set([ax1,ax4,ax7],'Position',[.10 .11 .685 .815]);
            axis([ax1,ax4,ax7],[-126 90 -72 108]);


            %ax8 for positioning the coronal (=frontal) slice position
            ax8 = axes('parent',f1);
            %entries of Y only change along second index of Y
            %coron is position along P-A axis
            coron=Y(1,abs(jpage),1);
            scatter(ax8,ones(1,218)*coron,linspace(min(min(zdata)),max(max(zdata)),2*length(zdata(1,:))),1,'white','filled')
            %Link them together
            linkaxes([ax1,ax4,ax7,ax8])
            % Hide the top axes
            ax8.Visible = 'off';
            ax8.XTick = [];
            ax8.YTick = [];
            %position of axes
            set([ax1,ax4,ax7,ax8],'Position',[.10 .11 .685 .815]);
            axis([ax1,ax4,ax7,ax8],[-126 90 -72 108]);
        end
        %---------------------------------------------%
        
        if SavePic
            print(gcf, '-dpng', ['' SaveDir '' filesep 'i' num2str(islice) '_j' num2str(jslice) '_k' num2str(kslice) '_sagittal.png']);
        end
        
    end

  
    if jpage>0
        %-----FIGURE 2 (axes 2,5,9,10)-----------------------------%
        if ~isequal(h0bject,'none')%If createplots.m was called as callback fct from sliders
            f2=figure(parent.Number);
            ax2=axes('parent',f2);
            set(f2,'name',['' fig_names '']);
        else
            f2=figure(fig2);
            ax2=axes('parent',f2);
            movegui(f2,'north');
            set(f2,'name',['' fig_names ', coronal view']);
        end
        
        for page=1:pages
            xdata2(:,page)=X(:,jpage,page);
            zdata2(:,page)=Z(:,jpage,page);
            imgdata2(:,page)=data(:,jpage,page);
            plot_imgdata2(:,page)=data_backdrop(:,jpage,page);
        end
        %orientate image like in fslview
        surf(ax2,xdata2,zdata2,plot_imgdata2,'EdgeColor','None');
        xlabel(ax2,'I');
        ylabel(ax2,'R');
        set(ax2,'Xdir','reverse');
        %text(24,66,'\leftarrow x=24,z=66');
        colormap(ax2,'gray');
        %axis equal
        if plotting==1
            if slicepos
                title(ax2,['Lag map with seed ' num2str(SeedRoiIntens) '. (i=' num2str(abs(islice)) ',j=' num2str(abs(jslice)) ',k=' num2str(abs(kslice)) ')']);
            else
                title(ax2,['Lag map with seed ' num2str(SeedRoiIntens) '']);
            end
        elseif plotting==0
            if slicepos
                title(ax2,['Lag Thread No. ' num2str(seedthread) '. (i=' num2str(abs(islice)) ',j=' num2str(abs(jslice)) ',k=' num2str(abs(kslice)) ')']);
            else
                title(ax2,['Lag Thread No. ' num2str(seedthread) '']);
            end
        elseif plotting==2 %lag projection
            if slicepos
                title(ax2,['Lag projection (i=' num2str(abs(islice)) ',j=' num2str(abs(jslice)) ',k=' num2str(abs(kslice)) ')']);
            else
                title(ax2,'Lag Projection');
            end
        end
        view(ax2,2);


        %Find the points in "imgdata" (=mask) that correspond to the nodes (=rois)
        I2=[];
        color2=[];
        %create colormap
        for iter=1:length(rois)
            I=find(imgdata2==rois(iter))';
            I2=[I2, I];
            color=ones(1,length(I))*lagmap(iter);
            color2=[color2, color];
        end


        [rows,columns]=size(imgdata2);

        %The point with index I is in the following column in imgdata
        %Set col=[] and row =[] in case I1 is empty, in order to delete the old
        %values;
        col=[];
        row=[];
        for iter=1:length(I2)
            col(iter)=ceil(I2(iter)/rows);
            row(iter)=I2(iter)-(col(iter)-1)*rows;
        end

        %Take any column in zdata, as all columns are the same
        %perhaps put +1 for better visual representation
        row_I2=xdata2(row,1)-1; %-1 -1,since x-axis reversed
        %Take any row in ydata, as all rows are the same
        col_I2=zdata2(1,col)+1; %+1

        ax5 = axes('parent',f2);
        scatter(ax5,row_I2,col_I2,18,color2,'filled')
        set(ax5,'Xdir','reverse');
        %Link them together
        linkaxes([ax2,ax5])
        % Hide the top axes
        ax5.Visible = 'off';
        ax5.XTick = [];
        ax5.YTick = [];
        %Give each one its own colormap
        colormap(ax5,'jet')
        % Then add colorbars and get everything lined up
        set([ax2,ax5],'Position',[.10 .11 .685 .815]);
        cb2 = colorbar(ax5,'Position',[.81 .11 .0675 .815]);
        if ~SetRange
            caxis(ax5,[floor(min(lagmap)),ceil(max(lagmap))]);
        else % manually set color range
            caxis(ax5,[RangeFrom,RangeTo]);
        end
        axis([ax2,ax5],[min(min(xdata2)) max(max(xdata2)) min(min(zdata2)) max(max(zdata2))]);
        %axis equal
        %text(-20,76,'\leftarrow y=-20,z=76');
        

        slider_j=uicontrol('Style','Slider','Min',0,'Max',108,'SliderStep',[1 1]/109,'Position',[270,6,250,16],'Value',jslice);
        j_text=uicontrol('Style','Text','String',num2str(abs(jslice)),'Position',[530,10,26,14]);
        set(slider_j,'Callback',{@createplots,file,-abs(islice),slider_j.Value,-abs(kslice),rois,lagmap,plotting,seedthread,slicepos,CRANGE,mtfile,SAVEPNG});
        
        

        if slicepos
            %ax9 for positioning the sagittal slice position
            ax9 = axes('parent',f2);
            %entries of X only change along first index of X
            %sagitt is position along R-L axis
            sagitt=X(abs(ipage),1,1);
            scatter(ax9,ones(1,2*length(zdata2(1,:)))*sagitt,linspace(min(min(zdata2)),max(max(zdata2)),2*length(zdata2(1,:))),1,'white','filled')
            set(ax9,'Xdir','reverse');
            %Link them together
            linkaxes([ax2,ax5,ax9])
            % Hide the top axes
            ax9.Visible = 'off';
            ax9.XTick = [];
            ax9.YTick = [];
            %position of axes
            set([ax2,ax5,ax9],'Position',[.10 .11 .685 .815]);
            axis([ax2,ax5,ax9],[min(min(xdata2)) max(max(xdata2)) min(min(zdata2)) max(max(zdata2))]);

            %ax10 for positioning the transverse slice position
            ax10 = axes('parent',f2);
            %entries of Z only change along third index of Z
            %transv is position along I-S axis
            transv=Z(1,1,abs(kpage));
            scatter(ax10,linspace(min(min(xdata2)),max(max(xdata2)),2*length(xdata2(1,:))),ones(1,2*length(xdata2(1,:)))*transv,1,'white','filled')
            %Link them together
            linkaxes([ax2,ax5,ax9,ax10])
            % Hide the top axes
            ax10.Visible = 'off';
            ax10.XTick = [];
            ax10.YTick = [];
            %position of axes
            set([ax2,ax5,ax9,ax10],'Position',[.10 .11 .685 .815]);
            axis([ax2,ax5,ax9,ax10],[min(min(xdata2)) max(max(xdata2)) min(min(zdata2)) max(max(zdata2))]);
        end

        %---------------------------------------------------------%
        if SavePic
            print(gcf, '-dpng', ['' SaveDir '' filesep 'i' num2str(islice) '_j' num2str(jslice) '_k' num2str(kslice) '_coronal.png']);
        end
        
    end

    if kpage>0
        %-----FIGURE 3 (axes 3,6,11,12)---------------------------%
        if ~isequal(h0bject,'none') %If createplots.m was called as callback fct from sliders
            f3=figure(parent.Number);
            ax3=axes('parent',f3);
            set(f3,'name',['' fig_names '']);
        else
            f3=figure(fig3);
            ax3=axes('parent',f3);
            movegui(f3,'northeast');
            set(f3,'name',['' fig_names ', transverse view']);
        end

        surf(ax3,X(:,:,kpage),Y(:,:,kpage),data_backdrop(:,:,kpage),'EdgeColor','None');
        xlabel(ax3,'P');
        ylabel(ax3,'R');
        set(ax3,'Xdir','reverse');
        colormap(ax3,'gray');
        if plotting==1
            if slicepos
                title(ax3,['Lag map with seed ' num2str(SeedRoiIntens) '. (i=' num2str(abs(islice)) ',j=' num2str(abs(jslice)) ',k=' num2str(abs(kslice)) ')']);
            else
                title(ax3,['Lag map with seed ' num2str(SeedRoiIntens) '']);
            end
        elseif plotting==0
            if slicepos
                title(ax3,['Lag Thread No. ' num2str(seedthread) '. (i=' num2str(abs(islice)) ',j=' num2str(abs(jslice)) ',k=' num2str(abs(kslice)) ')']);
            else
                title(ax3,['Lag Thread No. ' num2str(seedthread) '']);
            end
        elseif plotting==2 %lag projection
            if slicepos
                title(ax3,['Lag projection (i=' num2str(abs(islice)) ',j=' num2str(abs(jslice)) ',k=' num2str(abs(kslice)) ')']);
            else
                title(ax3,'Lag Projection');
            end
        end
        
        view(2);


        % Find the points in "imgdata" (=mask) that correspond to the nodes (=rois)
        I3=[];
        % Create colormap
        color3=[];
        for iter=1:length(rois)
            I=find(data(:,:,kpage)==rois(iter))';
            I3=[I3, I];
            color=ones(1,length(I))*lagmap(iter);
            color3=[color3, color];
        end


        [rows,columns]=size(data(:,:,kpage));

        % The point with index I is in the following column in imgdata
        % Set col=[] and row =[] in case I1 is empty, in order to delete the old
        % values;
        col=[];
        row=[];
        for iter=1:length(I3)
            col(iter)=ceil(I3(iter)/rows);
            row(iter)=I3(iter)-(col(iter)-1)*rows;
        end

        % Take any column in zdata, as all columns are the same
        % perhaps put +1 for better visual representation
        row_I3=X(row,1,kpage)-1; %-1, as x-axis reversed
        %Take any row in ydata, as all rows are the same
        col_I3=Y(1,col,kpage)+1; %+1

        ax6 = axes('parent',f3);
        scatter(ax6,row_I3,col_I3,18,color3,'filled')
        set(ax6,'Xdir','reverse');
        %Link them together
        linkaxes([ax3,ax6])
        % Hide the top axes
        ax6.Visible = 'off';
        ax6.XTick = [];
        ax6.YTick = [];
        %Give each one its own colormap
        colormap(ax6,'jet')
        % Then add colorbars and get everything lined up
        set([ax3,ax6],'Position',[.10 .11 .685 .815]);
        cb2 = colorbar(ax6,'Position',[.81 .11 .0675 .815]);
        if ~SetRange
            caxis(ax6,[floor(min(lagmap)),ceil(max(lagmap))]);
        else
            caxis(ax6,[RangeFrom,RangeTo]);
        end
        axis([ax3,ax6],[min(min(X(:,:,kpage))) max(max(X(:,:,kpage))) min(min(Y(:,:,kpage))) max(max(Y(:,:,kpage)))]);

        slider_k=uicontrol('Style','Slider','Min',0,'Max',90,'SliderStep',[1 1]/91,'Position',[270,6,250,16],'Value',kslice);
        k_text=uicontrol('Style','Text','String',num2str(abs(kslice)),'Position',[530,10,26,14]);
        set(slider_k,'Callback',{@createplots,file,-abs(islice),-abs(jslice),slider_k.Value,rois,lagmap,plotting,seedthread,slicepos,CRANGE,mtfile,SAVEPNG});
        
        if slicepos
            %ax11 for positioning the coronal (=frontal) slice position
            ax11 = axes('parent',f3);
            %entries of Y only change along second index of Y
            %coron is position along P-A axis
            coron=Y(1,abs(jpage),1);
            scatter(ax11,linspace(min(min(X(:,:,kpage))),max(max(X(:,:,kpage))),2*length(X(:,:,kpage))),ones(1,2*length(X(:,:,kpage)))*coron,1,'white','filled')
            %hold on
            %scatter(ax11,linspace(min(min(X(:,:,kpage))),max(max(X(:,:,kpage))),2*length(X(:,:,kpage))),ones(1,2*length(X(:,:,kpage)))*(coron+1),1,'white','filled')
            %Link them together
            linkaxes([ax3,ax6,ax11])
            % Hide the top axes
            ax11.Visible = 'off';
            ax11.XTick = [];
            ax11.YTick = [];
            %position of axes
            set([ax3,ax6,ax11],'Position',[.10 .11 .685 .815]);
            axis([ax3,ax6,ax11],[min(min(X(:,:,kpage))) max(max(X(:,:,kpage))) min(min(Y(:,:,kpage))) max(max(Y(:,:,kpage)))]);

            %ax12 for positioning the sagittal slice position
            ax12 = axes('parent',f3);
            %entries of X only change along first index of X
            %sagitt is position along R-L axis
            sagitt=X(abs(ipage),1,1);
            scatter(ax12,ones(1,2*length(Y(:,:,kpage)))*sagitt,linspace(min(min(Y(:,:,kpage))),max(max(Y(:,:,kpage))),2*length(Y(:,:,kpage))),1,'white','filled')
            set(ax12,'Xdir','reverse');
            %Link them together
            linkaxes([ax3,ax6,ax11,ax12])
            % Hide the top axes
            ax12.Visible = 'off';
            ax12.XTick = [];
            ax12.YTick = [];
            %position of axes
            set([ax3,ax6,ax11,ax12],'Position',[.10 .11 .685 .815]);
            axis([ax3,ax6,ax11,ax12],[min(min(X(:,:,kpage))) max(max(X(:,:,kpage))) min(min(Y(:,:,kpage))) max(max(Y(:,:,kpage)))]);
        end

    
        %---------------------------------------------%
        if SavePic
            print(gcf, '-dpng', ['' SaveDir '' filesep 'i' num2str(islice) '_j' num2str(jslice) '_k' num2str(kslice) '_tramsverse.png']);
        end
        
    end

    
        
    
end

end