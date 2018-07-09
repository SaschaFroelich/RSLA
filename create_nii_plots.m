function []=create_nii_plots(varargin)
%"rois" are the corresponding intensity values of the ROIs, "lagmap" the
%corresponding lag-vector. these two vectors have to be of the same length

% The last entry in vararginis 0 or 1 and determines whether or not to reverse plot

file=varargin{1};
mask=varargin{2};
SlicePos=varargin{3};
PlotTF=varargin{4};
islice=SlicePos{1};
jslice=SlicePos{2};
kslice=SlicePos{3};
ploti=PlotTF{1};
plotj=PlotTF{2};
plotk=PlotTF{3};
slicepos=varargin{5};
brainmask_backdrop=varargin{6};
ColorRange=varargin{7};
RevPlot=varargin{8};
PNGDir=varargin{9};

% Set Color Range if desired
if ColorRange{3}
    colorrange_min=ColorRange{1};
    colorrange_max=ColorRange{2};
end

if slicepos && (isnan(islice) || isnan(jslice) || isnan(kslice))
    msgbox('Slice-Position values of i,j or k not valid.','Errror','error');
    error('Slice-Position values of i,j or k not valid.');
end

[niftiplot,~,~,~,~]=readnifti(['' file '']);
if RevPlot==1
    % revert plot
   niftiplot=niftiplot*(-1); 
end

[niftiplotmask,~,~,~,~]=readnifti(['' mask '']);

niftiplot_masked=niftiplot.*niftiplotmask;

%Determine how to name the figure windows

fig_names=strsplit(file,filesep);
fig_names=strcat(fig_names{end-2},filesep,fig_names{end-1},filesep,fig_names{end});

tic
%R->L: x+ (91)
%P->A: y+ (109)
%I->S: z+ (91)

[data,~,~,~,~]=readnifti(file);
[data_backdrop,~,~,~,~]=readnifti(brainmask_backdrop);
header=readniftifileheader(file);

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
    
    %fig 1 (riding along i = left -> right)
        ipage=islice+1;
        
    %fig 2 (riding along j = anterior -> posterior)
        jpage=jslice+1;
        
    %fig 3 (riding along k = inferior -> superior)
        kpage=kslice+1;
    
    if ploti
        %-----FIGURE 1 (axes 1,4,7,8)----------------------%
            f1=figure(fig1);
            ax1=axes('parent',f1);
            movegui(f1,'northwest');
            set(f1,'name',['' fig_names ', sagittal view']);
        
        %109 in P->A Richtung
        for page=1:pages
            % ydata,zdata,imgdata etc are 2D meshgrids, while X,Y,Z, and
            % data re 3D matrices
            ydata(page,:)=Y(ipage,:,page);
            zdata(page,:)=Z(ipage,:,page);
            %imgdata is the data of the maskfile
            imgdata(page,:)=data(ipage,:,page);
            plot_imgdata(page,:)=data_backdrop(ipage,:,page);
            niftiplot_masked_imgdata(page,:)=niftiplot_masked(ipage,:,page);
        end
        % Plot gray background
        surf(ax1,ydata,zdata,plot_imgdata,'EdgeColor','None');
        
        colormap(ax1,'gray');
        xlabel(ax1,'I');
        ylabel(ax1,'P');
        %axis(ax1,[-126 90 -72 108]);
        %axis equal tight
        view(ax1,2);


        %Now use the data in niftiplot_masked_imgdata(page,:), assign the
        %correct "x-" and "y-" values and then scatter
        
        [rows,columns]=size(niftiplot_masked_imgdata);
        
        xvalues=[];
        yvalues=[];
        scatter_dat=[];
        for row=1:rows
            for col=1:columns
                if niftiplot_masked_imgdata(row,col)~=0
                    xvalues=[xvalues ydata(row,col)];
                    yvalues=[yvalues zdata(row,col)];
                    scatter_dat=[scatter_dat niftiplot_masked_imgdata(row,col)];
                end
            end
        end

        
        
        ax4 =axes('parent',f1);
        %REMEMBER: surf works with meshgrid function, while scatter does
        %not (scatter only places single points)
        
        %scatter onto the SAME points as specified by the maskfile, only
        %with color as specified by color1, which is defined via the lagmap
        scatter(ax4,xvalues,yvalues,18,scatter_dat,'filled')
        
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
        axis([ax1,ax4],[-126 90 -72 108]);
        if ~ColorRange{3}
            
        elseif ColorRange{3}
            caxis(ax4,[colorrange_min,colorrange_max]);
        end
       
        if slicepos
            % ax7 for positioning the transverse slice position
            ax7 =axes('parent',f1);
            %entries of Z only change along third index of Z
            %transv is position along I-S axis
            transv=Z(1,1,kpage);
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


            % ax8 for positioning the coronal (=frontal) slice position
            ax8 =axes('parent',f1);
            %entries of Y only change along second index of Y
            %coron is position along P-A axis
            coron=Y(1,jpage,1);
            scatter(ax8,ones(1,218)*coron,linspace(min(min(zdata)),max(max(zdata)),2*length(zdata(1,:))),1,'white','filled')
            linkaxes([ax1,ax4,ax7,ax8])
            % Hide the top axes
            ax8.Visible = 'off';
            ax8.XTick = [];
            ax8.YTick = [];
            % position of axes
            set([ax1,ax4,ax7,ax8],'Position',[.10 .11 .685 .815]);
            axis([ax1,ax4,ax7,ax8],[-126 90 -72 108]);
        end
        
        
        %---------------------------------------------%
        
        if isequal(PNGDir{1},1) && ~isempty(PNGDir{2})
            print(gcf, '-dpng', ['' PNGDir{2} '' filesep 'i' num2str(islice) '_j' num2str(jslice) '_k' num2str(kslice) '_sag.png']);
        end
        
    end
    
    if plotj
        %-----FIGURE 2 (axes 2,5,9,10)-----------------------------%
            f2=figure(fig2);
            ax2=axes('parent',f2);
            movegui(f2,'north');
            set(f2,'name',['' fig_names ', coronal view']);
        
        clear niftiplot_masked_imgdata
        for page=1:pages
            xdata2(:,page)=X(:,jpage,page);
            zdata2(:,page)=Z(:,jpage,page);
            imgdata2(:,page)=data(:,jpage,page);
            plot_imgdata2(:,page)=data_backdrop(:,jpage,page);
            niftiplot_masked_imgdata(:,page)=niftiplot_masked(:,jpage,page);
        end
        %orientate image like in fslview
        surf(ax2,xdata2,zdata2,plot_imgdata2,'EdgeColor','None');
        xlabel(ax2,'I');
        ylabel(ax2,'R');
        set(ax2,'Xdir','reverse');
        %text(24,66,'\leftarrow x=24,z=66');
        colormap(ax2,'gray');
        %axis equal
        view(ax2,2);

        xvalues=[];
        yvalues=[];
        scatter_dat=[];
        [rows,columns]=size(niftiplot_masked_imgdata);
        for row=1:rows
            for col=1:columns
                if niftiplot_masked_imgdata(row,col)~=0
                    xvalues=[xvalues xdata2(row,col)];
                    yvalues=[yvalues zdata2(row,col)];
                    scatter_dat=[scatter_dat niftiplot_masked_imgdata(row,col)];
                end
            end
        end


        ax5 =axes('parent',f2);
        scatter(ax5,xvalues,yvalues,18,scatter_dat,'filled')
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
        axis([ax2,ax5],[min(min(xdata2)) max(max(xdata2)) min(min(zdata2)) max(max(zdata2))]);
        %axis equal
        %text(-20,76,'\leftarrow y=-20,z=76');
        if ~ColorRange{3}
            
        elseif ColorRange{3}
            caxis(ax5,[colorrange_min,colorrange_max]);
        end

        if slicepos
            %ax9 for positioning the sagittal slice position
            ax9 =axes('parent',f2);
            %entries of X only change along first index of X
            %sagitt is position along R-L axis
            sagitt=X(ipage,1,1);
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
            ax10 =axes('parent',f2);
            %entries of Z only change along third index of Z
            %transv is position along I-S axis
            transv=Z(1,1,kpage);
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
        
        if isequal(PNGDir{1},1) && ~isempty(PNGDir{2})
            print(gcf, '-dpng', ['' PNGDir{2} '' filesep 'i' num2str(islice) '_j' num2str(jslice) '_k' num2str(kslice) '_cor.png']);
        end
        
    end

    
    if kpage
        %-----FIGURE 3 (axes 3,6,11,12)---------------------------%
            f3=figure(fig3);
            ax3=axes('parent',f3);
            movegui(f3,'northeast');
            set(f3,'name',['' fig_names ', transverse view']);

        surf(ax3,X(:,:,kpage),Y(:,:,kpage),data_backdrop(:,:,kpage),'EdgeColor','None');
        xlabel(ax3,'P');
        ylabel(ax3,'R');
        set(ax3,'Xdir','reverse');
        %text(-6,-54,'\leftarrow x=-6,y=-54');
        colormap(ax3,'gray');
        %axis equal
        view(ax3,2);

        clear niftiplot_masked_imgdata
        niftiplot_masked_imgdata=niftiplot_masked(:,:,kpage);
        xvalues=[];
        yvalues=[];
        scatter_dat=[];
        [rows,columns]=size(niftiplot_masked_imgdata);
        for row=1:rows
            for col=1:columns
                if niftiplot_masked_imgdata(row,col)~=0
                    xvalues=[xvalues X(row,col,kpage)];
                    yvalues=[yvalues Y(row,col,kpage)];
                    scatter_dat=[scatter_dat niftiplot_masked_imgdata(row,col)];
                end
            end
        end
        
        
        ax6 =axes('parent',f3);
        scatter(ax6,xvalues,yvalues,18,scatter_dat,'filled')
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
        %axis equal
        axis([ax3,ax6],[min(min(X(:,:,kpage))) max(max(X(:,:,kpage))) min(min(Y(:,:,kpage))) max(max(Y(:,:,kpage)))]);
        %text(-20,76,'\leftarrow y=-20,z=76');
        if ~ColorRange{3}
            
        elseif ColorRange{3}
            caxis(ax6,[colorrange_min,colorrange_max]);
        end

        if slicepos
            %ax11 for positioning the coronal (=frontal) slice position
            ax11 =axes('parent',f3);
            %entries of Y only change along second index of Y
            %coron is position along P-A axis
            coron=Y(1,jpage,1);
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
            ax12 =axes('parent',f3);
            %entries of X only change along first index of X
            %sagitt is position along R-L axis
            sagitt=X(ipage,1,1);
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
        
        
        if isequal(PNGDir{1},1) && ~isempty(PNGDir{2})
            print(gcf, '-dpng', ['' PNGDir{2} '' filesep 'i' num2str(islice) '_j' num2str(jslice) '_k' num2str(kslice) '_trans.png']);
        end
        
    end
  
end