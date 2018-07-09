function [gs,GSregressed_data] = performgsr(data,brainmask)
% Performs GSR on data. Every voxel in the brainmask-file whose intensity is
% Larger than 0 is considered part of the mask
% --------------------------
%
% Output:
%   gs : gs is the global signal (GS)
%   GSregressed_data : The data with the GS regressed.
%
% --------------------------

[maskdata,~,~,~,~] = readnifti(brainmask);

% If maskdata is 4D, take the first volumne
if length(size(maskdata))==4
   maskdata=maskdata(:,:,:,1);
end

% Check, if BOLDfile and maskfile are of same size (e.g. 91x109x91)
for i=1:3
   if size(data,i) ~= size(maskdata,i)
      error('BOLD-file and maskfile are of different matrix sizes'); 
   end
end

% gs is the average signal over the brainmask for each timestep (i.e. the Global Signal)
gs = zeros(1,size(data,4));

% maskvox is the number of voxels that have a value different from zero
maskvox = length(find(maskdata>0));
disp(['Number of non-zero voxels in maskfile: ' num2str(maskvox) '.']);

% data is ixjxkxl matrix
for l=1:size(data,4) % size(data,4) equals no. of timesteps
    % disp(['Entering step' num2str(l) '']);
    for i=1:size(data,1)
        for j=1:size(data,2)
            for k=1:size(data,3)
                if maskdata(i,j,k) > 0
                   if isnan(maskdata(i,j,k))
                       error('The maskfile contains NaN-entries.');
                   end
                   if isnan(data(i,j,k,l))
                       error('The imaging file contains NaN-entries.');
                   end
                   gs(1,l)=gs(1,l)+data(i,j,k,l)/maskvox;
                end
            end
        end
    end
end

% reshape data
data_res=reshape(data,[size(data,1)*size(data,2)*size(data,3) size(data,4)]);
   
disp(['Computing ' num2str(size(data,1)*size(data,2)*size(data,3)) ' regressions...']);

GSregressed_data_res=zeros(size(data,1)*size(data,2)*size(data,3),length(gs));

for i=1:size(data,1)*size(data,2)*size(data,3)
    if mod(i,10000)==0
        % disp(['Doing regression no. ' num2str(i) ' after ' num2str(toc) ' seconds.']);
    end
    
    if ~isequal(data_res(i,:),zeros(1,length(data_res(i,:))))
    
        % r is the vector of residuals
        [~,~,r]=regress(data_res(i,:)',[gs' ones(length(gs),1)]);
                   
        % GSregressed_data is the data vector with GS regressed out
        GSregressed_data_res(i,:)=r';
    end
        
end

GSregressed_data=reshape(GSregressed_data_res,[size(data,1), size(data,2), size(data,3), size(data,4)]);

end