function [] = performpca(matfile,storeredund,escape,calc_lagthreads)
fprintf('\n');
disp('PERFORMING PCA (performpca.m)');
% Performs PCA on TD.
m=load(['' matfile '']);

% Now compute TDz
disp('Executing computeTDz()');
[TDz,ColumnMeans,RowMeans]=computeTDz(m.TD);

if storeredund
    save(matfile,'TDz','-append');
end
save(matfile,'ColumnMeans','-append');
save(matfile,'RowMeans','-append');


%PERFORM PCA AND SAVE---------------------------------%
disp('Performing PCA on TDz.');
% Use pcacov on covariance matrix, which then returns the percentage
% of variance explained (which function pca() does not).

% Number of columns corresponds to number of variables.
% Number of rows corresponds to number of observations.
[COMPONENTS,EIGENVALUES,explained]=pcacov(cov(TDz));

save(matfile,'COMPONENTS','-append');
save(matfile,'EIGENVALUES','-append');
save(matfile,'explained','-append');

%----------------------------------------------------%


if calc_lagthreads
    %COMPUTE LAG THREAD TOPOGRAPHY-----------------------%
    [rows, ~] = size(TDz);

    %L contains lag thread topographies, ordered according to variance
    %explained.

    %TD(2,1)<0 if 1->2. This way, Column Means is the actual lag
    %projection map. Thus we have to multiply by negative 1:

    L=(-1)*(TDz*COMPONENTS)/sqrt(rows);

    if isfield(m,'comment')
        comment=strcat(m.comment,';','L=(-1)*(TDz*COMPONENTS)/sqrt(rows)');
    else
        comment='L=(-1)*(TDz*COMPONENTS)/sqrt(rows)';
    end

    save(matfile,'L','-append');
    save(matfile,'comment','-append');

    %%----------------------------------------------------%
end
disp('PCA Computation (performpca.m) done!');

if escape
    quit
end
end
