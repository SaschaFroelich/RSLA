function []=man(topic)

clc;

switch topic
    case 'rsla'
        disp('This toolbox consists of three separate GUIs and some useful stand-alone ->functions. The GUIs are:');
        disp('(1) The TD-generation GUI "rslastarter" which generates TD-matrices from nifti-files. This GUI also provides a batch-creator (only for linux). Type "rslastarter" or "rsla init" to start the GUI and "man rslastarter" for help.');
        disp('(2) The correlation viewer "CorrViewer". It lets you display the correlations used to compute TD-matrices TD. Type "CorrViewer" to start the GUI and "man CorrViewer" for help.');
        disp('(3) The analysis GUI "rslaanalyzer". It lets you analyse the generated .mat-files containing the matrices TD and provides handy functionalities. Type "rslaanalyzer" or "rsla analyze" to start the GUI and "man rslaanalyzer" for help.');
    case 'rslastarter'
        fprintf('This GUI lets you compute the time-delay matrix TD from a 4D nifti file. It also allows you to create a bash-script for batch execution of a number of 4D nifti files. Type "rslastarter" or "rsla init" to open.\n Hover over fields and buttons in the GUI for further explanation.')
    case 'CorrViewer'
        fprintf('This GUI lets you analyse the lagged cross-covariance curves for every pair of ROIs. Type ''CorrViewer'' to open.\n')
    case 'rslaanalyzer'
        fprintf('This GUI lets you analyse and manipulate mat-files containing TD in various ways. Take not of the tooltipstring of the vrious buttons which appear when hovering above them with the mouse cursor.\n')
    case 'functions'
        fprintf('FUNCTION extract_psychometrics.m\n\nNo input needed\n\nExtracts the measures as saved in columns of a csv-file (from excel-file) and writes them in a cell array which is then saved as a mat-file under direction "dir".')
        fprintf('\n\nFUNCTION combineCSV.m\n\nNo input needed\n\nExtracts the measures as saved in columns of a csv-file (from excel-file) and writes them in a cell array which is then saved as a mat-file under direction "dir".')
    case 'PVL'
        fprintf(['\n\n\nPSYCHOMETRICS VS LAG provides the possibility to correlate (pairwise pearson correlation (->PEARSON) or partial correlation (->PARTIAL, see below)) any numerical measures\nwith the latency values of ROIs. Therefore, a csv-file is needed (can be created from excel-file) containing the numeric measures. \n\nGENERAL INFO',...
        '\n\nCSV-FILE\nThe csv-file is created from the corresponding excel-file (Save As -> .csv. Use field delimiter '','' (comma)). The first column has to contain the identifications of the subjects\n(e.g. the names) which need to be in accordance with the selected ->NII-FILES. NB: The data has to begin in the SECOND row! The first row contains the title of the columns,\nspecifying the names of the numeric measures. The following columns contain the numeric measures. There must not be gaps (i.e. blank cells) for the\ncorresponding measures. If a measure is not available for a subject, then this cell can contain NaN. NaN-entries are identified by the GUI and subsequently ignored for correlation.\n',...
        '\nNII-FILES\nThe corresponding latency values are extracted from the selected nii-files, using the ROIs as specified by the ->MASK. To correctly associate the .nii-file with the corresponding\nnumeric measure as listed in the ->csv-file, the subject identifier in the first column of the nii-file has to be identically present in the subject''s corresponding nii-file.\nExample: If a subject in the csv-file is called "subjM006", then the corresponding nifty-file could be "swarsubjM006lag_projection.nii" as it contains the string "subjM006".\n',...
        '\nMASK\nNeeds to be a .nii-file compatible with the selected ->NII-FILES. All voxels in the mask with values >0 are considered ROI-voxels. The values of those\nvoxels in the subject-level nii-files (->NII-FILES) are then averaged (arithemtic mean) for each subject to give a single "latency value". The obtained latency values\nare then correlated with their corresponding numeric measures as extracted from the csv-file.\n',...
        '\nPEARSON\nWhen choosing this option, all chosen measures will be correlated with the corresponding latency value of the ROI as defined by the ->MASK. All correlations will\nbe stored as scatterplots in the directory of your choosing.\n',...
        '\nPARTIAL\nIf this option is chosen, all selected measures will be pairwise correlated with the latency value of the ROI as defined by the ->MASK, controlling for all other selected measures.\n']);
        fprintf('At the end, three new variables will stored in the workspace: \n1) AllSbjcts: Contains the subjects that were used in the previous correlation analysis.\n')
        fprintf('2) latency: Contains the latency values of the subjects in the chosen cluster. The order of latency values corresponds to the subjects as given to the new variable AllSbjcts.\n')
        fprintf('3) psychometric_score: Contains the psychometric scores of the subjects (only the last analysed score). The order of latency values corresponds to the subjects as given to the new variable AllSbjcts.\n')
        fprintf('\n\n PSYCHOMETRICS VS LAG was tested for the .nii-file option of pearson correlation and works correctly (mat-groupfile option not yet tested for correctness, neither was partialcorr option). Test: An excel-file was created with 10 example subjects (nii-files). The psychometric measure "test" was added. This test-measure was made equal to the latency in a chosen cluster (right cerebellum). Said latency was computed via a different matlab-script and verified by using fslstats. The resulting correlation between the lag and that cluster using "Psychometrics vs Lag" was then computed to be r=1 with p=7e-28.')
    case 'DataStructure'
        fprintf('All .mat-files containing TD-matrices have to be stored in the SAME folder. That is, if it is desired to perform tests on two groups and then compare these two groups, the .mat-files containing TD have to be in the same folder (not in subfolders). Reason: The matfile for both groups contains just one variable "pathtofiles" applying to the files of group 1 and group 2.');
    case 'ToDo'
        fprintf('(1) Allow more user-friendliness in handling data structure (see ->DataStructure). \n');
        fprintf('(2) Make it possible to move files ("pathtofiles variable in mat-groupfiles").');
        fprintf('(3) Test correctness of PSYCHOMETRICS VS LAG for nii-files for partial correlation.');
    case 'PLOTSLICE'
        fprintf('The text windows next to the characters i,j and k specify the 2D-slice of the brain which has to be plotted. Choosing MATLAB notation a:b:c will plot all 2D-slices from slice a to slice c in steps of b. \n ');
    case 'PLOTNII'
        fprintf('Plot a nifti file. The plotted slices are those chosen in the text-wdit windows i,j and k. If dialog window ''pick a mask file'' appears: Mask file specifies which regions of the nifti latency file shall be plotted and thus overlaid over the anatomical background.');
    case 'CRTLAGMAPS'
        fprintf('Before hitting the button ''Create nifti'', check the following: \n');
        fprintf('i) As this will create unsmoothed as well as smoothed versions of the Lag-Map, make sure the FWHM of the Gaussian smoothing kernel in each direction is set to the desired value (the three text boxes next to ''FWHM [mm]''). \n');
        fprintf('ii) If the box ''Create Group-level nifti file'' is checked, this will create a group-level nifti file containing the average of all Lag-Maps created in this run. That groupfile will have the prefix ''GROUPFILE''. \n');
        fprintf('Upon hitting the button ''Create nifti'', these are the steps to follow: \n');
        fprintf('1) Select the mat-files containing TD-matrices from which to create a Lag-Map (i.e. seeded cascade). For every mat-file, a nifti-file containing the seeded cascade will be created (the corresponding seeds will have to be chosen later on.). \n');
        fprintf('2) Select the directory where to generaed Lag-Map nifti files shall be stored. \n');
        fprintf('3) Select 3D nifti-file containing the seed for the Lag-Map (this file is called the ''mask''). Every voxel of intensity~=0 (i.e. not equal to zero) in the mask will be used as part of the seed-region. Voxels containing a single seed-region may be disjoint. NB: The intensity of the voxels in the mask does not matter as long as the intensities are ~=0, i.e. if voxel A has intensity of 2, voxel B has intensity of 1000, and voxel C intensity of -10, they will constitute the same seed region. \n');
    case 'maskfile'
        fprintf('The time delay matrix TD (n by m matrix) has n rows and m columns. Each row corresponds to a brain region as defined in the first maskfile,\nwhile each column corresponds to a brain region as defined by the second maskfile. All non-zero voxels with the same intensity value within a maskfile are considered as one ROI, or region. \nIf no second maskfile is selected, then TD will be of size nxn and skew-symmetric.\n');
    case 'parallelize'
        fprintf('Only for UNIX systems. This function enables the parallel computation of TD on a number of servers. This can be helpful as TD computation for a large number of ROIs can be time-consuming.');
        fprintf('\nIn order to use this option, the username under which the computation is performed has to be identical on all servers. Furthermore, if login to the servers requires a password certificates must be provided to the servers in order to circumvent entry request for password.');
    case 'PERMTEST'
        fprintf('If the original brainmask for TD is not continuous (i.e. it consists of separate regions), statistical signficance testing using SPM is not possible. Use a permutation test instead.')
    case 'Alert 1'
        disp('...');
    case 'Error 1'
        disp('...');
    otherwise
        disp('Not a valid option. Type "man rsla" for a general introduction.');
end