addpath '~/matlab_toolboxes/fieldtrip/'
addpath '~/matlab_toolboxes/shrinkIt/'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GET LIST OF SUBJECTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get subject IDs
cd('~/HCP/data/hpc/HCP500_Parcellation_Timeseries_Netmats/Results/node_timeseries/3T_Q1-Q6related468_MSMsulc_d300_ts2');
fnames = struct2table(dir()); %get list of file names
fnames = fnames(3:end,1); %remove . and .., remove other columns
fnames = table2array(fnames);
subjects = regexp(fnames, '^[0-9]*', 'match'); %remove ".txt".  becomes array of cells
subjects = table2array(cell2table(subjects));
S = numel(subjects);

cd('~/HCP/data/hpc/')
subjects2 = zeros(0);
disks = zeros(0);
for idisk=1:5
    
   d = strcat('disk',num2str(idisk));
   s = dir(d);
   
   %get only directories
   s = struct2table(s);
   dirs = table2array(s(:,4));
   names = table2array(s(dirs==1,1));
   %remove . and ..
   names = names(3:end); 
   subjects2 = [subjects2, names'];
   disks = [disks, repmat({d},[1 numel(names)])];
   
end

%remove subjects not included in ICA
[int, ia, ib] = intersect(subjects, subjects2');
disks = disks(ib)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ IN PARCELLATIONS AND CREATE INDICATOR ARRAY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_labels = '~/HCP/data/groupICA/groupICA_3T_Q1-Q6related468_MSMsulc_d300.ica/melodic_IC_ftb.dlabel.nii';
parcels = ft_read_cifti(fname_labels);  %structure with cell 'indexmax'
labels = parcels.indexmax;
csvwrite('~/srinc/labels.csv', labels);

%% CREATE INDICATOR MATRIX

%determine number of parcels in  X by Y by Z matrix named roi_mask
rois = 1:nrois_m;
V = size(parcels.indexmax, 1);
% VxQ binary indicator of parcel membership
parcels_ind = double(repmat(parcels.indexmax,1,nrois_m) == repmat(rois,V,1));
%pad with NaN for Q < 300
if(nrois_m < 300) parcels_ind(:,(nrois_m+1):300) = NaN; end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% READ IN ICA TIME COURSES (LR/RL/LR/RL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ICA_tcs = zeros(4800, 300, S); %ICA time courses resulting from dual regression

for isubj=1:S
    
    isubj

    cd('~/HCP/data/ICA_time_series/ICA300');

    %read in text file of ICA time courses, add to parcel_tcs
    fname = fnames(isubj);
    ICA_tcsi = readtable(fname{1}, 'Delimiter',' ','ReadVariableNames', false);
    ICA_tcs(:,:,isubj) = zscore(table2array(ICA_tcsi));
      
end

save('~/HCP/data/ICA_tcs300.mat','ICA_tcs','-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE CONNECTIVITY ESTIMATES & RELIABILITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% LOOP THROUGH DIFFERENT SCAN LENGTHS
%% COMPUTE RAW AND SHRINKAGE ESTIMATES OF CONNECTIVITY MATRIX
%% LATER: COMPUTE SHRINKAGE ESTIMATES OF CONNECTIVITY MATRIX (ONCE WE HAVE THE RELIABILITY THIS IS EASIER)
%% RELIABILITY = ABSOLUTE PERCENT ERROR & ICC

%% FOR Q=300 ONLY:
%% 1) COMPUTE SEED-LEVEL INTER-SESSION RELIABILITY 
%% 2) COMPUTE EDGE-LEVEL INTER-SESSION RELIABILITY 

lengths = [300, 600, 900, 1200, 1500, 1800, 2100, 2400];
L = numel(lengths);
varU = zeros(300, 300, L);
varX = zeros(300, 300, L);

%compute final estimate within each day
C_end = zeros(300,300,S,2); %dim4 = day 1, day 2
for isubj=1:S  
    C_end(:,:,isubj,1) = corrcoef(ICA_tcs(1:2400,:,isubj)); 
    C_end(:,:,isubj,2) = corrcoef(ICA_tcs(2401:4800,:,isubj)); 
end

for ii=1:L

    len = lengths(ii)

    %compute day 1 and day 2 connectivity estimates
    [C1, C2] = deal(zeros(300,300,S));
    part1 = 1:len;
    part2 = 2401:(2401+len-1);
    for isubj=1:S
        C1(:,:,isubj) = corrcoef(ICA_tcs(part1,:,isubj));
        C2(:,:,isubj) = corrcoef(ICA_tcs(part2,:,isubj));
    end

    varU_ii = (1/2)*var((C2 - C1), 0, 3); %computes var over 3rd dimension, second argument indicates that the denominator should be n-1
    varW_ii = (var(C1, 0, 3) + var(C2, 0, 3))/2;
    varX_ii = varW_ii - varU_ii;

    varU(:,:,ii) = varU_ii;
    varX(:,:,ii) = varX_ii;

end

varX(varX < 0) = 0; %threshold at zero

ICC = varX./(varU+varX);
for ii=1:L csvwrite(strcat('~/HCP/app/ICC',num2str(lengths(ii)),'.csv'), ICC(:,:,ii)); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISUALIZE OVERLAP BY NETWORK (WRITE CIFTI)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

overlap = csvread('~/reliability_app/overlap.csv');

fname_labels = strcat('~/HCP/data/groupICA/groupICA_3T_Q1-Q6related468_MSMsulc_d300.ica/melodic_IC_ftb.dlabel.nii');
parcels = ft_read_cifti(fname_labels);  %structure with cell 'indexmax'

cifti = parcels;
for ii = 1:300
    cifti.indexmax(parcels.indexmax==ii) = overlap(ii,2);
end
ft_write_cifti('~/reliability_app/overlap', cifti, 'parameter', 'indexmax');

%do by network
for k = 1:7
    overlap_k = csvread(strcat('~/reliability_app/overlap_net',num2str(k),'.csv'));
    cifti = parcels;
    for ii = 1:300
        cifti.indexmax(parcels.indexmax==ii) = overlap_k(ii,2);
    end
    fname = strcat('~/reliability_app/overlap_net',num2str(k));
    ft_write_cifti(fname, cifti, 'parameter', 'indexmax');
end






