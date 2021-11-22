%% Eigenvector and IND (1000 permutations; 74 subjects)
%calculate actual correlations
cd /Users/Qinggang1/Documents/MATLAB/MR/MR_files_from_linux/
load('Eigen_243_bipos0130.mat')
load('rsfmri_cov.mat')
cd community_new2/
load('community_Louvain_243pos.mat')
Eigen_FPCN = Eigen_243_bipos0130(community_243(:,6) == 1, :, :);
Eigen_FPCN_avg = mean(Eigen_FPCN, 1);
Eigen_FPCN_reshaped = reshape(Eigen_FPCN_avg, 74,30);
r_actual = zeros(30,1);
for thresh = 1:size(Eigen_FPCN_reshaped, 2)
    matrix = [rsfmricov(:,19) Eigen_FPCN_reshaped(:,thresh) rsfmricov(:, 2:8) rsfmricov(:,14)];
    [rho, p] = partialcorr(matrix);
    r_actual(thresh) = rho(1,2);
end

%permute IV
permuted_IND = zeros(74, 1000);
IND = rsfmricov(:,19);
for perm = 1:1000
    curr_perm = IND(randperm(length(IND)));
    permuted_IND(:, perm) = curr_perm;
end
    

%calculate permuted correlations
r_permuted = zeros(30,1000);
for run = 1:1000
    for thresh = 1:size(Eigen_FPCN_reshaped, 2)
        matrix = [permuted_IND(:, run) Eigen_FPCN_reshaped(:, thresh) rsfmricov(:, 2:8) rsfmricov(:, 14)];
        [rho, p] = partialcorr(matrix);
        r_permuted(thresh, run) = rho(1,2);
    end
end

%find correlation at 95th percentile
[r_max_permuted, i_max_permuted] = max(r_permuted);
r_sorted = sort(r_max_permuted, 'descend');
r_crit = r_sorted(50);

%calculate permuted AUC for the critical value
AUC_crit_list = [];
included_line = [];
for i = 1:size(r_permuted, 2)
    curr_r_p = r_permuted(:, i);
    if max(curr_r_p) > r_crit
        curr_r_p(curr_r_p <= r_crit) = r_crit;
        curr_r_p = curr_r_p - r_crit;
        curr_AUC = trapz(curr_r_p);
        included_line = [included_line curr_r_p];
        AUC_crit_list = [AUC_crit_list; curr_AUC];
    end
end
AUC_crit = mean(AUC_crit_list);
        
%calculate actual AUC
AUC_actual = 0;
if max(r_actual) > r_crit
    r_adjusted = r_actual;
    r_adjusted(r_adjusted <= r_crit) = r_crit;
    r_adjusted = r_adjusted - r_crit;
    AUC_actual = trapz(r_adjusted);
end
%% Eigenvector and IND (1000 permutations; 110 subjects)
%calculate actual correlations
cd /Users/Qinggang1/Documents/MATLAB/MR/MR_files_from_linux/
load('Eigen_208_bipos0130.mat')
load('rsfmri110_cov.mat')
cd community_new/
load('community_Louvain_208pos.mat')
Eigen_FPCN = Eigen_208_bipos0130(group_assign(:,6) == 3, :, :);
Eigen_FPCN_avg = mean(Eigen_FPCN, 1);
Eigen_FPCN_reshaped = reshape(Eigen_FPCN_avg, 110,30);
r_actual = zeros(30,1);
for thresh = 1:size(Eigen_FPCN_reshaped, 2)
    matrix = [rsfmri110cov(:,19) Eigen_FPCN_reshaped(:,thresh) rsfmri110cov(:, 2:8) rsfmri110cov(:,14)];
    [rho, p] = partialcorr(matrix);
    r_actual(thresh) = rho(1,2);
end

%permute IV
permuted_IND = zeros(110, 1000);
IND = rsfmri110cov(:,19);
for perm = 1:1000
    curr_perm = IND(randperm(length(IND)));
    permuted_IND(:, perm) = curr_perm;
end
    

%calculate permuted correlations
r_permuted = zeros(30,1000);
for run = 1:1000
    disp(run)
    for thresh = 1:size(Eigen_FPCN_reshaped, 2)
        matrix = [permuted_IND(:, run) Eigen_FPCN_reshaped(:, thresh) rsfmri110cov(:, 2:8) rsfmri110cov(:, 14)];
        [rho, p] = partialcorr(matrix);
        r_permuted(thresh, run) = rho(1,2);
    end
end

%find correlation at 95th percentile
[r_max_permuted, i_max_permuted] = max(r_permuted);
r_sorted = sort(r_max_permuted, 'descend');
r_crit = r_sorted(50);

%calculate permuted AUC for the critical value
AUC_crit_list = [];
included_line = [];
for i = 1:size(r_permuted, 2)
    curr_r_p = r_permuted(:, i);
    if max(curr_r_p) > r_crit
        curr_r_p(curr_r_p <= r_crit) = r_crit;
        curr_r_p = curr_r_p - r_crit;
        curr_AUC = trapz(curr_r_p);
        included_line = [included_line curr_r_p];
        AUC_crit_list = [AUC_crit_list; curr_AUC];
    end
end
AUC_crit = mean(AUC_crit_list);
        
%calculate actual AUC
AUC_actual = 0;
if max(r_actual) > r_crit
    r_adjusted = r_actual;
    r_adjusted(r_adjusted <= r_crit) = r_crit;
    r_adjusted = r_adjusted - r_crit;
    AUC_actual = trapz(r_adjusted);
end
%% Eigenvector and IND (1000 permutations; 74 subjects)
% Each ROI separately

%calculate actual correlations
cd /Users/Qinggang1/Documents/MATLAB/MR/MR_files_from_linux/
load('Eigen_243_bipos0130.mat')
load('rsfmri_cov.mat')
cd community_new2/
load('community_Louvain_243pos.mat')
r_actual = zeros(243,30);
for thresh = 1:size(Eigen_243_bipos0130, 3)
    for roi = 1:size(Eigen_243_bipos0130, 1)
        matrix = [rsfmricov(:,19) Eigen_243_bipos0130(roi,:,thresh)' rsfmricov(:, 2:8) rsfmricov(:,14)];
        [rho, p] = partialcorr(matrix);
        r_actual(roi,thresh) = rho(1,2);
    end
end

%permute IV
permuted_IND = zeros(74, 1000);
IND = rsfmricov(:,19);
for perm = 1:1000
    curr_perm = IND(randperm(length(IND)));
    permuted_IND(:, perm) = curr_perm;
end

%calculate permuted correlations (WARNING: takes very long!)
r_permuted = zeros(243,30,1000);
for roi = 1:size(Eigen_243_bipos0130, 1)
    disp(roi)
    for run = 1:1000
        for thresh = 1:size(Eigen_243_bipos0130, 3)
            matrix = [permuted_IND(:, run) Eigen_243_bipos0130(roi, :, thresh)' rsfmricov(:, 2:8) rsfmricov(:, 14)];
            [rho, p] = partialcorr(matrix);
            r_permuted(roi, thresh, run) = rho(1,2);
        end
    end
end

%find correlation at 95th percentile
[r_max_permuted, i_max_permuted] = max(r_permuted, [], [1 2], 'linear');
r_max_permuted = reshape(r_max_permuted, 1000, 1);
i_max_permuted = reshape(i_max_permuted, 1000, 1);
r_sorted = sort(r_max_permuted, 'descend');
r_crit = r_sorted(50);
[roi_max_p, thresh_max_p, run_max_p] = ind2sub([243,30,1000],i_max_permuted);

%calculate permuted AUC for the critical value
AUC_crit_list = [];
included_line = [];
for run = 1:size(r_permuted, 3)
    for roi = 1:size(r_permuted, 1)
        curr_r_p = r_permuted(roi, :, run);
        if max(curr_r_p) > r_crit
            curr_r_p(curr_r_p <= r_crit) = r_crit;
            curr_r_p = curr_r_p - r_crit;
            curr_AUC = trapz(curr_r_p);
            included_line = [included_line curr_r_p];
            AUC_crit_list = [AUC_crit_list; curr_AUC];
        end
    end
end
AUC_crit = mean(AUC_crit_list);

%calculate actual AUC and compare it with critical AUC
sig_roi_list = [];
for roi = 1:size(Eigen_243_bipos0130, 1)
    AUC_actual = 0;
    if max(r_actual(roi,:)) > r_crit
        r_adjusted = r_actual;
        r_adjusted(r_adjusted <= r_crit) = r_crit;
        r_adjusted = r_adjusted - r_crit;
        AUC_actual = trapz(r_adjusted);
        if AUC_actual > AUC_crit
            sig_roi_list = [sig_roi_list; community_243(roi,:)];
        end
    end
end
%% Biomarkers and modularity (1000 permutations; 74 subjects)

cd /Users/Qinggang1/Documents/MATLAB/MR/MR_files_from_linux/
load('Qlouvain_243_bipos0110.mat')
load('rsfmri_cov.mat')
r_actual = zeros(10,1);
med = zeros(74,1);
med(rsfmricov(:, 11) == 1 | rsfmricov(:,12) == 1 | rsfmricov(:,13) == 1) = 1;
for thresh = 1:size(Qlouvain_243_bipos0110, 2)
    curr_Q = Qlouvain_243_bipos0110(:,thresh) .^ 4;
    %curr_Q = Qlouvain_243_bipos0110(:,thresh);
    %matrix = [rsfmricov(:,17) curr_Q rsfmricov(:, 2:7) rsfmricov(:,9:15)];
    matrix = [rsfmricov(:,17) curr_Q rsfmricov(:, 2:7) med rsfmricov(:,9:10) rsfmricov(:,14:15)];
    [rho, p] = partialcorr(matrix);
    r_actual(thresh) = rho(1,2);
end

%permute IV
permuted_bio = zeros(74, 1000);
biomarker = rsfmricov(:,17);
for perm = 1:1000
    curr_perm = biomarker(randperm(length(biomarker)));
    permuted_bio(:, perm) = curr_perm;
end

%calculate permuted correlations
r_permuted = zeros(10,1000);
for run = 1:1000
    disp(run)
    for thresh = 1:size(Qlouvain_243_bipos0110, 2)
        curr_Q = Qlouvain_243_bipos0110(:, thresh) .^ 4;
        matrix = [permuted_bio(:, run) curr_Q rsfmricov(:, 2:7) med rsfmricov(:,9:10) rsfmricov(:,14:15)];
        [rho, p] = partialcorr(matrix);
        r_permuted(thresh, run) = rho(1,2);
    end
end

%find correlation at 95th percentile
[r_min_permuted, i_min_permuted] = min(r_permuted);
r_sorted = sort(r_min_permuted, 'ascend');
r_crit = r_sorted(50);

%calculate permuted AUC for the critical value
AUC_crit_list = [];
included_line = [];
for i = 1:size(r_permuted, 2)
    curr_r_p = r_permuted(:, i);
    if min(curr_r_p) < r_crit
        curr_r_p(curr_r_p >= r_crit) = r_crit;
        curr_r_p = curr_r_p - r_crit;
        curr_AUC = abs(trapz(curr_r_p));
        included_line = [included_line curr_r_p];
        AUC_crit_list = [AUC_crit_list; curr_AUC];
    end
end
AUC_crit = mean(AUC_crit_list);

%calculate actual AUC
AUC_actual = 0;
if min(r_actual) < r_crit
    r_adjusted = r_actual;
    r_adjusted(r_adjusted >= r_crit) = r_crit;
    r_adjusted = r_adjusted - r_crit;
    AUC_actual = abs(trapz(r_adjusted));
end
%% Biomarkers and modularity (1000 permutations; 74 subjects)
cd /Users/Qinggang1/Documents/MATLAB/MR/MR_files_from_linux/
load('Qlouvain_243_bipos0130.mat')
load('rsfmri_cov.mat')
r_actual = zeros(30,1);
med = zeros(74,1);
med(rsfmricov(:, 11) == 1 | rsfmricov(:,12) == 1 | rsfmricov(:,13) == 1) = 1;
for thresh = 1:size(Q_243_bipos0130, 2)
    curr_Q = Q_243_bipos0130(:,thresh) .^ 4;
    %curr_Q = Qlouvain_243_bipos0110(:,thresh);
    %matrix = [rsfmricov(:,17) curr_Q rsfmricov(:, 2:7) rsfmricov(:,9:15)];
    matrix = [rsfmricov(:,17) curr_Q rsfmricov(:, 2:7) med rsfmricov(:,9:10) rsfmricov(:,14:15)];
    [rho, p] = partialcorr(matrix);
    r_actual(thresh) = rho(1,2);
end

%permute IV
permuted_bio = zeros(74, 1000);
biomarker = rsfmricov(:,17);
for perm = 1:1000
    curr_perm = biomarker(randperm(length(biomarker)));
    permuted_bio(:, perm) = curr_perm;
end

%calculate permuted correlations
r_permuted = zeros(30,1000);
for run = 1:1000
    disp(run)
    for thresh = 1:size(Q_243_bipos0130, 2)
        curr_Q = Q_243_bipos0130(:, thresh) .^ 4;
        matrix = [permuted_bio(:, run) curr_Q rsfmricov(:, 2:7) med rsfmricov(:,9:10) rsfmricov(:,14:15)];
        [rho, p] = partialcorr(matrix);
        r_permuted(thresh, run) = rho(1,2);
    end
end

%find correlation at 95th percentile
[r_min_permuted, i_min_permuted] = min(r_permuted);
r_sorted = sort(r_min_permuted, 'ascend');
r_crit = r_sorted(50);

%calculate permuted AUC for the critical value
AUC_crit_list = [];
included_line = [];
for i = 1:size(r_permuted, 2)
    curr_r_p = r_permuted(:, i);
    if min(curr_r_p) < r_crit
        curr_r_p(curr_r_p >= r_crit) = r_crit;
        curr_r_p = curr_r_p - r_crit;
        curr_AUC = abs(trapz(curr_r_p));
        included_line = [included_line curr_r_p];
        AUC_crit_list = [AUC_crit_list; curr_AUC];
    end
end
AUC_crit = mean(AUC_crit_list);

%calculate actual AUC
AUC_actual = 0;
if min(r_actual) < r_crit
    r_adjusted = r_actual;
    r_adjusted(r_adjusted >= r_crit) = r_crit;
    r_adjusted = r_adjusted - r_crit;
    AUC_actual = abs(trapz(r_adjusted));
end
