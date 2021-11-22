%% Local efficiency across [0.01-0.10]
cd /home/qinggang/research/MIDUS_refresher/resting_conn
load('power_e74_thresh50_30zed.mat')
%cd community
%load('community_Louvain_209pos.mat')
%ci = com_Louvain_209pos(:,3);
threshold = 0.01:0.01:0.10;
Eloc_243_bipos0110 = [];
for cost = threshold
    disp(cost)
    Eloc = [];
    for i = 1:size(power_e74_thresh50_30zed, 3)
        curr_sub = power_e74_thresh50_30zed(:, :, i);
        curr_sub_clean = weight_conversion(threshold_proportional(curr_sub, cost),'binarize');
        %Eloc = [Eloc efficiency_bin(curr_sub_clean, 1)];
        Eloc = [Eloc; mean(efficiency_bin(curr_sub_clean, 1))];
    end
    %Eloc_list = cat(3, Eloc_list, Eloc);
    Eloc_243_bipos0110 = [Eloc_243_bipos0110 Eloc];
end
%% Correlation with local efficiency by community
cd /home/qinggang/research/MIDUS_refresher/resting_conn
load('MR_resting_biomarker.mat')
load('Eloc_209_bipos0110.mat')
cd community
load('community_Louvain_209pos.mat')
r_Eloc_0110 = [];
p_Eloc_0110 = [];
for i = 1:size(Eloc_list, 3)
    curr_thresh = Eloc_list(:, :, i);
    r_Eloc = [];
    p_Eloc = [];
    for sys = 1:5
        curr_sys = curr_thresh(com_Louvain_209pos(:, 3) == sys, :);
        curr_sys_avg = mean(curr_sys, 1)';
        MR_AGB = [MR_biomarker(:,12) [curr_sys_avg(1:112); curr_sys_avg(114:119)] MR_biomarker(:,2) MR_biomarker(:,3)];
        [rho, p] = partialcorr(MR_AGB);
        r_Eloc = [r_Eloc; rho(1,2)];
        p_Eloc = [p_Eloc; p(1,2)];
    end
    r_Eloc_0110 = [r_Eloc_0110 r_Eloc];
    p_Eloc_0110 = [p_Eloc_0110 p_Eloc];
end

%% Correlation with local efficiency across all nodes

cd /home/qinggang/research/MIDUS_refresher/resting_conn
load('MR_resting_biomarker.mat')
load('Eloc_209_bipos0110.mat')
r_Eloc_0110 = [];
p_Eloc_0110 = [];
for i = 1:size(Eloc_list, 3)
    curr_thresh = Eloc_list(:, :, i);
    Eloc = mean(curr_thresh, 1)';
    MR_AGB = [MR_biomarker(:,12) [Eloc(1:112); Eloc(114:119)] MR_biomarker(:,2) MR_biomarker(:,3)];
    [rho, p] = partialcorr(MR_AGB);
    r_Eloc_0110 = [r_Eloc_0110 rho(1,2)];
    p_Eloc_0110 = [p_Eloc_0110 p(1,2)];
end

%% Global efficiency all nodes across [0.01-0.10]
cd /home/qinggang/research/MIDUS_refresher/resting_conn
load('power_thresh50_25zed.mat')
threshold = 0.01:0.01:0.10;
Eglob_list = [];
for cost = threshold
    Eglob = [];
    for i = 1:size(power_e74_thresh50_30zed, 3)
        curr_sub = power_e74_thresh50_30zed(:, :, i);
        curr_sub_clean = weight_conversion(threshold_proportional(curr_sub, cost),'binarize');
        Eglob = [Eglob; efficiency_bin(curr_sub_clean)];
    end
    Eglob_list = [Eglob_list Eglob];
end
%% Eglob of communities across [0.01-0.10]
cd /home/qinggang/research/MIDUS_refresher/resting_conn/community/
load('com1_209pos.mat')
load('com2_209pos.mat')
load('com3_209pos.mat')
load('com4_209pos.mat')
load('com5_209pos.mat')
communities = {com1; com2; com3; com4; com5};
Eglob_com = [];
for g = 1:size(communities, 1)
    curr_graph = communities{g, 1};
    Eloc_list = [];
    threshold = 0.01:0.01:0.10;
    for cost = threshold
        Eloc = [];
        for i = 1:size(curr_graph, 3)
            curr_sub = curr_graph(:, :, i);
            curr_sub_clean = weight_conversion(threshold_proportional(curr_sub, cost),'binarize');
            Eloc = [Eloc; efficiency_bin(curr_sub_clean)];
        end
        Eloc_list = [Eloc_list Eloc];
    end
    Eglob_com = cat(3, Eglob_com, Eloc_list);
end
%% Eloc across communities

cd /home/qinggang/research/MIDUS_refresher/resting_conn/community_new2/
load('com1_243pos.mat')
load('com2_243pos.mat')
load('com3_243pos.mat')
load('com4_243pos.mat')
load('com5_243pos.mat')
communities = {com1; com2; com3; com4; com5};
Eloc_com = [];
for g = 1:size(communities, 1)
    curr_graph = communities{g, 1};
    Eloc_list = [];
    threshold = 0.01:0.01:0.10;
    for cost = threshold
        Eloc = [];
        for i = 1:size(curr_graph, 3)
            curr_sub = curr_graph(:, :, i);
            curr_sub_clean = weight_conversion(threshold_proportional(curr_sub, cost),'binarize');
            Eloc = [Eloc; mean(efficiency_bin(curr_sub_clean, 1))];
        end
        Eloc_list = [Eloc_list Eloc];
    end
    Eloc_com = cat(3, Eloc_com, Eloc_list);
end
%% Eloc across 0.01 to 0.30 (243 nodes) all nodes

%cd /home/qinggang/research/MIDUS_refresher/resting_conn
cd /Users/Qinggang1/Documents/MATLAB/MR/MR_files_from_linux/community_new_fd/
load('powerfd_e72_th50_10.mat')
threshold = 0.01:0.01:0.30;
Eloc_243_bipos0130 = [];
for cost = threshold
    disp(cost)
    Eloc = [];
    for i = 1:size(powerfd_e72_th50, 3)
        curr_sub = powerfd_e72_th50(:, :, i);
        curr_sub_clean = weight_conversion(threshold_proportional(curr_sub, cost),'binarize');
        Eloc = [Eloc efficiency_bin(curr_sub_clean, 1)];
        %Eigen = [Eigen; mean(efficiency_bin(curr_sub_clean, 1))];
    end
    Eloc_243_bipos0130 = cat(3, Eloc_243_bipos0130, Eloc);
    %Eigen_243_bipos0110 = [Eigen_243_bipos0110 Eloc];
end

%cd community_new2/
load('community_Louvain231_own.mat')
Eloc_243_m = [];
for i = 1:5
    curr_m = Eloc_243_bipos0130(community_Louvain_231pos(:, 6) == i, :, 1:20);
    disp(size(curr_m,1))
    curr_m_0210 = mean(curr_m, 3);
    curr_m_avg = mean(curr_m_0210, 1);
    Eloc_243_m = [Eloc_243_m curr_m_avg'];
end

Eloc_243 = mean(Eloc_243_bipos0130(:,:, 1:20), 3);
Eloc_243 = mean(Eloc_243, 1)';
%% Correlation with global efficiency all nodes
cd /home/qinggang/research/MIDUS_refresher/resting_conn
load('MR_resting_biomarker.mat')
load('Eglob_209_bipos0110.mat')
r_Eglob_0110 = [];
for i = 1:size(Eloc_list, 2)
    MR_AGB = [MR_biomarker(:,12) [Eloc_list(1:112, i); ...
        Eloc_list(114:119, i)] MR_biomarker(:,2) MR_biomarker(:,3)];
    rho = partialcorr(MR_AGB);
    r_Eglob_0110 = [r_Eglob_0110; rho(1,2)];
end

%% Correlation with Eglob of communities
cd /home/qinggang/research/MIDUS_refresher/resting_conn
load('MR_resting_biomarker.mat')
cd community/
load('Eglob_communities.mat')
r_Eglob_com = [];
for com = 1:size(Eglob_com, 3)
    curr_com = Eglob_com(:, :, com);
    r_Eglob_0110 = [];
    for i = 1:size(curr_com, 2)
        MR_AGB = [MR_biomarker(:,12) [curr_com(1:112, i); curr_com(114:119, i)] MR_biomarker(:,2) MR_biomarker(:,3)];
        rho = partialcorr(MR_AGB);
        r_Eglob_0110 = [r_Eglob_0110; rho(1,2)];
    end
    r_Eglob_com = [r_Eglob_com r_Eglob_0110];
end
%% Eglob and Eloc of random graphs (0.2) 100 iterations
cd /home/qinggang/research/MIDUS_refresher/resting_conn
load('power_thresh50_25zed.mat')
threshold_range = [0.49 0.50];
Eglob_rand = [];
Eloc_rand = [];
for thres = threshold_range
    disp(thres)
    allsub_thres = create_rand_G(power_thresh50_25, thres);
    Eglob_list = [];
    Eloc_list = [];
    for perm = 1:100
        Eglob = [];
        Eloc = [];
        r = randi([1 50], 1, 119);
        for i = 1:size(r, 2)
            curr_network = allsub_thres(:,:,r(i),i);
            Eglob = [Eglob; efficiency_bin(curr_network)];
            Eloc = [Eloc; mean(efficiency_bin(curr_network, 1))];
        end
        Eglob_list = [Eglob_list; mean(Eglob)];
        Eloc_list = [Eloc_list; mean(Eloc)];
    end
    Eglob_rand = [Eglob_rand mean(Eglob_list)];
    Eloc_rand = [Eloc_rand mean(Eloc_list)];
end

%% Get component (number of disconnected nodes)
%cd /home/qinggang/research/MIDUS_refresher/resting_conn
cd /Users/Qinggang1/Documents/MATLAB/MR/MR_files_from_linux/
load('power_e74_thresh50_30zed.mat')
threshold = [0.01: 0.01: 0.20];
c_list = [];
for sub = 1:size(power_e74_thresh50_30zed, 3)
    curr_sub = power_e74_thresh50_30zed(:,:,sub);
    c_sub = [];
    for cost = threshold
        curr_sub_clean = weight_conversion(threshold_proportional(curr_sub, cost),'binarize');
        [c, c_size] = get_components(curr_sub_clean);
        tbl = tabulate(c_size);
        c_sub = [c_sub tbl(1, 2)];
    end
    c_list = [c_list; c_sub];
end
c_avg = mean(c_list, 1);

%% Average pos and neg weights across participants
cd /home/qinggang/research/MIDUS_refresher/resting_conn
load('power_thresh50_25zed.mat')
pos_list = [];
neg_list = [];
for sub = 1:size(power_thresh50_25, 3)
    curr_sub = power_thresh50_25(:,:,sub);
    curr_pos = curr_sub(curr_sub > 0);
    curr_neg = curr_sub(curr_sub < 0);
    pos_list = [pos_list; mean(curr_pos)];
    neg_list = [neg_list; mean(curr_neg)];
end

%% Modularity across gamme value (209 nodes binary positive only)
cd /Users/Qinggang1/Documents/MATLAB/MR/MR_files_from_linux/
load('power_e74_thresh50_30zed.mat')
cd community_new_fd/
load('powerfd_e72_th50_10.mat')
cd ..

g_threshold = 0.75:0.05:1.25;
Q243_gammavar = [];
for gamma = g_threshold
    disp(gamma)
    Q_gamma = cal_Q_louvain(power_e74_thresh50_30zed, [0.01:0.01:0.30], gamma);
    Q243_gammavar = [Q243_gammavar mean(mean(Q_gamma, 2))];
end
save('Qlouvain243_gammavar.mat', 'Q243_gammavar')
%% Modularity at 1 gamma value (243 nodes)
cd /Users/Qinggang1/Documents/MATLAB/MR/MR_files_from_linux/
load('power_e74_thresh50_30zed.mat')
Q_gamma_1 = cal_Q_louvain(power_e74_thresh50_30zed, [0.01:0.01:0.30], 1);

%% New modularity and biomarkers correlation
cd /Users/Qinggang1/Documents/MATLAB/MR/MR_files_from_linux/
load('Q_243_bipos0130_gamma125.mat')
load('rsfmri_cov.mat')
Q_210 = Q_gamma_125(:, 2:10);
Q_210_avg = mean(Q_210, 2);
Q_210_normed = Q_210_avg .^ 4;
med = zeros(74, 1);
med(rsfmricov(:,11) == 1 | rsfmricov(:,12) == 1 | rsfmricov(:,13) == 1) = 1;
matrix = [Q_210_normed rsfmricov(:,17) rsfmricov(:, 2:7) rsfmricov(:, 9:10) med rsfmricov(:, 14:15)];
[rho, p] = partialcorr(matrix);



%% Brain signal variability whole brain
cd /home/qinggang/research/MIDUS_refresher/resting_conn/'conn_MR_new (another copy)'/results/preprocessing/
GM_times = [];
GM_times_sd = [];
for sub = 1:119
    filename = sprintf('ROI_Subject%03d_Condition001.mat', sub);
    timeseries = load(filename, 'data');
    GM_times = [GM_times timeseries.data{1,1}];
    GM_times_sd = [GM_times_sd; std(timeseries.data{1,1})];
end
cd /home/qinggang/research/MIDUS_refresher/resting_conn/
load('excluded_scan002.mat')
GMtimes_sd_74 = GM_times_sd(excluded_scan < 118);
%% Brain signal variability by power network (FD < 0.2 mm)
cd /home/qinggang/research/MIDUS_refresher/resting_conn/'conn_MR_new (another copy)'/results/preprocessing/
GM_times = zeros(236,264,119);
GM_times_sd = zeros(119,264);
for sub = 1:119
    disp(sub)
    filename = sprintf('ROI_Subject%03d_Condition001.mat', sub);
    timeseries = load(filename, 'data');
    roi_times = zeros(236, 264);
    roi_times_sd = zeros(1,264);
    for roi = 168:431
        roi_times(:, roi-167) = timeseries.data{1,roi};
        roi_times_sd(roi-167) = std(timeseries.data{1,roi});
    end
    GM_times(:,:,sub) = roi_times;
    GM_times_sd(sub, :) = roi_times_sd;
end
cd /home/qinggang/research/MIDUS_refresher/resting_conn/
load('excluded_scan002.mat')
cd /home/qinggang/research/MIDUS_refresher/resting_conn/community_new2/
load('power_e2_valid50.mat')
load('community_Louvain_243pos.mat');
GMtimes_sd_filtered = GM_times_sd(excluded_scan < 118, valid_vx(:,2) > 29);
GMtimes_sd_network = zeros(74,5);
for nw = 1:5
    curr_nw = GMtimes_sd_filtered(:, community_243(:,6) == nw);
    disp(size(curr_nw, 2))
    GMtimes_sd_network(:, nw) = mean(curr_nw, 2);
end

%BNV assignment
load('Power243_BNV_assignment.mat');
GMtimes_sd_network = zeros(74,5);
for nw = 1:5
    curr_nw = GMtimes_sd_filtered(:, Power243_BrainNetViewer(:,6) == nw);
    disp(size(curr_nw, 2))
    GMtimes_sd_network(:, nw) = mean(curr_nw, 2);
end
%% Brain signal variability by power network (FD < 0.5 mm)
cd /home/qinggang/research/MIDUS_refresher/resting_conn/'conn_MR_new (copy)'/results/preprocessing/
GM_times = zeros(236,264,119);
GM_times_sd = zeros(119,264);
for sub = 1:119
    disp(sub)
    filename = sprintf('ROI_Subject%03d_Condition001.mat', sub);
    timeseries = load(filename, 'data');
    roi_times = zeros(236, 264);
    roi_times_sd = zeros(1,264);
    for roi = 168:431
        roi_times(:, roi-167) = timeseries.data{1,roi};
        roi_times_sd(roi-167) = std(timeseries.data{1,roi});
    end
    GM_times(:,:,sub) = roi_times;
    GM_times_sd(sub, :) = roi_times_sd;
end
cd /home/qinggang/research/MIDUS_refresher/resting_conn/
load('excluded_scan.mat', 'data')
cd /home/qinggang/research/MIDUS_refresher/resting_conn/community_new/
load('power_valid_50percent.mat')
load('community_Louvain_208pos.mat');
%also clears out the sub with missing data
GMtimes_sd_filtered = GM_times_sd(data < 74 & data ~= 60, valid_vx(:,2) > 24);
GMtimes_sd_network = zeros(110,5);
for nw = 1:5
    curr_nw = GMtimes_sd_filtered(:, group_assign(:,6) == nw);
    disp(size(curr_nw, 2))
    GMtimes_sd_network(:, nw) = mean(curr_nw, 2);
end

%BNV assignment
load('Power208_BNV_assignment.mat');
GMtimes_sd_network = zeros(110,5);
for nw = 1:5
    curr_nw = GMtimes_sd_filtered(:, Power208_BrainNetViewer(:,6) == nw);
    disp(size(curr_nw, 2))
    GMtimes_sd_network(:, nw) = mean(curr_nw, 2);
end

%% Correlation Brain signal var & biomarkers by region
cd /home/qinggang/research/MIDUS_refresher/resting_conn/
load('GMtimes_sd_74_243.mat');
load('rsfmri_cov.mat');
results = zeros(2, 243);
for i = 1:243
    [r, p] = corrcoef(rsfmricov(:,17), GMtimes_sd_filtered(:, i));
    results(1, i) = r(1,2);
    results(2, i) = p(1,2);
end
results_cov = zeros(2, 243);
for i = 1:243
    matrix = [GMtimes_sd_filtered(:,i) rsfmricov(:,17) rsfmricov(:, 2:7) rsfmricov(:, 9:15)];
    [rho, p] = partialcorr(matrix);
    results_cov(1, i) = rho(1,2);
    results_cov(2, i) = p(1,2);
end
%% Correlation Brain signal var & INT by network
cd /home/qinggang/research/MIDUS_refresher/resting_conn/
%load('rsfmri_cov.mat');
load('rsfmri110_cov.mat');
cd community_new/
%cd community_new2/
load('GMtimes_sd_network.mat');
results = zeros(2, 5);
for i = 1:5
    [r, p] = corrcoef(rsfmri110cov(:,18), GMtimes_sd_network(:,i));
    results(1, i) = r(1,2);
    results(2, i) = p(1,2);
end
results_cov = zeros(2, 5);
for i = 1:5
    matrix = [rsfmri110cov(:,18) GMtimes_sd_network(:,i) rsfmri110cov(:, 2:8) rsfmri110cov(:, 14)];
    [rho, p] = partialcorr(matrix);
    results_cov(1, i) = rho(1,2);
    results_cov(2, i) = p(1,2);
end

%% Ploting Gephi
%select two subs
cd /home/qinggang/research/MIDUS_refresher/resting_conn/
load('/home/qinggang/research/MIDUS_refresher/resting_conn/power_e74_thresh50_30zed.mat')
load('/home/qinggang/research/MIDUS_refresher/resting_conn/community_new2/community_Louvain_243pos.mat')
M = power_e74_thresh50_30zed(:,:,[42 67]);
M_20 = zeros(243,243,2);
M_20(:,:,1) = weight_conversion(threshold_proportional(M(:,:,1), 0.2),'binarize');
M_20(:,:,2) = weight_conversion(threshold_proportional(M(:,:,2), 0.2),'binarize');
node_list = community_243(:, [1 6]);
edge_list1 = [];
edge_list2 = [];
curr_sub = M_20(:,:,1);
for i = 1:243
    for j = i:243
        if curr_sub(i, j) == 1
            edge_list1 = [edge_list1; i j];
        end
    end
end
curr_sub = M_20(:,:,2);
for i = 1:243
    for j = i:243
        if curr_sub(i, j) == 1
            edge_list2 = [edge_list2; i j];
        end
    end
end
%% Ploting Gephi
% 1SD above and below biomarker
cd /home/qinggang/research/MIDUS_refresher/resting_conn/
load('/home/qinggang/research/MIDUS_refresher/resting_conn/power_e74_thresh50_30zed.mat')
load('/home/qinggang/research/MIDUS_refresher/resting_conn/community_new2/community_Louvain_243pos.mat')
load('rsfmri_cov.mat')
M1 = power_e74_thresh50_30zed(:,:,rsfmricov(:,17) < -1);
M2 = power_e74_thresh50_30zed(:,:,rsfmricov(:,17) > 1);
M1_avg = mean(M1, 3);
M2_avg = mean(M2, 3);
M1_10 = weight_conversion(threshold_proportional(M1_avg, 0.1),'binarize');
M2_10 = weight_conversion(threshold_proportional(M2_avg, 0.1),'binarize');
edge_list1 = [];
edge_list2 = [];
for i = 1:243
    for j = i:243
        if M1_10(i, j) == 1
            edge_list1 = [edge_list1; i j];
        end
        if M2_10(i, j) == 1
            edge_list2 = [edge_list2; i j];
        end
    end
end
bm_low_10 = edge_list1;
bm_high_10 = edge_list2;

%% Ploting Gephi
% 1SD above and below modularity
cd /home/qinggang/research/MIDUS_refresher/resting_conn/
load('/home/qinggang/research/MIDUS_refresher/resting_conn/power_e74_thresh50_30zed.mat')
load('Qlouvain_243_bipos0110.mat')
M1 = power_e74_thresh50_30zed(:,:,Qlouvain_243_bipos0110(:,10) > 0.494);
M2 = power_e74_thresh50_30zed(:,:,Qlouvain_243_bipos0110(:,10) < 0.421);
M1_avg = mean(M1, 3);
M2_avg = mean(M2, 3);
M1_10 = weight_conversion(threshold_proportional(M1_avg, 0.1),'binarize');
M2_10 = weight_conversion(threshold_proportional(M2_avg, 0.1),'binarize');
edge_list1 = [];
edge_list2 = [];
for i = 1:243
    for j = i:243
        if M1_10(i, j) == 1
            edge_list1 = [edge_list1; i j];
        end
        if M2_10(i, j) == 1
            edge_list2 = [edge_list2; i j];
        end
    end
end
Q_high_10 = edge_list1;
Q_low_10 = edge_list2;
%% %% Ploting Gephi
% all subjects
cd /home/qinggang/research/MIDUS_refresher/resting_conn/
load('/home/qinggang/research/MIDUS_refresher/resting_conn/power_e74_thresh50_30zed.mat')
M_allsub = mean(power_e74_thresh50_30zed, 3);
M_allsub_20 = weight_conversion(threshold_proportional(M_allsub, 0.2),'binarize');
edge_list = [];
for i = 1:243
    for j = i:243
        if M_allsub_20(i, j) == 1
            edge_list = [edge_list; i j];
        end
    end
end
edge_allsub = edge_list;
%% Eigenvector centrality across 0.01 to 0.30 (243 nodes)

cd /home/qinggang/research/MIDUS_refresher/resting_conn
load('power_e74_thresh50_30zed.mat')
threshold = 0.01:0.01:0.30;
Eigen_243_bipos0130 = [];
for cost = threshold
    disp(cost)
    Eigen = [];
    for i = 1:size(power_e74_thresh50_30zed, 3)
        curr_sub = power_e74_thresh50_30zed(:, :, i);
        curr_sub_clean = weight_conversion(threshold_proportional(curr_sub, cost),'binarize');
        Eigen = [Eigen eigenvector_centrality_und(curr_sub_clean)];
        %Eigen = [Eigen; mean(efficiency_bin(curr_sub_clean, 1))];
    end
    Eigen_243_bipos0130 = cat(3, Eigen_243_bipos0130, Eigen);
    %Eigen_243_bipos0110 = [Eigen_243_bipos0110 Eloc];
end
%% Eigenvector centrality across 0.01 to 0.30 (208 nodes)

cd /home/qinggang/research/MIDUS_refresher/resting_conn
load('power_e110_thresh50_25zed.mat')
threshold = 0.01:0.01:0.30;
Eigen_208_bipos0130 = [];
for cost = threshold
    disp(cost)
    Eigen = [];
    for i = 1:size(power_e110_thresh50_25zed, 3)
        curr_sub = power_e110_thresh50_25zed(:, :, i);
        curr_sub_clean = weight_conversion(threshold_proportional(curr_sub, cost),'binarize');
        Eigen = [Eigen eigenvector_centrality_und(curr_sub_clean)];
        %Eigen = [Eigen; mean(efficiency_bin(curr_sub_clean, 1))];
    end
    Eigen_208_bipos0130 = cat(3, Eigen_208_bipos0130, Eigen);
    %Eigen_243_bipos0110 = [Eigen_243_bipos0110 Eloc];
end
%% Correlation between Eigenvector centrality and IND across all regions
cd /home/qinggang/research/MIDUS_refresher/resting_conn
load('Eigen_243_bipos0110.mat')
Eigen_avg = mean(Eigen_243_bipos0110(:,:,2:10), 3);
load('rsfmri_cov.mat')
results = zeros(2, 243);
for i = 1:243
    [r, p] = corrcoef(rsfmricov(:,19), Eigen_avg(i,:)');
    results(1, i) = r(1,2);
    results(2, i) = p(1,2);
end
results_cov = zeros(2, 243);
for i = 1:243
    matrix = [rsfmricov(:,19) Eigen_avg(i,:)' rsfmricov(:, 2:3) rsfmricov(:, 14)];
    [rho, p] = partialcorr(matrix);
    results_cov(1, i) = rho(1,2);
    results_cov(2, i) = p(1,2);
end
%% Correlation between Eigenvector centrality and IND by network (own assignment)
%cd /home/qinggang/research/MIDUS_refresher/resting_conn
cd /Users/Qinggang1/Documents/MATLAB/MR/MR_files_from_linux/
%load('Eigen_243_bipos0110.mat')
load('Eigen_243_bipos0130.mat')
%load('Eigen_208_bipos0130.mat')
%load('rsfmri110_cov.mat')
load('rsfmri_cov.mat')
cd community_new2/
%cd community_new/
%load('community_Louvain_208pos.mat')
load('community_Louvain_243pos.mat')
Eigen_avg = mean(Eigen_208_bipos0130(:,:,1:30), 3);
results = zeros(2, 5);
for i = 1:5
    Eigen_network = Eigen_avg(community_243(:,6) == i, :);
    Eigen_v = mean(Eigen_network, 1);
    [r, p] = corrcoef(rsfmricov(:,19), Eigen_v');
    results(1, i) = r(1,2);
    results(2, i) = p(1,2);
end
results_cov = zeros(2, 5);
for i = 1:5
    Eigen_network = Eigen_avg(community_243(:,6) == i, :);
    Eigen_v = mean(Eigen_network, 1);
    matrix = [rsfmricov(:,19) Eigen_v' rsfmricov(:, 2:8) rsfmricov(:, 14)];
    [rho, p] = partialcorr(matrix);
    results_cov(1, i) = rho(1,2);
    results_cov(2, i) = p(1,2);
end
%% Correlation between Eigenvector centrality and IND by network (BNV assignment)
cd /home/qinggang/research/MIDUS_refresher/resting_conn
load('Eigen_208_bipos0130.mat')
%load('Eigen_243_bipos0130.mat')
Eigen_avg = mean(Eigen_208_bipos0130(:,:,1:30), 3);
load('rsfmri110_cov.mat')
%load('rsfmri_cov.mat')
cd community_new/
%cd community_new2/
load('Power208_BNV_assignment.mat')
results = zeros(2, 5);
for i = 1:5
    Eigen_network = Eigen_avg(Power208_BrainNetViewer(:,6) == i, :);
    Eigen_v = mean(Eigen_network, 1);
    [r, p] = corrcoef(rsfmri110cov(:,19), Eigen_v');
    results(1, i) = r(1,2);
    results(2, i) = p(1,2);
end
results_cov = zeros(2, 5);
for i = 1:5
    Eigen_network = Eigen_avg(Power208_BrainNetViewer(:,6) == i, :);
    Eigen_v = mean(Eigen_network, 1);
    matrix = [rsfmri110cov(:,19) Eigen_v' rsfmri110cov(:, 2:8) rsfmri110cov(:, 14)];
    [rho, p] = partialcorr(matrix);
    results_cov(1, i) = rho(1,2);
    results_cov(2, i) = p(1,2);
end
%% Plotting eigenvector centrality
Eigen_allnode = mean(Eigen_208_bipos0130, 3);
Eigen_fpcn = Eigen_allnode(group_assign(:,6) == 3, :);
Eigen_fpcn = mean(Eigen_fpcn, 1);
Eigen_fpcn = Eigen_fpcn';

tbl = table(rsfmri110cov(:, 19), rsfmri110cov(:, 2), rsfmri110cov(:, 3), ...
    rsfmri110cov(:, 4), rsfmri110cov(:, 5), rsfmri110cov(:, 6), rsfmri110cov(:, 7), ... 
    rsfmri110cov(:, 8), rsfmri110cov(:,14), Eigen_fpcn, 'VariableNames', ... 
{'IND', 'sex', 'age', 'race1', 'race2', 'race3', 'race4', 'BMI', 'Education', 'fpcn_c'});

mdl = fitlm(tbl);

tbl1 = table(rsfmri110cov(:, 2), rsfmri110cov(:, 3), ...
    rsfmri110cov(:, 4), rsfmri110cov(:, 5), rsfmri110cov(:, 6), rsfmri110cov(:, 7), ... 
    rsfmri110cov(:, 8), rsfmri110cov(:,14), Eigen_fpcn, 'VariableNames', ... 
{'sex', 'age', 'race1', 'race2', 'race3', 'race4', 'BMI', 'Education', 'fpcn_c'});
mdl1 = fitlm(tbl1);
res1 = mdl1.Residuals.Raw;

tbl2 = table(rsfmri110cov(:, 2), rsfmri110cov(:, 3), ...
    rsfmri110cov(:, 4), rsfmri110cov(:, 5), rsfmri110cov(:, 6), rsfmri110cov(:, 7), ... 
    rsfmri110cov(:, 8), rsfmri110cov(:,14), rsfmri110cov(:, 19),  'VariableNames', ... 
{'sex', 'age', 'race1', 'race2', 'race3', 'race4', 'BMI', 'Education', 'IND'});
mdl2 = fitlm(tbl2);
res2 = mdl2.Residuals.Raw;
%% Eigen centrality of FPCN

Eigen_FPCN = mean(Eigen_243_bipos0130(community_243(:,6) == 1, :, :), 1);
Eigen_FPCN_avg = reshape(Eigen_FPCN, 74, 30);
