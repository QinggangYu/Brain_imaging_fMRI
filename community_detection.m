%% group partition
ciu_list = [];
Q_list = [];
for index = 1:size(power_clean_zeroed, 3)
	disp(index)
	current_sub = power_clean_zeroed(:,:,index);
    current_sub_clean = threshold_absolute(current_sub, 0);
    %current_sub_clean = weight_conversion(threshold_proportional(current_sub, 0.15), 'binarize');
    module = [];
    Q = [];
    for run = 1:500
		node = size(current_sub_clean, 1);
		M = 1:node;
		Q0 = -1;
		Q1 = 0;
		while Q1 - Q0 > 0.00001
			Q0 = Q1;
			[M, Q1] = community_louvain(current_sub_clean, [], M);
		end
		module = [module M];
        Q = [Q Q1];
    end
    D = agreement(module);
    D = D/500;
    ciu = consensus_und(D, 0.5, 500);
    ciu_list = [ciu_list ciu];
    Q_list = [Q_list; mean(Q)];
end
D_list = agreement(ciu_list);
D_list = D_list/106;
group_assign = consensus_und(D_list, 0.5, 500);
save('Qlouvain_219_posonly.mat', 'Q_list')
save('community_Louvain_219pos.mat', 'group_assign')