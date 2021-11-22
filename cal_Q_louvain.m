function Q = cal_Q_louvain(matrix, threshold, gamma)

%threshold = [0.11: 0.01: 0.20];
%threshold = 0.014;
%threshold = 0.019;
Q = [];
for index = 1:size(matrix, 3)
	disp(index)
	current_sub = matrix(:,:,index);
    %for i = 1:numel(current_sub)
        %if current_sub(i) < 0
            %current_sub(i) = current_sub(i) * (-1);
        %end
    %end
	Q_currsub = [];
	for cost = threshold
		current_sub_clean = weight_conversion(threshold_proportional(current_sub, cost),'binarize');
		Q_list = [];
		for run = 1:500
			node = size(current_sub_clean, 1);
			M = 1:node;
			Q0 = -1;
			Q1 = 0;
			while Q1 - Q0 > 0.00001
				Q0 = Q1;
				[M, Q1] = community_louvain(current_sub_clean, gamma, M);
			end
			Q_list = [Q_list Q1];
		end
		Q_currsub = [Q_currsub mean(Q_list)];
	end
	Q = [Q; Q_currsub];
end

