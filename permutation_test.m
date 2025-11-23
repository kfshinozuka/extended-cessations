function p_val = permutation_test(NS_cond, control_cond, n_perm)

num_epochs_NScond = length(NS_cond);
true_diff = mean(NS_cond) - mean(control_cond);
all_data = [NS_cond; control_cond];
nTotal = length(all_data);
perm_diffs = zeros(n_perm,1);
for i = 1:n_perm
    % Shuffle the labels
    shuffled_idx = randperm(nTotal);
    groupNS_idx = shuffled_idx(1:num_epochs_NScond);
    groupcond_idx = shuffled_idx(num_epochs_NScond+1:end);
    
    % Compute permuted means
    meanNS = mean(all_data(groupNS_idx));
    meancond = mean(all_data(groupcond_idx));
    
    % Store difference
    perm_diffs(i) = meanNS - meancond;
end

p_val = mean(abs(perm_diffs) >= abs(true_diff));