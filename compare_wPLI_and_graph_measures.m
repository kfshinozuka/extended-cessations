%% wPLI - Yeo networks - source-space

yeo_indices = {[(1:9)'; (51:58)'], [(10:15)'; (59:66)'], [(16:23)'; (67:73)'], [(24:30)'; (74:78)'], [(31:33)'; (79:80)'], [(34:37)'; (81:89)'], [(38:50)'; (90:100)']};
wPLI_yeo_epoch = cell(31,6);
for f = 1:31 % number of subjects
    for fr = 1:6 % number of frequency bands
        wPLI_yeo_epoch{f,fr} = zeros(7,7,size(wPLI_epoch{f},3));
        for ep = 1:size(wPLI_epoch{f},3)
            for i = 1:7
                for j = 1:7      
                    net_i = yeo_indices{i};
                    net_j = yeo_indices{j};  
                    subMat = wPLI_epoch{f,fr}(net_i, net_j, ep);
                    wPLI_yeo_epoch{f,fr}(i,j,ep) = mean(subMat(:));
                end
            end
        end
    end
end

%% Statistical analysis - source-space - Yeo - within-subject

n_perm = 5000;
pval_sub003 = ones(7,7,12); % network x network x (frequency band * condition)
pval_sub001 = ones(7,7,12);
pval_sub010 = ones(7,7,12);
pval_sub029 = ones(7,7,12);
pval_sub034 = ones(7,7,12);

for fr = 1:6
    for y1 = 1:7
        for y2 = y1:7
            sub001_NS = squeeze(wPLI_yeo_epoch{11,fr}(y1,y2,:));
            sub001_count = [squeeze(wPLI_yeo_epoch{7,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{9,fr}(y1,y2,:))];
            pval_sub001(y1,y2,(fr-1)*2+1) = permutation_test(sub001_NS, sub001_count, n_perm);
            pval_sub001(y2,y1,(fr-1)*2+1) = pval_sub001(y1,y2,(fr-1)*2+1);
            sub001_memory = [squeeze(wPLI_yeo_epoch{8,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{10,fr}(y1,y2,:))];
            pval_sub001(y1,y2,fr*2) = permutation_test(sub001_NS, sub001_memory, n_perm);
            pval_sub001(y2,y1,fr*2) = pval_sub001(y1,y2,fr*2);

            sub003_NS = [squeeze(wPLI_yeo_epoch{3,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{4,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{5,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{6,fr}(y1,y2,:))];
            sub003_count = squeeze(wPLI_yeo_epoch{1,fr}(y1,y2,:));
            pval_sub003(y1,y2,(fr-1)*2+1) = permutation_test(sub003_NS, sub003_count, n_perm);
            pval_sub003(y2,y1,(fr-1)*2+1) = pval_sub003(y1,y2,(fr-1)*2+1);
            sub003_memory = squeeze(wPLI_yeo_epoch{2,fr}(y1,y2,:));
            pval_sub003(y1,y2,fr*2) = permutation_test(sub003_NS, sub003_memory, n_perm);
            pval_sub003(y2,y1,fr*2) = pval_sub003(y1,y2,fr*2);
            
            sub010_NS = [squeeze(wPLI_yeo_epoch{17,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{18,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{19,fr}(y1,y2,:))];
            sub010_count = [squeeze(wPLI_yeo_epoch{12,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{14,fr}(y1,y2,:))];
            pval_sub010(y1,y2,(fr-1)*2+1) = permutation_test(sub010_NS, sub010_count, n_perm);
            pval_sub010(y2,y1,(fr-1)*2+1) = pval_sub010(y1,y2,(fr-1)*2+1);
            sub010_memory = [squeeze(wPLI_yeo_epoch{13,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{15,fr}(y1,y2,:))];
            pval_sub010(y1,y2,fr*2) = permutation_test(sub010_NS, sub010_memory, n_perm);
            pval_sub010(y2,y1,fr*2) = pval_sub010(y1,y2,fr*2);

            sub029_NS = [squeeze(wPLI_yeo_epoch{24,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{25,fr}(y1,y2,:))];
            sub029_count = [squeeze(wPLI_yeo_epoch{20,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{22,fr}(y1,y2,:))];
            pval_sub029(y1,y2,(fr-1)*2+1) = permutation_test(sub029_NS, sub029_count, n_perm);
            pval_sub029(y2,y1,(fr-1)*2+1) = pval_sub029(y1,y2,(fr-1)*2+1);
            sub029_memory = [squeeze(wPLI_yeo_epoch{21,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{23,fr}(y1,y2,:))];
            pval_sub029(y1,y2,fr*2) = permutation_test(sub029_NS, sub029_memory, n_perm);
            pval_sub029(y2,y1,fr*2) = pval_sub029(y1,y2,fr*2);  

            sub034_NS = [squeeze(wPLI_yeo_epoch{30,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{31,fr}(y1,y2,:))];
            sub034_count = [squeeze(wPLI_yeo_epoch{26,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{28,fr}(y1,y2,:))];
            pval_sub034(y1,y2,(fr-1)*2+1) = permutation_test(sub034_NS, sub034_count, n_perm);
            pval_sub034(y2,y1,(fr-1)*2+1) = pval_sub034(y1,y2,(fr-1)*2+1);
            sub034_memory = [squeeze(wPLI_yeo_epoch{27,fr}(y1,y2,:)); squeeze(wPLI_yeo_epoch{29,fr}(y1,y2,:))];
            pval_sub034(y1,y2,fr*2) = permutation_test(sub034_NS, sub034_memory, n_perm);
            pval_sub034(y2,y1,fr*2) = pval_sub034(y1,y2,fr*2);             
        end
    end
end

[~,~,~,adj_pval_sub001] = fdr_bh_2015(pval_sub001(:));
[~,~,~,adj_pval_sub003] = fdr_bh_2015(pval_sub003(:));
[~,~,~,adj_pval_sub010] = fdr_bh_2015(pval_sub010(:));
[~,~,~,adj_pval_sub029] = fdr_bh_2015(pval_sub029(:));
[~,~,~,adj_pval_sub034] = fdr_bh_2015(pval_sub034(:));

adj_pval_sub001 = reshape(adj_pval_sub001, size(pval_sub001));
adj_pval_sub003 = reshape(adj_pval_sub003, size(pval_sub003));
adj_pval_sub010 = reshape(adj_pval_sub010, size(pval_sub010));
adj_pval_sub029 = reshape(adj_pval_sub029, size(pval_sub029));
adj_pval_sub034 = reshape(adj_pval_sub034, size(pval_sub034));

A=adj_pval_sub001(:,:,1);
sub001_delta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub001(:,:,2);
sub001_delta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub003(:,:,1);
sub003_delta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub003(:,:,2);
sub003_delta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub010(:,:,1);
sub010_delta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub010(:,:,2);
sub010_delta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub029(:,:,1);
sub029_delta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub029(:,:,2);
sub029_delta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub034(:,:,1);
sub034_delta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub034(:,:,2);
sub034_delta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);

vectors_delta = {sub001_delta_count_sig, sub001_delta_memory_sig, sub003_delta_count_sig, sub003_delta_memory_sig, sub010_delta_count_sig, sub010_delta_memory_sig, sub029_delta_count_sig, sub029_delta_memory_sig, sub034_delta_count_sig, sub034_delta_memory_sig};

all_values_delta = vertcat(vectors_delta{:}); 
[unique_vals_delta, ~, idx] = unique(all_values_delta);
counts = accumarray(idx, 1);
[max_count, max_idx] = max(counts);
most_common_value_delta = unique_vals_delta(max_idx);

A=adj_pval_sub001(:,:,3);
sub001_theta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub001(:,:,4);
sub001_theta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub003(:,:,3);
sub003_theta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub003(:,:,4);
sub003_theta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub010(:,:,3);
sub010_theta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub010(:,:,4);
sub010_theta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub029(:,:,3);
sub029_theta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub029(:,:,4);
sub029_theta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub034(:,:,3);
sub034_theta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub034(:,:,4);
sub034_theta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);

vectors_theta = {sub001_theta_count_sig, sub001_theta_memory_sig, sub003_theta_count_sig, sub003_theta_memory_sig, sub010_theta_count_sig, sub010_theta_memory_sig, sub029_theta_count_sig, sub029_theta_memory_sig, sub034_theta_count_sig, sub034_theta_memory_sig};

all_values_theta = vertcat(vectors_theta{:});
[unique_vals_theta, ~, idx] = unique(all_values_theta);
counts = accumarray(idx, 1);
[max_count, max_idx] = max(counts);
most_common_value_theta = unique_vals_theta(max_idx);

A=adj_pval_sub001(:,:,5);
sub001_alpha_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub001(:,:,6);
sub001_alpha_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub003(:,:,5);
sub003_alpha_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub003(:,:,6);
sub003_alpha_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub010(:,:,5);
sub010_alpha_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub010(:,:,6);
sub010_alpha_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub029(:,:,5);
sub029_alpha_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub029(:,:,6);
sub029_alpha_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub034(:,:,5);
sub034_alpha_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub034(:,:,6);
sub034_alpha_memory_sig = find(A(triu(true(size(A)), 1))<0.05);

vectors_alpha = {sub001_alpha_count_sig, sub001_alpha_memory_sig, sub003_alpha_count_sig, sub003_alpha_memory_sig, sub010_alpha_count_sig, sub010_alpha_memory_sig, sub029_alpha_count_sig, sub029_alpha_memory_sig, sub034_alpha_count_sig, sub034_alpha_memory_sig};

all_values_alpha = vertcat(vectors_alpha{:});  
[unique_vals_alpha, ~, idx] = unique(all_values_alpha);
counts = accumarray(idx, 1);
[max_count, max_idx] = max(counts);
most_common_value_alpha = unique_vals_alpha(max_idx);

A=adj_pval_sub001(:,:,7);
sub001_lowbeta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub001(:,:,8);
sub001_lowbeta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub003(:,:,7);
sub003_lowbeta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub003(:,:,8);
sub003_lowbeta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub010(:,:,7);
sub010_lowbeta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub010(:,:,8);
sub010_lowbeta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub029(:,:,7);
sub029_lowbeta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub029(:,:,8);
sub029_lowbeta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub034(:,:,7);
sub034_lowbeta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub034(:,:,8);
sub034_lowbeta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);

vectors_lowbeta = {sub001_lowbeta_count_sig, sub001_lowbeta_memory_sig, sub003_lowbeta_count_sig, sub003_lowbeta_memory_sig, sub010_lowbeta_count_sig, sub010_lowbeta_memory_sig, sub029_lowbeta_count_sig, sub029_lowbeta_memory_sig, sub034_lowbeta_count_sig, sub034_lowbeta_memory_sig};

all_values_lowbeta = vertcat(vectors_lowbeta{:}); 
[unique_vals_lowbeta, ~, idx] = unique(all_values_lowbeta);
counts = accumarray(idx, 1);
[max_count, max_idx] = max(counts);
most_common_value_lowbeta = unique_vals_lowbeta(max_idx);

A=adj_pval_sub001(:,:,9);
sub001_highbeta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub001(:,:,10);
sub001_highbeta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub003(:,:,9);
sub003_highbeta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub003(:,:,10);
sub003_highbeta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub010(:,:,9);
sub010_highbeta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub010(:,:,10);
sub010_highbeta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub029(:,:,9);
sub029_highbeta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub029(:,:,10);
sub029_highbeta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub034(:,:,9);
sub034_highbeta_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub034(:,:,10);
sub034_highbeta_memory_sig = find(A(triu(true(size(A)), 1))<0.05);

vectors_highbeta = {sub001_highbeta_count_sig, sub001_highbeta_memory_sig, sub003_highbeta_count_sig, sub003_highbeta_memory_sig, sub010_highbeta_count_sig, sub010_highbeta_memory_sig, sub029_highbeta_count_sig, sub029_highbeta_memory_sig, sub034_highbeta_count_sig, sub034_highbeta_memory_sig};

all_values_highbeta = vertcat(vectors_highbeta{:});  
[unique_vals_highbeta, ~, idx] = unique(all_values_highbeta);
counts = accumarray(idx, 1);
[max_count, max_idx] = max(counts);
most_common_value_highbeta = unique_vals_highbeta(max_idx);

A=adj_pval_sub001(:,:,11);
sub001_gamma_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub001(:,:,12);
sub001_gamma_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub003(:,:,11);
sub003_gamma_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub003(:,:,12);
sub003_gamma_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub010(:,:,11);
sub010_gamma_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub010(:,:,12);
sub010_gamma_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub029(:,:,11);
sub029_gamma_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub029(:,:,12);
sub029_gamma_memory_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub034(:,:,11);
sub034_gamma_count_sig = find(A(triu(true(size(A)), 1))<0.05);
A=adj_pval_sub034(:,:,12);
sub034_gamma_memory_sig = find(A(triu(true(size(A)), 1))<0.05);

vectors_gamma = {sub001_gamma_count_sig, sub001_gamma_memory_sig, sub003_gamma_count_sig, sub003_gamma_memory_sig, sub010_gamma_count_sig, sub010_gamma_memory_sig, sub029_gamma_count_sig, sub029_gamma_memory_sig, sub034_gamma_count_sig, sub034_gamma_memory_sig};

all_values_gamma = vertcat(vectors_gamma{:}); 
[unique_vals_gamma, ~, idx] = unique(all_values_gamma);
counts = accumarray(idx, 1);
[max_count, max_idx] = max(counts);
most_common_value_gamma = unique_vals_gamma(max_idx);

%% wPLI plotting parameters

% Row/column indices where to draw lines
idx = [17, 31, 46, 58, 63, 76];

% Corresponding labels
labels = {'visual', 'SMN', 'DAN', 'VAN', 'limbic', ...
          'FPN', 'DMN'};

% Original LH and RH indices for each network
visual_LH = 1:9;
SMN_LH    = 10:15;
DAN_LH    = 16:23;
VAN_LH    = 24:30;
limbic_LH = 31:33;
FPN_LH    = 34:37;
DMN_LH    = 38:50;

visual_RH = 51:58;
SMN_RH    = 59:66;
DAN_RH    = 67:73;
VAN_RH    = 74:78;
limbic_RH = 79:80;
FPN_RH    = 81:89;
DMN_RH    = 90:100;

% New ordering: LH then RH for each network
new_order = [visual_LH, visual_RH, ...
             SMN_LH, SMN_RH, ...
             DAN_LH, DAN_RH, ...
             VAN_LH, VAN_RH, ...
             limbic_LH, limbic_RH, ...
             FPN_LH, FPN_RH, ...
             DMN_LH, DMN_RH];

%% Plot wPLI - Yeo - delta

% Custom tick labels
xlabels = {'Visual', 'SMN', 'DAN', 'VAN', 'Limbic', 'FPN', 'DMN'};
ylabels = xlabels;

wPLI_yeo_sub001_count_delta = {mean(wPLI_yeo_epoch{7,1},3), mean(wPLI_yeo_epoch{9,1},3)};
wPLI_yeo_sub001_count_delta_avg = zeros(size(wPLI_yeo_epoch{7,1},1));
for i = 1:size(wPLI_yeo_epoch{7,1},1)
    for j = 1:size(wPLI_yeo_epoch{7,1},1)
        wPLI_yeo_sub001_count_delta_avg(i,j) = mean([wPLI_yeo_sub001_count_delta{1}(i,j) wPLI_yeo_sub001_count_delta{2}(i,j)]);
    end
end

wPLI_yeo_sub001_memory_delta = {mean(wPLI_yeo_epoch{8,1},3), mean(wPLI_yeo_epoch{10,1},3)};
wPLI_yeo_sub001_memory_delta_avg = zeros(size(wPLI_yeo_epoch{8,1},1));
for i = 1:size(wPLI_yeo_epoch{8,1},1)
    for j = 1:size(wPLI_yeo_epoch{8,1},1)
        wPLI_yeo_sub001_memory_delta_avg(i,j) = mean([wPLI_yeo_sub001_memory_delta{1}(i,j) wPLI_yeo_sub001_memory_delta{2}(i,j)]);
    end
end

wPLI_yeo_sub001_NS_delta_avg = mean(wPLI_yeo_epoch{11,1},3);

min_val = min(min([wPLI_yeo_sub001_count_delta_avg(:); wPLI_yeo_sub001_memory_delta_avg(:); wPLI_yeo_sub001_NS_delta_avg(:)]));
max_val = max(max([wPLI_yeo_sub001_count_delta_avg(:); wPLI_yeo_sub001_memory_delta_avg(:); wPLI_yeo_sub001_NS_delta_avg(:)]));

figure; 
subplot(3,5,1)
imagesc(wPLI_yeo_sub001_count_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,6)
imagesc(wPLI_yeo_sub001_memory_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,11)
imagesc(wPLI_yeo_sub001_NS_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub003_count_delta_avg = mean(wPLI_yeo_epoch{1,1},3);

wPLI_yeo_sub003_memory_delta_avg = mean(wPLI_yeo_epoch{2,1},3);

wPLI_yeo_sub003_NS_delta = {mean(wPLI_yeo_epoch{3,1},3), mean(wPLI_yeo_epoch{4,1},3), mean(wPLI_yeo_epoch{5,1},3), mean(wPLI_yeo_epoch{6,1},3)};
wPLI_yeo_sub003_NS_delta_avg = zeros(size(wPLI_yeo_epoch{3,1},1));
for i = 1:size(wPLI_yeo_epoch{3,1},1)
    for j = 1:size(wPLI_yeo_epoch{3,1},1)
        wPLI_yeo_sub003_NS_delta_avg(i,j) = mean([wPLI_yeo_sub003_NS_delta{1}(i,j) wPLI_yeo_sub003_NS_delta{2}(i,j) wPLI_yeo_sub003_NS_delta{3}(i,j) wPLI_yeo_sub003_NS_delta{4}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub003_count_delta_avg; wPLI_yeo_sub003_memory_delta_avg; wPLI_yeo_sub003_NS_delta_avg]));
max_val = max(max([wPLI_yeo_sub003_count_delta_avg; wPLI_yeo_sub003_memory_delta_avg; wPLI_yeo_sub003_NS_delta_avg]));

subplot(3,5,2)
imagesc(wPLI_yeo_sub003_count_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,7)
imagesc(wPLI_yeo_sub003_memory_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,12)
imagesc(wPLI_yeo_sub003_NS_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub010_count_delta = {mean(wPLI_yeo_epoch{12,1},3), mean(wPLI_yeo_epoch{14,1},3)};
wPLI_yeo_sub010_count_delta_avg = zeros(size(wPLI_yeo_epoch{12,1},1));
for i = 1:size(wPLI_yeo_epoch{12,1},1)
    for j = 1:size(wPLI_yeo_epoch{12,1},1)
        wPLI_yeo_sub010_count_delta_avg(i,j) = mean([wPLI_yeo_sub010_count_delta{1}(i,j) wPLI_yeo_sub010_count_delta{2}(i,j)]);
    end
end

wPLI_yeo_sub010_memory_delta = {mean(wPLI_yeo_epoch{13,1},3), mean(wPLI_yeo_epoch{15,1},3)};
wPLI_yeo_sub010_memory_delta_avg = zeros(size(wPLI_yeo_epoch{13,1},1));
for i = 1:size(wPLI_yeo_epoch{13,1},1)
    for j = 1:size(wPLI_yeo_epoch{13,1},1)
        wPLI_yeo_sub010_memory_delta_avg(i,j) = mean([wPLI_yeo_sub010_memory_delta{1}(i,j) wPLI_yeo_sub010_memory_delta{2}(i,j)]);
    end
end

wPLI_yeo_sub010_NS_delta = {mean(wPLI_yeo_epoch{17,1},3), mean(wPLI_yeo_epoch{18,1},3), mean(wPLI_yeo_epoch{19,1},3)};
wPLI_yeo_sub010_NS_delta_avg = zeros(size(wPLI_yeo_epoch{17,1},1));
for i = 1:size(wPLI_yeo_epoch{17,1},1)
    for j = 1:size(wPLI_yeo_epoch{17,1},1)
        wPLI_yeo_sub010_NS_delta_avg(i,j) = mean([wPLI_yeo_sub010_NS_delta{1}(i,j) wPLI_yeo_sub010_NS_delta{2}(i,j) wPLI_yeo_sub010_NS_delta{3}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub010_count_delta_avg(:); wPLI_yeo_sub010_memory_delta_avg(:); wPLI_yeo_sub010_NS_delta_avg(:)]));
max_val = max(max([wPLI_yeo_sub010_count_delta_avg(:); wPLI_yeo_sub010_memory_delta_avg(:); wPLI_yeo_sub010_NS_delta_avg(:)]));

subplot(3,5,3)
imagesc(wPLI_yeo_sub010_count_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,8)
imagesc(wPLI_yeo_sub010_memory_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,13)
imagesc(wPLI_yeo_sub010_NS_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub029_count_delta = {mean(wPLI_yeo_epoch{20,1},3), mean(wPLI_yeo_epoch{22,1},3)};
wPLI_yeo_sub029_count_delta_avg = zeros(size(wPLI_yeo_epoch{20,1},1));
for i = 1:size(wPLI_yeo_epoch{20,1},1)
    for j = 1:size(wPLI_yeo_epoch{20,1},1)
        wPLI_yeo_sub029_count_delta_avg(i,j) = mean([wPLI_yeo_sub029_count_delta{1}(i,j) wPLI_yeo_sub029_count_delta{2}(i,j)]);
    end
end

wPLI_yeo_sub029_memory_delta = {mean(wPLI_yeo_epoch{21,1},3), mean(wPLI_yeo_epoch{23,1},3)};
wPLI_yeo_sub029_memory_delta_avg = zeros(size(wPLI_yeo_epoch{21,1},1));
for i = 1:size(wPLI_yeo_epoch{21,1},1)
    for j = 1:size(wPLI_yeo_epoch{21,1},1)
        wPLI_yeo_sub029_memory_delta_avg(i,j) = mean([wPLI_yeo_sub029_memory_delta{1}(i,j) wPLI_yeo_sub029_memory_delta{2}(i,j)]);
    end
end

wPLI_yeo_sub029_NS_delta = {mean(wPLI_yeo_epoch{24,1},3), mean(wPLI_yeo_epoch{25,1},3)};
wPLI_yeo_sub029_NS_delta_avg = zeros(size(wPLI_yeo_epoch{24,1},1));
for i = 1:size(wPLI_yeo_epoch{24,1},1)
    for j = 1:size(wPLI_yeo_epoch{24,1},1)
        wPLI_yeo_sub029_NS_delta_avg(i,j) = mean([wPLI_yeo_sub029_NS_delta{1}(i,j) wPLI_yeo_sub029_NS_delta{2}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub029_count_delta_avg(:); wPLI_yeo_sub029_memory_delta_avg(:); wPLI_yeo_sub029_NS_delta_avg(:)]));
max_val = max(max([wPLI_yeo_sub029_count_delta_avg(:); wPLI_yeo_sub029_memory_delta_avg(:); wPLI_yeo_sub029_NS_delta_avg(:)]));

subplot(3,5,4)
imagesc(wPLI_yeo_sub029_count_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,9)
imagesc(wPLI_yeo_sub029_memory_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,14)
imagesc(wPLI_yeo_sub029_NS_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub034_count_delta = {mean(wPLI_yeo_epoch{26,1},3), mean(wPLI_yeo_epoch{28,1},3)};
wPLI_yeo_sub034_count_delta_avg = zeros(size(wPLI_yeo_epoch{7,1},1));
for i = 1:size(wPLI_yeo_epoch{26,1},1)
    for j = 1:size(wPLI_yeo_epoch{26,1},1)
        wPLI_yeo_sub034_count_delta_avg(i,j) = mean([wPLI_yeo_sub034_count_delta{1}(i,j) wPLI_yeo_sub034_count_delta{2}(i,j)]);
    end
end

wPLI_yeo_sub034_memory_delta = {mean(wPLI_yeo_epoch{27,1},3), mean(wPLI_yeo_epoch{29,1},3)};
wPLI_yeo_sub034_memory_delta_avg = zeros(size(wPLI_yeo_epoch{8,1},1));
for i = 1:size(wPLI_yeo_epoch{27,1},1)
    for j = 1:size(wPLI_yeo_epoch{27,1},1)
        wPLI_yeo_sub034_memory_delta_avg(i,j) = mean([wPLI_yeo_sub034_memory_delta{1}(i,j) wPLI_yeo_sub034_memory_delta{2}(i,j)]);
    end
end

wPLI_yeo_sub034_NS_delta = {mean(wPLI_yeo_epoch{30,1},3), mean(wPLI_yeo_epoch{31,1},3)};
wPLI_yeo_sub034_NS_delta_avg = zeros(size(wPLI_yeo_epoch{30,1},1));
for i = 1:size(wPLI_yeo_epoch{30,1},1)
    for j = 1:size(wPLI_yeo_epoch{30,1},1)
        wPLI_yeo_sub034_NS_delta_avg(i,j) = mean([wPLI_yeo_sub034_NS_delta{1}(i,j) wPLI_yeo_sub034_NS_delta{2}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub034_count_delta_avg(:); wPLI_yeo_sub034_memory_delta_avg(:); wPLI_yeo_sub034_NS_delta_avg(:)]));
max_val = max(max([wPLI_yeo_sub034_count_delta_avg(:); wPLI_yeo_sub034_memory_delta_avg(:); wPLI_yeo_sub034_NS_delta_avg(:)]));

subplot(3,5,5)
imagesc(wPLI_yeo_sub034_count_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,10)
imagesc(wPLI_yeo_sub034_memory_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,15)
imagesc(wPLI_yeo_sub034_NS_delta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

%% Plot wPLI - yeo - theta

wPLI_yeo_sub001_count_theta = {mean(wPLI_yeo_epoch{7,2},3), mean(wPLI_yeo_epoch{9,2},3)};
wPLI_yeo_sub001_count_theta_avg = zeros(size(wPLI_yeo_epoch{7,2},1));
for i = 1:size(wPLI_yeo_epoch{7,2},1)
    for j = 1:size(wPLI_yeo_epoch{7,2},1)
        wPLI_yeo_sub001_count_theta_avg(i,j) = mean([wPLI_yeo_sub001_count_theta{1}(i,j) wPLI_yeo_sub001_count_theta{2}(i,j)]);
    end
end

wPLI_yeo_sub001_memory_theta = {mean(wPLI_yeo_epoch{8,2},3), mean(wPLI_yeo_epoch{10,2},3)};
wPLI_yeo_sub001_memory_theta_avg = zeros(size(wPLI_yeo_epoch{8,2},1));
for i = 1:size(wPLI_yeo_epoch{8,2},1)
    for j = 1:size(wPLI_yeo_epoch{8,2},1)
        wPLI_yeo_sub001_memory_theta_avg(i,j) = mean([wPLI_yeo_sub001_memory_theta{1}(i,j) wPLI_yeo_sub001_memory_theta{2}(i,j)]);
    end
end

wPLI_yeo_sub001_NS_theta_avg = mean(wPLI_yeo_epoch{11,2},3);

min_val = min(min([wPLI_yeo_sub001_count_theta_avg(:); wPLI_yeo_sub001_memory_theta_avg(:); wPLI_yeo_sub001_NS_theta_avg(:)]));
max_val = max(max([wPLI_yeo_sub001_count_theta_avg(:); wPLI_yeo_sub001_memory_theta_avg(:); wPLI_yeo_sub001_NS_theta_avg(:)]));

figure; 
subplot(3,5,1)
imagesc(wPLI_yeo_sub001_count_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,6)
imagesc(wPLI_yeo_sub001_memory_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,11)
imagesc(wPLI_yeo_sub001_NS_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub003_count_theta_avg = mean(wPLI_yeo_epoch{1,2},3);

wPLI_yeo_sub003_memory_theta_avg = mean(wPLI_yeo_epoch{2,2},3);

wPLI_yeo_sub003_NS_theta = {mean(wPLI_yeo_epoch{3,2},3), mean(wPLI_yeo_epoch{4,2},3), mean(wPLI_yeo_epoch{5,2},3), mean(wPLI_yeo_epoch{6,2},3)};
wPLI_yeo_sub003_NS_theta_avg = zeros(size(wPLI_yeo_epoch{3,2},1));
for i = 1:size(wPLI_yeo_epoch{3,2},1)
    for j = 1:size(wPLI_yeo_epoch{3,2},1)
        wPLI_yeo_sub003_NS_theta_avg(i,j) = mean([wPLI_yeo_sub003_NS_theta{1}(i,j) wPLI_yeo_sub003_NS_theta{2}(i,j) wPLI_yeo_sub003_NS_theta{3}(i,j) wPLI_yeo_sub003_NS_theta{4}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub003_count_theta_avg; wPLI_yeo_sub003_memory_theta_avg; wPLI_yeo_sub003_NS_theta_avg]));
max_val = max(max([wPLI_yeo_sub003_count_theta_avg; wPLI_yeo_sub003_memory_theta_avg; wPLI_yeo_sub003_NS_theta_avg]));

subplot(3,5,2)
imagesc(wPLI_yeo_sub003_count_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,7)
imagesc(wPLI_yeo_sub003_memory_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,12)
imagesc(wPLI_yeo_sub003_NS_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub010_count_theta = {mean(wPLI_yeo_epoch{12,2},3), mean(wPLI_yeo_epoch{14,2},3)};
wPLI_yeo_sub010_count_theta_avg = zeros(size(wPLI_yeo_epoch{12,2},1));
for i = 1:size(wPLI_yeo_epoch{12,2},1)
    for j = 1:size(wPLI_yeo_epoch{12,2},1)
        wPLI_yeo_sub010_count_theta_avg(i,j) = mean([wPLI_yeo_sub010_count_theta{1}(i,j) wPLI_yeo_sub010_count_theta{2}(i,j)]);
    end
end

wPLI_yeo_sub010_memory_theta = {mean(wPLI_yeo_epoch{13,2},3), mean(wPLI_yeo_epoch{15,2},3)};
wPLI_yeo_sub010_memory_theta_avg = zeros(size(wPLI_yeo_epoch{13,2},1));
for i = 1:size(wPLI_yeo_epoch{13,2},1)
    for j = 1:size(wPLI_yeo_epoch{13,2},1)
        wPLI_yeo_sub010_memory_theta_avg(i,j) = mean([wPLI_yeo_sub010_memory_theta{1}(i,j) wPLI_yeo_sub010_memory_theta{2}(i,j)]);
    end
end

wPLI_yeo_sub010_NS_theta = {mean(wPLI_yeo_epoch{17,2},3), mean(wPLI_yeo_epoch{18,2},3), mean(wPLI_yeo_epoch{19,2},3)};
wPLI_yeo_sub010_NS_theta_avg = zeros(size(wPLI_yeo_epoch{17,2},1));
for i = 1:size(wPLI_yeo_epoch{17,2},1)
    for j = 1:size(wPLI_yeo_epoch{17,2},1)
        wPLI_yeo_sub010_NS_theta_avg(i,j) = mean([wPLI_yeo_sub010_NS_theta{1}(i,j) wPLI_yeo_sub010_NS_theta{2}(i,j) wPLI_yeo_sub010_NS_theta{3}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub010_count_theta_avg(:); wPLI_yeo_sub010_memory_theta_avg(:); wPLI_yeo_sub010_NS_theta_avg(:)]));
max_val = max(max([wPLI_yeo_sub010_count_theta_avg(:); wPLI_yeo_sub010_memory_theta_avg(:); wPLI_yeo_sub010_NS_theta_avg(:)]));

subplot(3,5,3)
imagesc(wPLI_yeo_sub010_count_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,8)
imagesc(wPLI_yeo_sub010_memory_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,13)
imagesc(wPLI_yeo_sub010_NS_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub029_count_theta = {mean(wPLI_yeo_epoch{20,2},3), mean(wPLI_yeo_epoch{22,2},3)};
wPLI_yeo_sub029_count_theta_avg = zeros(size(wPLI_yeo_epoch{20,2},1));
for i = 1:size(wPLI_yeo_epoch{20,2},1)
    for j = 1:size(wPLI_yeo_epoch{20,2},1)
        wPLI_yeo_sub029_count_theta_avg(i,j) = mean([wPLI_yeo_sub029_count_theta{1}(i,j) wPLI_yeo_sub029_count_theta{2}(i,j)]);
    end
end

wPLI_yeo_sub029_memory_theta = {mean(wPLI_yeo_epoch{21,2},3), mean(wPLI_yeo_epoch{23,2},3)};
wPLI_yeo_sub029_memory_theta_avg = zeros(size(wPLI_yeo_epoch{21,2},1));
for i = 1:size(wPLI_yeo_epoch{21,2},1)
    for j = 1:size(wPLI_yeo_epoch{21,2},1)
        wPLI_yeo_sub029_memory_theta_avg(i,j) = mean([wPLI_yeo_sub029_memory_theta{1}(i,j) wPLI_yeo_sub029_memory_theta{2}(i,j)]);
    end
end

wPLI_yeo_sub029_NS_theta = {mean(wPLI_yeo_epoch{24,2},3), mean(wPLI_yeo_epoch{25,2},3)};
wPLI_yeo_sub029_NS_theta_avg = zeros(size(wPLI_yeo_epoch{24,2},1));
for i = 1:size(wPLI_yeo_epoch{24,2},1)
    for j = 1:size(wPLI_yeo_epoch{24,2},1)
        wPLI_yeo_sub029_NS_theta_avg(i,j) = mean([wPLI_yeo_sub029_NS_theta{1}(i,j) wPLI_yeo_sub029_NS_theta{2}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub029_count_theta_avg(:); wPLI_yeo_sub029_memory_theta_avg(:); wPLI_yeo_sub029_NS_theta_avg(:)]));
max_val = max(max([wPLI_yeo_sub029_count_theta_avg(:); wPLI_yeo_sub029_memory_theta_avg(:); wPLI_yeo_sub029_NS_theta_avg(:)]));

subplot(3,5,4)
imagesc(wPLI_yeo_sub029_count_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,9)
imagesc(wPLI_yeo_sub029_memory_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,14)
imagesc(wPLI_yeo_sub029_NS_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub034_count_theta = {mean(wPLI_yeo_epoch{26,2},3), mean(wPLI_yeo_epoch{28,2},3)};
wPLI_yeo_sub034_count_theta_avg = zeros(size(wPLI_yeo_epoch{7,2},1));
for i = 1:size(wPLI_yeo_epoch{26,2},1)
    for j = 1:size(wPLI_yeo_epoch{26,2},1)
        wPLI_yeo_sub034_count_theta_avg(i,j) = mean([wPLI_yeo_sub034_count_theta{1}(i,j) wPLI_yeo_sub034_count_theta{2}(i,j)]);
    end
end

wPLI_yeo_sub034_memory_theta = {mean(wPLI_yeo_epoch{27,2},3), mean(wPLI_yeo_epoch{29,2},3)};
wPLI_yeo_sub034_memory_theta_avg = zeros(size(wPLI_yeo_epoch{8,2},1));
for i = 1:size(wPLI_yeo_epoch{27,2},1)
    for j = 1:size(wPLI_yeo_epoch{27,2},1)
        wPLI_yeo_sub034_memory_theta_avg(i,j) = mean([wPLI_yeo_sub034_memory_theta{1}(i,j) wPLI_yeo_sub034_memory_theta{2}(i,j)]);
    end
end

wPLI_yeo_sub034_NS_theta = {mean(wPLI_yeo_epoch{30,2},3), mean(wPLI_yeo_epoch{31,2},3)};
wPLI_yeo_sub034_NS_theta_avg = zeros(size(wPLI_yeo_epoch{30,2},1));
for i = 1:size(wPLI_yeo_epoch{30,2},1)
    for j = 1:size(wPLI_yeo_epoch{30,2},1)
        wPLI_yeo_sub034_NS_theta_avg(i,j) = mean([wPLI_yeo_sub034_NS_theta{1}(i,j) wPLI_yeo_sub034_NS_theta{2}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub034_count_theta_avg(:); wPLI_yeo_sub034_memory_theta_avg(:); wPLI_yeo_sub034_NS_theta_avg(:)]));
max_val = max(max([wPLI_yeo_sub034_count_theta_avg(:); wPLI_yeo_sub034_memory_theta_avg(:); wPLI_yeo_sub034_NS_theta_avg(:)]));

subplot(3,5,5)
imagesc(wPLI_yeo_sub034_count_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,10)
imagesc(wPLI_yeo_sub034_memory_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,15)
imagesc(wPLI_yeo_sub034_NS_theta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

%% Plot wPLI - yeo - alpha

wPLI_yeo_sub001_count_alpha = {mean(wPLI_yeo_epoch{7,3},3), mean(wPLI_yeo_epoch{9,3},3)};
wPLI_yeo_sub001_count_alpha_avg = zeros(size(wPLI_yeo_epoch{7,3},1));
for i = 1:size(wPLI_yeo_epoch{7,3},1)
    for j = 1:size(wPLI_yeo_epoch{7,3},1)
        wPLI_yeo_sub001_count_alpha_avg(i,j) = mean([wPLI_yeo_sub001_count_alpha{1}(i,j) wPLI_yeo_sub001_count_alpha{2}(i,j)]);
    end
end

wPLI_yeo_sub001_memory_alpha = {mean(wPLI_yeo_epoch{8,3},3), mean(wPLI_yeo_epoch{10,3},3)};
wPLI_yeo_sub001_memory_alpha_avg = zeros(size(wPLI_yeo_epoch{8,3},1));
for i = 1:size(wPLI_yeo_epoch{8,3},1)
    for j = 1:size(wPLI_yeo_epoch{8,3},1)
        wPLI_yeo_sub001_memory_alpha_avg(i,j) = mean([wPLI_yeo_sub001_memory_alpha{1}(i,j) wPLI_yeo_sub001_memory_alpha{2}(i,j)]);
    end
end

wPLI_yeo_sub001_NS_alpha_avg = mean(wPLI_yeo_epoch{11,3},3);

min_val = min(min([wPLI_yeo_sub001_count_alpha_avg(:); wPLI_yeo_sub001_memory_alpha_avg(:); wPLI_yeo_sub001_NS_alpha_avg(:)]));
max_val = max(max([wPLI_yeo_sub001_count_alpha_avg(:); wPLI_yeo_sub001_memory_alpha_avg(:); wPLI_yeo_sub001_NS_alpha_avg(:)]));

figure; 
subplot(3,5,1)
imagesc(wPLI_yeo_sub001_count_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,6)
imagesc(wPLI_yeo_sub001_memory_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,11)
imagesc(wPLI_yeo_sub001_NS_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub003_count_alpha_avg = mean(wPLI_yeo_epoch{1,3},3);

wPLI_yeo_sub003_memory_alpha_avg = mean(wPLI_yeo_epoch{2,3},3);

wPLI_yeo_sub003_NS_alpha = {mean(wPLI_yeo_epoch{3,3},3), mean(wPLI_yeo_epoch{4,3},3), mean(wPLI_yeo_epoch{5,3},3), mean(wPLI_yeo_epoch{6,3},3)};
wPLI_yeo_sub003_NS_alpha_avg = zeros(size(wPLI_yeo_epoch{3,3},1));
for i = 1:size(wPLI_yeo_epoch{3,3},1)
    for j = 1:size(wPLI_yeo_epoch{3,3},1)
        wPLI_yeo_sub003_NS_alpha_avg(i,j) = mean([wPLI_yeo_sub003_NS_alpha{1}(i,j) wPLI_yeo_sub003_NS_alpha{2}(i,j) wPLI_yeo_sub003_NS_alpha{3}(i,j) wPLI_yeo_sub003_NS_alpha{4}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub003_count_alpha_avg; wPLI_yeo_sub003_memory_alpha_avg; wPLI_yeo_sub003_NS_alpha_avg]));
max_val = max(max([wPLI_yeo_sub003_count_alpha_avg; wPLI_yeo_sub003_memory_alpha_avg; wPLI_yeo_sub003_NS_alpha_avg]));

subplot(3,5,2)
imagesc(wPLI_yeo_sub003_count_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,7)
imagesc(wPLI_yeo_sub003_memory_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,12)
imagesc(wPLI_yeo_sub003_NS_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub010_count_alpha = {mean(wPLI_yeo_epoch{12,3},3), mean(wPLI_yeo_epoch{14,3},3)};
wPLI_yeo_sub010_count_alpha_avg = zeros(size(wPLI_yeo_epoch{12,3},1));
for i = 1:size(wPLI_yeo_epoch{12,3},1)
    for j = 1:size(wPLI_yeo_epoch{12,3},1)
        wPLI_yeo_sub010_count_alpha_avg(i,j) = mean([wPLI_yeo_sub010_count_alpha{1}(i,j) wPLI_yeo_sub010_count_alpha{2}(i,j)]);
    end
end

wPLI_yeo_sub010_memory_alpha = {mean(wPLI_yeo_epoch{13,3},3), mean(wPLI_yeo_epoch{15,3},3)};
wPLI_yeo_sub010_memory_alpha_avg = zeros(size(wPLI_yeo_epoch{13,3},1));
for i = 1:size(wPLI_yeo_epoch{13,3},1)
    for j = 1:size(wPLI_yeo_epoch{13,3},1)
        wPLI_yeo_sub010_memory_alpha_avg(i,j) = mean([wPLI_yeo_sub010_memory_alpha{1}(i,j) wPLI_yeo_sub010_memory_alpha{2}(i,j)]);
    end
end

wPLI_yeo_sub010_NS_alpha = {mean(wPLI_yeo_epoch{17,3},3), mean(wPLI_yeo_epoch{18,3},3), mean(wPLI_yeo_epoch{19,3},3)};
wPLI_yeo_sub010_NS_alpha_avg = zeros(size(wPLI_yeo_epoch{17,3},1));
for i = 1:size(wPLI_yeo_epoch{17,3},1)
    for j = 1:size(wPLI_yeo_epoch{17,3},1)
        wPLI_yeo_sub010_NS_alpha_avg(i,j) = mean([wPLI_yeo_sub010_NS_alpha{1}(i,j) wPLI_yeo_sub010_NS_alpha{2}(i,j) wPLI_yeo_sub010_NS_alpha{3}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub010_count_alpha_avg(:); wPLI_yeo_sub010_memory_alpha_avg(:); wPLI_yeo_sub010_NS_alpha_avg(:)]));
max_val = max(max([wPLI_yeo_sub010_count_alpha_avg(:); wPLI_yeo_sub010_memory_alpha_avg(:); wPLI_yeo_sub010_NS_alpha_avg(:)]));

subplot(3,5,3)
imagesc(wPLI_yeo_sub010_count_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,8)
imagesc(wPLI_yeo_sub010_memory_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,13)
imagesc(wPLI_yeo_sub010_NS_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub029_count_alpha = {mean(wPLI_yeo_epoch{20,3},3), mean(wPLI_yeo_epoch{22,3},3)};
wPLI_yeo_sub029_count_alpha_avg = zeros(size(wPLI_yeo_epoch{20,3},1));
for i = 1:size(wPLI_yeo_epoch{20,3},1)
    for j = 1:size(wPLI_yeo_epoch{20,3},1)
        wPLI_yeo_sub029_count_alpha_avg(i,j) = mean([wPLI_yeo_sub029_count_alpha{1}(i,j) wPLI_yeo_sub029_count_alpha{2}(i,j)]);
    end
end

wPLI_yeo_sub029_memory_alpha = {mean(wPLI_yeo_epoch{21,3},3), mean(wPLI_yeo_epoch{23,3},3)};
wPLI_yeo_sub029_memory_alpha_avg = zeros(size(wPLI_yeo_epoch{21,3},1));
for i = 1:size(wPLI_yeo_epoch{21,3},1)
    for j = 1:size(wPLI_yeo_epoch{21,3},1)
        wPLI_yeo_sub029_memory_alpha_avg(i,j) = mean([wPLI_yeo_sub029_memory_alpha{1}(i,j) wPLI_yeo_sub029_memory_alpha{2}(i,j)]);
    end
end

wPLI_yeo_sub029_NS_alpha = {mean(wPLI_yeo_epoch{24,3},3), mean(wPLI_yeo_epoch{25,3},3)};
wPLI_yeo_sub029_NS_alpha_avg = zeros(size(wPLI_yeo_epoch{24,3},1));
for i = 1:size(wPLI_yeo_epoch{24,3},1)
    for j = 1:size(wPLI_yeo_epoch{24,3},1)
        wPLI_yeo_sub029_NS_alpha_avg(i,j) = mean([wPLI_yeo_sub029_NS_alpha{1}(i,j) wPLI_yeo_sub029_NS_alpha{2}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub029_count_alpha_avg(:); wPLI_yeo_sub029_memory_alpha_avg(:); wPLI_yeo_sub029_NS_alpha_avg(:)]));
max_val = max(max([wPLI_yeo_sub029_count_alpha_avg(:); wPLI_yeo_sub029_memory_alpha_avg(:); wPLI_yeo_sub029_NS_alpha_avg(:)]));

subplot(3,5,4)
imagesc(wPLI_yeo_sub029_count_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,9)
imagesc(wPLI_yeo_sub029_memory_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,14)
imagesc(wPLI_yeo_sub029_NS_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub034_count_alpha = {mean(wPLI_yeo_epoch{26,3},3), mean(wPLI_yeo_epoch{28,3},3)};
wPLI_yeo_sub034_count_alpha_avg = zeros(size(wPLI_yeo_epoch{7,3},1));
for i = 1:size(wPLI_yeo_epoch{26,3},1)
    for j = 1:size(wPLI_yeo_epoch{26,3},1)
        wPLI_yeo_sub034_count_alpha_avg(i,j) = mean([wPLI_yeo_sub034_count_alpha{1}(i,j) wPLI_yeo_sub034_count_alpha{2}(i,j)]);
    end
end

wPLI_yeo_sub034_memory_alpha = {mean(wPLI_yeo_epoch{27,3},3), mean(wPLI_yeo_epoch{29,3},3)};
wPLI_yeo_sub034_memory_alpha_avg = zeros(size(wPLI_yeo_epoch{8,3},1));
for i = 1:size(wPLI_yeo_epoch{27,3},1)
    for j = 1:size(wPLI_yeo_epoch{27,3},1)
        wPLI_yeo_sub034_memory_alpha_avg(i,j) = mean([wPLI_yeo_sub034_memory_alpha{1}(i,j) wPLI_yeo_sub034_memory_alpha{2}(i,j)]);
    end
end

wPLI_yeo_sub034_NS_alpha = {mean(wPLI_yeo_epoch{30,3},3), mean(wPLI_yeo_epoch{31,3},3)};
wPLI_yeo_sub034_NS_alpha_avg = zeros(size(wPLI_yeo_epoch{30,3},1));
for i = 1:size(wPLI_yeo_epoch{30,3},1)
    for j = 1:size(wPLI_yeo_epoch{30,3},1)
        wPLI_yeo_sub034_NS_alpha_avg(i,j) = mean([wPLI_yeo_sub034_NS_alpha{1}(i,j) wPLI_yeo_sub034_NS_alpha{2}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub034_count_alpha_avg(:); wPLI_yeo_sub034_memory_alpha_avg(:); wPLI_yeo_sub034_NS_alpha_avg(:)]));
max_val = max(max([wPLI_yeo_sub034_count_alpha_avg(:); wPLI_yeo_sub034_memory_alpha_avg(:); wPLI_yeo_sub034_NS_alpha_avg(:)]));

subplot(3,5,5)
imagesc(wPLI_yeo_sub034_count_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,10)
imagesc(wPLI_yeo_sub034_memory_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,15)
imagesc(wPLI_yeo_sub034_NS_alpha_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

%% Plot wPLI - yeo - low beta

wPLI_yeo_sub001_count_lowbeta = {mean(wPLI_yeo_epoch{7,4},3), mean(wPLI_yeo_epoch{9,4},3)};
wPLI_yeo_sub001_count_lowbeta_avg = zeros(size(wPLI_yeo_epoch{7,4},1));
for i = 1:size(wPLI_yeo_epoch{7,4},1)
    for j = 1:size(wPLI_yeo_epoch{7,4},1)
        wPLI_yeo_sub001_count_lowbeta_avg(i,j) = mean([wPLI_yeo_sub001_count_lowbeta{1}(i,j) wPLI_yeo_sub001_count_lowbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub001_memory_lowbeta = {mean(wPLI_yeo_epoch{8,4},3), mean(wPLI_yeo_epoch{10,4},3)};
wPLI_yeo_sub001_memory_lowbeta_avg = zeros(size(wPLI_yeo_epoch{8,4},1));
for i = 1:size(wPLI_yeo_epoch{8,4},1)
    for j = 1:size(wPLI_yeo_epoch{8,4},1)
        wPLI_yeo_sub001_memory_lowbeta_avg(i,j) = mean([wPLI_yeo_sub001_memory_lowbeta{1}(i,j) wPLI_yeo_sub001_memory_lowbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub001_NS_lowbeta_avg = mean(wPLI_yeo_epoch{11,4},3);

min_val = min(min([wPLI_yeo_sub001_count_lowbeta_avg(:); wPLI_yeo_sub001_memory_lowbeta_avg(:); wPLI_yeo_sub001_NS_lowbeta_avg(:)]));
max_val = max(max([wPLI_yeo_sub001_count_lowbeta_avg(:); wPLI_yeo_sub001_memory_lowbeta_avg(:); wPLI_yeo_sub001_NS_lowbeta_avg(:)]));

figure; 
subplot(3,5,1)
imagesc(wPLI_yeo_sub001_count_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,6)
imagesc(wPLI_yeo_sub001_memory_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,11)
imagesc(wPLI_yeo_sub001_NS_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub003_count_lowbeta_avg = mean(wPLI_yeo_epoch{1,4},3);

wPLI_yeo_sub003_memory_lowbeta_avg = mean(wPLI_yeo_epoch{2,4},3);

wPLI_yeo_sub003_NS_lowbeta = {mean(wPLI_yeo_epoch{3,4},3), mean(wPLI_yeo_epoch{4,4},3), mean(wPLI_yeo_epoch{5,4},3), mean(wPLI_yeo_epoch{6,4},3)};
wPLI_yeo_sub003_NS_lowbeta_avg = zeros(size(wPLI_yeo_epoch{3,4},1));
for i = 1:size(wPLI_yeo_epoch{3,4},1)
    for j = 1:size(wPLI_yeo_epoch{3,4},1)
        wPLI_yeo_sub003_NS_lowbeta_avg(i,j) = mean([wPLI_yeo_sub003_NS_lowbeta{1}(i,j) wPLI_yeo_sub003_NS_lowbeta{2}(i,j) wPLI_yeo_sub003_NS_lowbeta{3}(i,j) wPLI_yeo_sub003_NS_lowbeta{4}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub003_count_lowbeta_avg; wPLI_yeo_sub003_memory_lowbeta_avg; wPLI_yeo_sub003_NS_lowbeta_avg]));
max_val = max(max([wPLI_yeo_sub003_count_lowbeta_avg; wPLI_yeo_sub003_memory_lowbeta_avg; wPLI_yeo_sub003_NS_lowbeta_avg]));

subplot(3,5,2)
imagesc(wPLI_yeo_sub003_count_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,7)
imagesc(wPLI_yeo_sub003_memory_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,12)
imagesc(wPLI_yeo_sub003_NS_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub010_count_lowbeta = {mean(wPLI_yeo_epoch{12,4},3), mean(wPLI_yeo_epoch{14,4},3)};
wPLI_yeo_sub010_count_lowbeta_avg = zeros(size(wPLI_yeo_epoch{12,4},1));
for i = 1:size(wPLI_yeo_epoch{12,4},1)
    for j = 1:size(wPLI_yeo_epoch{12,4},1)
        wPLI_yeo_sub010_count_lowbeta_avg(i,j) = mean([wPLI_yeo_sub010_count_lowbeta{1}(i,j) wPLI_yeo_sub010_count_lowbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub010_memory_lowbeta = {mean(wPLI_yeo_epoch{13,4},3), mean(wPLI_yeo_epoch{15,4},3)};
wPLI_yeo_sub010_memory_lowbeta_avg = zeros(size(wPLI_yeo_epoch{13,4},1));
for i = 1:size(wPLI_yeo_epoch{13,4},1)
    for j = 1:size(wPLI_yeo_epoch{13,4},1)
        wPLI_yeo_sub010_memory_lowbeta_avg(i,j) = mean([wPLI_yeo_sub010_memory_lowbeta{1}(i,j) wPLI_yeo_sub010_memory_lowbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub010_NS_lowbeta = {mean(wPLI_yeo_epoch{17,4},3), mean(wPLI_yeo_epoch{18,4},3), mean(wPLI_yeo_epoch{19,4},3)};
wPLI_yeo_sub010_NS_lowbeta_avg = zeros(size(wPLI_yeo_epoch{17,4},1));
for i = 1:size(wPLI_yeo_epoch{17,4},1)
    for j = 1:size(wPLI_yeo_epoch{17,4},1)
        wPLI_yeo_sub010_NS_lowbeta_avg(i,j) = mean([wPLI_yeo_sub010_NS_lowbeta{1}(i,j) wPLI_yeo_sub010_NS_lowbeta{2}(i,j) wPLI_yeo_sub010_NS_lowbeta{3}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub010_count_lowbeta_avg(:); wPLI_yeo_sub010_memory_lowbeta_avg(:); wPLI_yeo_sub010_NS_lowbeta_avg(:)]));
max_val = max(max([wPLI_yeo_sub010_count_lowbeta_avg(:); wPLI_yeo_sub010_memory_lowbeta_avg(:); wPLI_yeo_sub010_NS_lowbeta_avg(:)]));

subplot(3,5,3)
imagesc(wPLI_yeo_sub010_count_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,8)
imagesc(wPLI_yeo_sub010_memory_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,13)
imagesc(wPLI_yeo_sub010_NS_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub029_count_lowbeta = {mean(wPLI_yeo_epoch{20,4},3), mean(wPLI_yeo_epoch{22,4},3)};
wPLI_yeo_sub029_count_lowbeta_avg = zeros(size(wPLI_yeo_epoch{20,4},1));
for i = 1:size(wPLI_yeo_epoch{20,4},1)
    for j = 1:size(wPLI_yeo_epoch{20,4},1)
        wPLI_yeo_sub029_count_lowbeta_avg(i,j) = mean([wPLI_yeo_sub029_count_lowbeta{1}(i,j) wPLI_yeo_sub029_count_lowbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub029_memory_lowbeta = {mean(wPLI_yeo_epoch{21,4},3), mean(wPLI_yeo_epoch{23,4},3)};
wPLI_yeo_sub029_memory_lowbeta_avg = zeros(size(wPLI_yeo_epoch{21,4},1));
for i = 1:size(wPLI_yeo_epoch{21,4},1)
    for j = 1:size(wPLI_yeo_epoch{21,4},1)
        wPLI_yeo_sub029_memory_lowbeta_avg(i,j) = mean([wPLI_yeo_sub029_memory_lowbeta{1}(i,j) wPLI_yeo_sub029_memory_lowbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub029_NS_lowbeta = {mean(wPLI_yeo_epoch{24,4},3), mean(wPLI_yeo_epoch{25,4},3)};
wPLI_yeo_sub029_NS_lowbeta_avg = zeros(size(wPLI_yeo_epoch{24,4},1));
for i = 1:size(wPLI_yeo_epoch{24,4},1)
    for j = 1:size(wPLI_yeo_epoch{24,4},1)
        wPLI_yeo_sub029_NS_lowbeta_avg(i,j) = mean([wPLI_yeo_sub029_NS_lowbeta{1}(i,j) wPLI_yeo_sub029_NS_lowbeta{2}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub029_count_lowbeta_avg(:); wPLI_yeo_sub029_memory_lowbeta_avg(:); wPLI_yeo_sub029_NS_lowbeta_avg(:)]));
max_val = max(max([wPLI_yeo_sub029_count_lowbeta_avg(:); wPLI_yeo_sub029_memory_lowbeta_avg(:); wPLI_yeo_sub029_NS_lowbeta_avg(:)]));

subplot(3,5,4)
imagesc(wPLI_yeo_sub029_count_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,9)
imagesc(wPLI_yeo_sub029_memory_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,14)
imagesc(wPLI_yeo_sub029_NS_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub034_count_lowbeta = {mean(wPLI_yeo_epoch{26,4},3), mean(wPLI_yeo_epoch{28,4},3)};
wPLI_yeo_sub034_count_lowbeta_avg = zeros(size(wPLI_yeo_epoch{7,4},1));
for i = 1:size(wPLI_yeo_epoch{26,4},1)
    for j = 1:size(wPLI_yeo_epoch{26,4},1)
        wPLI_yeo_sub034_count_lowbeta_avg(i,j) = mean([wPLI_yeo_sub034_count_lowbeta{1}(i,j) wPLI_yeo_sub034_count_lowbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub034_memory_lowbeta = {mean(wPLI_yeo_epoch{27,4},3), mean(wPLI_yeo_epoch{29,4},3)};
wPLI_yeo_sub034_memory_lowbeta_avg = zeros(size(wPLI_yeo_epoch{8,4},1));
for i = 1:size(wPLI_yeo_epoch{27,4},1)
    for j = 1:size(wPLI_yeo_epoch{27,4},1)
        wPLI_yeo_sub034_memory_lowbeta_avg(i,j) = mean([wPLI_yeo_sub034_memory_lowbeta{1}(i,j) wPLI_yeo_sub034_memory_lowbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub034_NS_lowbeta = {mean(wPLI_yeo_epoch{30,4},3), mean(wPLI_yeo_epoch{31,4},3)};
wPLI_yeo_sub034_NS_lowbeta_avg = zeros(size(wPLI_yeo_epoch{30,4},1));
for i = 1:size(wPLI_yeo_epoch{30,4},1)
    for j = 1:size(wPLI_yeo_epoch{30,4},1)
        wPLI_yeo_sub034_NS_lowbeta_avg(i,j) = mean([wPLI_yeo_sub034_NS_lowbeta{1}(i,j) wPLI_yeo_sub034_NS_lowbeta{2}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub034_count_lowbeta_avg(:); wPLI_yeo_sub034_memory_lowbeta_avg(:); wPLI_yeo_sub034_NS_lowbeta_avg(:)]));
max_val = max(max([wPLI_yeo_sub034_count_lowbeta_avg(:); wPLI_yeo_sub034_memory_lowbeta_avg(:); wPLI_yeo_sub034_NS_lowbeta_avg(:)]));

subplot(3,5,5)
imagesc(wPLI_yeo_sub034_count_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,10)
imagesc(wPLI_yeo_sub034_memory_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,15)
imagesc(wPLI_yeo_sub034_NS_lowbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

%% Plot wPLI - yeo - high beta

wPLI_yeo_sub001_count_highbeta = {mean(wPLI_yeo_epoch{7,5},3), mean(wPLI_yeo_epoch{9,5},3)};
wPLI_yeo_sub001_count_highbeta_avg = zeros(size(wPLI_yeo_epoch{7,5},1));
for i = 1:size(wPLI_yeo_epoch{7,5},1)
    for j = 1:size(wPLI_yeo_epoch{7,5},1)
        wPLI_yeo_sub001_count_highbeta_avg(i,j) = mean([wPLI_yeo_sub001_count_highbeta{1}(i,j) wPLI_yeo_sub001_count_highbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub001_memory_highbeta = {mean(wPLI_yeo_epoch{8,5},3), mean(wPLI_yeo_epoch{10,5},3)};
wPLI_yeo_sub001_memory_highbeta_avg = zeros(size(wPLI_yeo_epoch{8,5},1));
for i = 1:size(wPLI_yeo_epoch{8,5},1)
    for j = 1:size(wPLI_yeo_epoch{8,5},1)
        wPLI_yeo_sub001_memory_highbeta_avg(i,j) = mean([wPLI_yeo_sub001_memory_highbeta{1}(i,j) wPLI_yeo_sub001_memory_highbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub001_NS_highbeta_avg = mean(wPLI_yeo_epoch{11,5},3);

min_val = min(min([wPLI_yeo_sub001_count_highbeta_avg(:); wPLI_yeo_sub001_memory_highbeta_avg(:); wPLI_yeo_sub001_NS_highbeta_avg(:)]));
max_val = max(max([wPLI_yeo_sub001_count_highbeta_avg(:); wPLI_yeo_sub001_memory_highbeta_avg(:); wPLI_yeo_sub001_NS_highbeta_avg(:)]));

figure; 
subplot(3,5,1)
imagesc(wPLI_yeo_sub001_count_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,6)
imagesc(wPLI_yeo_sub001_memory_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,11)
imagesc(wPLI_yeo_sub001_NS_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub003_count_highbeta_avg = mean(wPLI_yeo_epoch{1,5},3);

wPLI_yeo_sub003_memory_highbeta_avg = mean(wPLI_yeo_epoch{2,5},3);

wPLI_yeo_sub003_NS_highbeta = {mean(wPLI_yeo_epoch{3,5},3), mean(wPLI_yeo_epoch{4,5},3), mean(wPLI_yeo_epoch{5,5},3), mean(wPLI_yeo_epoch{6,5},3)};
wPLI_yeo_sub003_NS_highbeta_avg = zeros(size(wPLI_yeo_epoch{3,5},1));
for i = 1:size(wPLI_yeo_epoch{3,5},1)
    for j = 1:size(wPLI_yeo_epoch{3,5},1)
        wPLI_yeo_sub003_NS_highbeta_avg(i,j) = mean([wPLI_yeo_sub003_NS_highbeta{1}(i,j) wPLI_yeo_sub003_NS_highbeta{2}(i,j) wPLI_yeo_sub003_NS_highbeta{3}(i,j) wPLI_yeo_sub003_NS_highbeta{4}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub003_count_highbeta_avg; wPLI_yeo_sub003_memory_highbeta_avg; wPLI_yeo_sub003_NS_highbeta_avg]));
max_val = max(max([wPLI_yeo_sub003_count_highbeta_avg; wPLI_yeo_sub003_memory_highbeta_avg; wPLI_yeo_sub003_NS_highbeta_avg]));

subplot(3,5,2)
imagesc(wPLI_yeo_sub003_count_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,7)
imagesc(wPLI_yeo_sub003_memory_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,12)
imagesc(wPLI_yeo_sub003_NS_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub010_count_highbeta = {mean(wPLI_yeo_epoch{12,5},3), mean(wPLI_yeo_epoch{14,5},3)};
wPLI_yeo_sub010_count_highbeta_avg = zeros(size(wPLI_yeo_epoch{12,5},1));
for i = 1:size(wPLI_yeo_epoch{12,5},1)
    for j = 1:size(wPLI_yeo_epoch{12,5},1)
        wPLI_yeo_sub010_count_highbeta_avg(i,j) = mean([wPLI_yeo_sub010_count_highbeta{1}(i,j) wPLI_yeo_sub010_count_highbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub010_memory_highbeta = {mean(wPLI_yeo_epoch{13,5},3), mean(wPLI_yeo_epoch{15,5},3)};
wPLI_yeo_sub010_memory_highbeta_avg = zeros(size(wPLI_yeo_epoch{13,5},1));
for i = 1:size(wPLI_yeo_epoch{13,5},1)
    for j = 1:size(wPLI_yeo_epoch{13,5},1)
        wPLI_yeo_sub010_memory_highbeta_avg(i,j) = mean([wPLI_yeo_sub010_memory_highbeta{1}(i,j) wPLI_yeo_sub010_memory_highbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub010_NS_highbeta = {mean(wPLI_yeo_epoch{17,5},3), mean(wPLI_yeo_epoch{18,5},3), mean(wPLI_yeo_epoch{19,5},3)};
wPLI_yeo_sub010_NS_highbeta_avg = zeros(size(wPLI_yeo_epoch{17,5},1));
for i = 1:size(wPLI_yeo_epoch{17,5},1)
    for j = 1:size(wPLI_yeo_epoch{17,5},1)
        wPLI_yeo_sub010_NS_highbeta_avg(i,j) = mean([wPLI_yeo_sub010_NS_highbeta{1}(i,j) wPLI_yeo_sub010_NS_highbeta{2}(i,j) wPLI_yeo_sub010_NS_highbeta{3}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub010_count_highbeta_avg(:); wPLI_yeo_sub010_memory_highbeta_avg(:); wPLI_yeo_sub010_NS_highbeta_avg(:)]));
max_val = max(max([wPLI_yeo_sub010_count_highbeta_avg(:); wPLI_yeo_sub010_memory_highbeta_avg(:); wPLI_yeo_sub010_NS_highbeta_avg(:)]));

subplot(3,5,3)
imagesc(wPLI_yeo_sub010_count_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,8)
imagesc(wPLI_yeo_sub010_memory_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,13)
imagesc(wPLI_yeo_sub010_NS_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub029_count_highbeta = {mean(wPLI_yeo_epoch{20,5},3), mean(wPLI_yeo_epoch{22,5},3)};
wPLI_yeo_sub029_count_highbeta_avg = zeros(size(wPLI_yeo_epoch{20,5},1));
for i = 1:size(wPLI_yeo_epoch{20,5},1)
    for j = 1:size(wPLI_yeo_epoch{20,5},1)
        wPLI_yeo_sub029_count_highbeta_avg(i,j) = mean([wPLI_yeo_sub029_count_highbeta{1}(i,j) wPLI_yeo_sub029_count_highbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub029_memory_highbeta = {mean(wPLI_yeo_epoch{21,5},3), mean(wPLI_yeo_epoch{23,5},3)};
wPLI_yeo_sub029_memory_highbeta_avg = zeros(size(wPLI_yeo_epoch{21,5},1));
for i = 1:size(wPLI_yeo_epoch{21,5},1)
    for j = 1:size(wPLI_yeo_epoch{21,5},1)
        wPLI_yeo_sub029_memory_highbeta_avg(i,j) = mean([wPLI_yeo_sub029_memory_highbeta{1}(i,j) wPLI_yeo_sub029_memory_highbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub029_NS_highbeta = {mean(wPLI_yeo_epoch{24,5},3), mean(wPLI_yeo_epoch{25,5},3)};
wPLI_yeo_sub029_NS_highbeta_avg = zeros(size(wPLI_yeo_epoch{24,5},1));
for i = 1:size(wPLI_yeo_epoch{24,5},1)
    for j = 1:size(wPLI_yeo_epoch{24,5},1)
        wPLI_yeo_sub029_NS_highbeta_avg(i,j) = mean([wPLI_yeo_sub029_NS_highbeta{1}(i,j) wPLI_yeo_sub029_NS_highbeta{2}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub029_count_highbeta_avg(:); wPLI_yeo_sub029_memory_highbeta_avg(:); wPLI_yeo_sub029_NS_highbeta_avg(:)]));
max_val = max(max([wPLI_yeo_sub029_count_highbeta_avg(:); wPLI_yeo_sub029_memory_highbeta_avg(:); wPLI_yeo_sub029_NS_highbeta_avg(:)]));

subplot(3,5,4)
imagesc(wPLI_yeo_sub029_count_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,9)
imagesc(wPLI_yeo_sub029_memory_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,14)
imagesc(wPLI_yeo_sub029_NS_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub034_count_highbeta = {mean(wPLI_yeo_epoch{26,5},3), mean(wPLI_yeo_epoch{28,5},3)};
wPLI_yeo_sub034_count_highbeta_avg = zeros(size(wPLI_yeo_epoch{7,5},1));
for i = 1:size(wPLI_yeo_epoch{26,5},1)
    for j = 1:size(wPLI_yeo_epoch{26,5},1)
        wPLI_yeo_sub034_count_highbeta_avg(i,j) = mean([wPLI_yeo_sub034_count_highbeta{1}(i,j) wPLI_yeo_sub034_count_highbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub034_memory_highbeta = {mean(wPLI_yeo_epoch{27,5},3), mean(wPLI_yeo_epoch{29,5},3)};
wPLI_yeo_sub034_memory_highbeta_avg = zeros(size(wPLI_yeo_epoch{8,5},1));
for i = 1:size(wPLI_yeo_epoch{27,5},1)
    for j = 1:size(wPLI_yeo_epoch{27,5},1)
        wPLI_yeo_sub034_memory_highbeta_avg(i,j) = mean([wPLI_yeo_sub034_memory_highbeta{1}(i,j) wPLI_yeo_sub034_memory_highbeta{2}(i,j)]);
    end
end

wPLI_yeo_sub034_NS_highbeta = {mean(wPLI_yeo_epoch{30,5},3), mean(wPLI_yeo_epoch{31,5},3)};
wPLI_yeo_sub034_NS_highbeta_avg = zeros(size(wPLI_yeo_epoch{30,5},1));
for i = 1:size(wPLI_yeo_epoch{30,5},1)
    for j = 1:size(wPLI_yeo_epoch{30,5},1)
        wPLI_yeo_sub034_NS_highbeta_avg(i,j) = mean([wPLI_yeo_sub034_NS_highbeta{1}(i,j) wPLI_yeo_sub034_NS_highbeta{2}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub034_count_highbeta_avg(:); wPLI_yeo_sub034_memory_highbeta_avg(:); wPLI_yeo_sub034_NS_highbeta_avg(:)]));
max_val = max(max([wPLI_yeo_sub034_count_highbeta_avg(:); wPLI_yeo_sub034_memory_highbeta_avg(:); wPLI_yeo_sub034_NS_highbeta_avg(:)]));

subplot(3,5,5)
imagesc(wPLI_yeo_sub034_count_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,10)
imagesc(wPLI_yeo_sub034_memory_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,15)
imagesc(wPLI_yeo_sub034_NS_highbeta_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

%% Plot wPLI - yeo - gamma

wPLI_yeo_sub001_count_gamma = {mean(wPLI_yeo_epoch{7,6},3), mean(wPLI_yeo_epoch{9,6},3)};
wPLI_yeo_sub001_count_gamma_avg = zeros(size(wPLI_yeo_epoch{7,6},1));
for i = 1:size(wPLI_yeo_epoch{7,6},1)
    for j = 1:size(wPLI_yeo_epoch{7,6},1)
        wPLI_yeo_sub001_count_gamma_avg(i,j) = mean([wPLI_yeo_sub001_count_gamma{1}(i,j) wPLI_yeo_sub001_count_gamma{2}(i,j)]);
    end
end

wPLI_yeo_sub001_memory_gamma = {mean(wPLI_yeo_epoch{8,6},3), mean(wPLI_yeo_epoch{10,6},3)};
wPLI_yeo_sub001_memory_gamma_avg = zeros(size(wPLI_yeo_epoch{8,6},1));
for i = 1:size(wPLI_yeo_epoch{8,6},1)
    for j = 1:size(wPLI_yeo_epoch{8,6},1)
        wPLI_yeo_sub001_memory_gamma_avg(i,j) = mean([wPLI_yeo_sub001_memory_gamma{1}(i,j) wPLI_yeo_sub001_memory_gamma{2}(i,j)]);
    end
end

wPLI_yeo_sub001_NS_gamma_avg = mean(wPLI_yeo_epoch{11,6},3);

min_val = min(min([wPLI_yeo_sub001_count_gamma_avg(:); wPLI_yeo_sub001_memory_gamma_avg(:); wPLI_yeo_sub001_NS_gamma_avg(:)]));
max_val = max(max([wPLI_yeo_sub001_count_gamma_avg(:); wPLI_yeo_sub001_memory_gamma_avg(:); wPLI_yeo_sub001_NS_gamma_avg(:)]));

figure; 
subplot(3,5,1)
imagesc(wPLI_yeo_sub001_count_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,6)
imagesc(wPLI_yeo_sub001_memory_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,11)
imagesc(wPLI_yeo_sub001_NS_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub003_count_gamma_avg = mean(wPLI_yeo_epoch{1,6},3);

wPLI_yeo_sub003_memory_gamma_avg = mean(wPLI_yeo_epoch{2,6},3);

wPLI_yeo_sub003_NS_gamma = {mean(wPLI_yeo_epoch{3,6},3), mean(wPLI_yeo_epoch{4,6},3), mean(wPLI_yeo_epoch{5,6},3), mean(wPLI_yeo_epoch{6,6},3)};
wPLI_yeo_sub003_NS_gamma_avg = zeros(size(wPLI_yeo_epoch{3,6},1));
for i = 1:size(wPLI_yeo_epoch{3,6},1)
    for j = 1:size(wPLI_yeo_epoch{3,6},1)
        wPLI_yeo_sub003_NS_gamma_avg(i,j) = mean([wPLI_yeo_sub003_NS_gamma{1}(i,j) wPLI_yeo_sub003_NS_gamma{2}(i,j) wPLI_yeo_sub003_NS_gamma{3}(i,j) wPLI_yeo_sub003_NS_gamma{4}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub003_count_gamma_avg; wPLI_yeo_sub003_memory_gamma_avg; wPLI_yeo_sub003_NS_gamma_avg]));
max_val = max(max([wPLI_yeo_sub003_count_gamma_avg; wPLI_yeo_sub003_memory_gamma_avg; wPLI_yeo_sub003_NS_gamma_avg]));

subplot(3,5,2)
imagesc(wPLI_yeo_sub003_count_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,7)
imagesc(wPLI_yeo_sub003_memory_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,12)
imagesc(wPLI_yeo_sub003_NS_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub010_count_gamma = {mean(wPLI_yeo_epoch{12,6},3), mean(wPLI_yeo_epoch{14,6},3)};
wPLI_yeo_sub010_count_gamma_avg = zeros(size(wPLI_yeo_epoch{12,6},1));
for i = 1:size(wPLI_yeo_epoch{12,6},1)
    for j = 1:size(wPLI_yeo_epoch{12,6},1)
        wPLI_yeo_sub010_count_gamma_avg(i,j) = mean([wPLI_yeo_sub010_count_gamma{1}(i,j) wPLI_yeo_sub010_count_gamma{2}(i,j)]);
    end
end

wPLI_yeo_sub010_memory_gamma = {mean(wPLI_yeo_epoch{13,6},3), mean(wPLI_yeo_epoch{15,6},3)};
wPLI_yeo_sub010_memory_gamma_avg = zeros(size(wPLI_yeo_epoch{13,6},1));
for i = 1:size(wPLI_yeo_epoch{13,6},1)
    for j = 1:size(wPLI_yeo_epoch{13,6},1)
        wPLI_yeo_sub010_memory_gamma_avg(i,j) = mean([wPLI_yeo_sub010_memory_gamma{1}(i,j) wPLI_yeo_sub010_memory_gamma{2}(i,j)]);
    end
end

wPLI_yeo_sub010_NS_gamma = {mean(wPLI_yeo_epoch{17,6},3), mean(wPLI_yeo_epoch{18,6},3), mean(wPLI_yeo_epoch{19,6},3)};
wPLI_yeo_sub010_NS_gamma_avg = zeros(size(wPLI_yeo_epoch{17,6},1));
for i = 1:size(wPLI_yeo_epoch{17,6},1)
    for j = 1:size(wPLI_yeo_epoch{17,6},1)
        wPLI_yeo_sub010_NS_gamma_avg(i,j) = mean([wPLI_yeo_sub010_NS_gamma{1}(i,j) wPLI_yeo_sub010_NS_gamma{2}(i,j) wPLI_yeo_sub010_NS_gamma{3}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub010_count_gamma_avg(:); wPLI_yeo_sub010_memory_gamma_avg(:); wPLI_yeo_sub010_NS_gamma_avg(:)]));
max_val = max(max([wPLI_yeo_sub010_count_gamma_avg(:); wPLI_yeo_sub010_memory_gamma_avg(:); wPLI_yeo_sub010_NS_gamma_avg(:)]));

subplot(3,5,3)
imagesc(wPLI_yeo_sub010_count_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,8)
imagesc(wPLI_yeo_sub010_memory_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,13)
imagesc(wPLI_yeo_sub010_NS_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub029_count_gamma = {mean(wPLI_yeo_epoch{20,6},3), mean(wPLI_yeo_epoch{22,6},3)};
wPLI_yeo_sub029_count_gamma_avg = zeros(size(wPLI_yeo_epoch{20,6},1));
for i = 1:size(wPLI_yeo_epoch{20,6},1)
    for j = 1:size(wPLI_yeo_epoch{20,6},1)
        wPLI_yeo_sub029_count_gamma_avg(i,j) = mean([wPLI_yeo_sub029_count_gamma{1}(i,j) wPLI_yeo_sub029_count_gamma{2}(i,j)]);
    end
end

wPLI_yeo_sub029_memory_gamma = {mean(wPLI_yeo_epoch{21,6},3), mean(wPLI_yeo_epoch{23,6},3)};
wPLI_yeo_sub029_memory_gamma_avg = zeros(size(wPLI_yeo_epoch{21,6},1));
for i = 1:size(wPLI_yeo_epoch{21,6},1)
    for j = 1:size(wPLI_yeo_epoch{21,6},1)
        wPLI_yeo_sub029_memory_gamma_avg(i,j) = mean([wPLI_yeo_sub029_memory_gamma{1}(i,j) wPLI_yeo_sub029_memory_gamma{2}(i,j)]);
    end
end

wPLI_yeo_sub029_NS_gamma = {mean(wPLI_yeo_epoch{24,6},3), mean(wPLI_yeo_epoch{25,6},3)};
wPLI_yeo_sub029_NS_gamma_avg = zeros(size(wPLI_yeo_epoch{24,6},1));
for i = 1:size(wPLI_yeo_epoch{24,6},1)
    for j = 1:size(wPLI_yeo_epoch{24,6},1)
        wPLI_yeo_sub029_NS_gamma_avg(i,j) = mean([wPLI_yeo_sub029_NS_gamma{1}(i,j) wPLI_yeo_sub029_NS_gamma{2}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub029_count_gamma_avg(:); wPLI_yeo_sub029_memory_gamma_avg(:); wPLI_yeo_sub029_NS_gamma_avg(:)]));
max_val = max(max([wPLI_yeo_sub029_count_gamma_avg(:); wPLI_yeo_sub029_memory_gamma_avg(:); wPLI_yeo_sub029_NS_gamma_avg(:)]));

subplot(3,5,4)
imagesc(wPLI_yeo_sub029_count_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,9)
imagesc(wPLI_yeo_sub029_memory_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,14)
imagesc(wPLI_yeo_sub029_NS_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

wPLI_yeo_sub034_count_gamma = {mean(wPLI_yeo_epoch{26,6},3), mean(wPLI_yeo_epoch{28,6},3)};
wPLI_yeo_sub034_count_gamma_avg = zeros(size(wPLI_yeo_epoch{7,6},1));
for i = 1:size(wPLI_yeo_epoch{26,6},1)
    for j = 1:size(wPLI_yeo_epoch{26,6},1)
        wPLI_yeo_sub034_count_gamma_avg(i,j) = mean([wPLI_yeo_sub034_count_gamma{1}(i,j) wPLI_yeo_sub034_count_gamma{2}(i,j)]);
    end
end

wPLI_yeo_sub034_memory_gamma = {mean(wPLI_yeo_epoch{27,6},3), mean(wPLI_yeo_epoch{29,6},3)};
wPLI_yeo_sub034_memory_gamma_avg = zeros(size(wPLI_yeo_epoch{8,6},1));
for i = 1:size(wPLI_yeo_epoch{27,6},1)
    for j = 1:size(wPLI_yeo_epoch{27,6},1)
        wPLI_yeo_sub034_memory_gamma_avg(i,j) = mean([wPLI_yeo_sub034_memory_gamma{1}(i,j) wPLI_yeo_sub034_memory_gamma{2}(i,j)]);
    end
end

wPLI_yeo_sub034_NS_gamma = {mean(wPLI_yeo_epoch{30,6},3), mean(wPLI_yeo_epoch{31,6},3)};
wPLI_yeo_sub034_NS_gamma_avg = zeros(size(wPLI_yeo_epoch{30,6},1));
for i = 1:size(wPLI_yeo_epoch{30,6},1)
    for j = 1:size(wPLI_yeo_epoch{30,6},1)
        wPLI_yeo_sub034_NS_gamma_avg(i,j) = mean([wPLI_yeo_sub034_NS_gamma{1}(i,j) wPLI_yeo_sub034_NS_gamma{2}(i,j)]);
    end
end

min_val = min(min([wPLI_yeo_sub034_count_gamma_avg(:); wPLI_yeo_sub034_memory_gamma_avg(:); wPLI_yeo_sub034_NS_gamma_avg(:)]));
max_val = max(max([wPLI_yeo_sub034_count_gamma_avg(:); wPLI_yeo_sub034_memory_gamma_avg(:); wPLI_yeo_sub034_NS_gamma_avg(:)]));

subplot(3,5,5)
imagesc(wPLI_yeo_sub034_count_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,10)
imagesc(wPLI_yeo_sub034_memory_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

subplot(3,5,15)
imagesc(wPLI_yeo_sub034_NS_gamma_avg); colorbar; clim([min_val max_val])
xticks(1:7);
yticks(1:7);
xticklabels(xlabels);
yticklabels(ylabels);
hold on;

%% Global and local efficiency

addpath(genpath("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Time/Scripts/BCT/2019_03_03_BCT"))

E_global = cell(length(d_combined),6);
E_local = cell(length(d_combined),6);

ip = java.net.InetAddress.getLocalHost.getHostAddress().string;
pctconfig('hostname',ip);
parpool;

parfor f = 1:length(d_combined)
    for fr = 1:6
        wPLI = wPLI_epoch{f,fr};
        E_global{f,fr} = zeros(size(wPLI,3),1);
        E_local{f,fr} = zeros(size(wPLI,3),num_reg);        
        for i = 1:length(E_global{f,fr})
            E_global{f,fr}(i) = efficiency_wei(squeeze(wPLI(:,:,i)));
            E_local{f,fr}(i,:) = efficiency_wei(squeeze(wPLI(:,:,i)),2);
        end
        fr
    end
    f
end


%% Modularity

modularity = cell(length(d_combined),6);
Ci = cell(length(d_combined),6);

for f = 1:length(d_combined)
    for fr = 1:6
        wPLI = wPLI_epoch{f,fr};
        modularity{f,fr} = zeros(size(wPLI,3),1);
        Ci{f,fr} = zeros(size(wPLI,3),num_reg);
        for i = 1:size(wPLI,3)
            [Ci{f,fr}(i,:), modularity{f,fr}(i)] = community_louvain(squeeze(wPLI(:,:,i)));
        end
        fr
    end
    f
end

%% Participation coefficient

participation = cell(length(d_combined),5);

for f = 1:length(d_combined)
    for fr = 1:6
        wPLI = wPLI_epoch{f,fr};
        participation{f,fr} = zeros(size(wPLI,3),num_reg);
        for i = 1:size(wPLI,3)
            participation{f,fr}(i,:) = participation_coef(squeeze(wPLI(:,:,i)), Ci{f,fr}(i,:));
        end
        fr
    end
    f
end

%% MST - diameter & eccentricity

diameter = cell(length(d_combined),5);
eccentricity = cell(length(d_combined),5);

for f = 1:length(d_combined)
    for fr = 1:6
        wPLI = wPLI_epoch{f,fr};
        diameter{f,fr} = zeros(size(wPLI,3),1);
        eccentricity{f,fr} = zeros(size(wPLI,3),100);
        for i = 1:size(wPLI,3)
            T = minspantree(graph(squeeze(wPLI(:,:,i))));
            [~,~,eccentricity{f,fr}(i,:),~,diameter{f,fr}(i)] = charpath(full(adjacency(T, 'weighted')));
        end
        fr
    end
    f
end

%% Statistical analysis - global efficiency - within-subject

n_perm = 5000;

pval_sub001_count = zeros(6,1);
pval_sub001_memory = zeros(6,1);
pval_sub003_count = zeros(6,1);
pval_sub003_memory = zeros(6,1);
pval_sub010_count = zeros(6,1);
pval_sub010_memory = zeros(6,1);
pval_sub029_count = zeros(6,1);
pval_sub029_memory = zeros(6,1);
pval_sub034_count = zeros(6,1);
pval_sub034_memory = zeros(6,1);

for fr = 1:6
    sub001_NS = E_global{11,fr};
    sub001_count = [E_global{7,fr}; E_global{9,fr}];
    pval_sub001_count(fr) = permutation_test(sub001_NS, sub001_count, n_perm);
    sub001_memory = [E_global{8,fr}; E_global{10,fr}];
    pval_sub001_memory(fr) = permutation_test(sub001_NS, sub001_memory, n_perm);

    sub003_NS = [E_global{3,fr}; E_global{4,fr}; E_global{5,fr}; E_global{6,fr}];
    sub003_count = E_global{1,fr};
    pval_sub003_count(fr) = permutation_test(sub003_NS, sub003_count, n_perm);
    sub003_memory = E_global{2,fr};
    pval_sub003_memory(fr) = permutation_test(sub003_NS, sub003_memory, n_perm);
    
    sub010_NS = [E_global{17,fr}; E_global{18,fr}; E_global{19,fr}];
    sub010_count = [E_global{12,fr}; E_global{14,fr}];
    pval_sub010_count(fr) = permutation_test(sub010_NS, sub010_count, n_perm);
    sub010_memory = [E_global{13,fr}; E_global{15,fr}];
    pval_sub010_memory(fr) = permutation_test(sub010_NS, sub010_memory, n_perm);

    sub029_NS = [E_global{24,fr}; E_global{25,fr}];
    sub029_count = [E_global{20,fr}; E_global{22,fr}];
    pval_sub029_count(fr) = permutation_test(sub029_NS, sub029_count, n_perm);
    sub029_memory = [E_global{21,fr}; E_global{23,fr}];
    pval_sub029_memory(fr) = permutation_test(sub029_NS, sub029_memory, n_perm);

    sub034_NS = [E_global{30,fr}; E_global{31,fr}];
    sub034_count = [E_global{26,fr}; E_global{28,fr}];
    pval_sub034_count(fr) = permutation_test(sub034_NS, sub034_count, n_perm);
    sub034_memory = [E_global{27,fr}; E_global{29,fr}];
    pval_sub034_memory(fr) = permutation_test(sub034_NS, sub034_memory, n_perm);
end

[~,~,~,adj_pval_sub001_E_global] = fdr_bh_2015([pval_sub001_count; pval_sub001_memory]);
[~,~,~,adj_pval_sub003_E_global] = fdr_bh_2015([pval_sub003_count; pval_sub003_memory]);
[~,~,~,adj_pval_sub010_E_global] = fdr_bh_2015([pval_sub010_count; pval_sub010_memory]);
[~,~,~,adj_pval_sub029_E_global] = fdr_bh_2015([pval_sub029_count; pval_sub029_memory]);
[~,~,~,adj_pval_sub034_E_global] = fdr_bh_2015([pval_sub034_count; pval_sub034_memory]);

%% Statistical analysis - mean local efficiency - within-subject

n_perm = 5000;

pval_sub001_count = zeros(6,1);
pval_sub001_memory = zeros(6,1);
pval_sub003_count = zeros(6,1);
pval_sub003_memory = zeros(6,1);
pval_sub010_count = zeros(6,1);
pval_sub010_memory = zeros(6,1);
pval_sub029_count = zeros(6,1);
pval_sub029_memory = zeros(6,1);
pval_sub034_count = zeros(6,1);
pval_sub034_memory = zeros(6,1);

for fr = 1:6
    sub001_NS = mean(E_local{11,fr},2);
    sub001_count = [mean(E_local{7,fr},2); mean(E_local{9,fr},2)];
    pval_sub001_count(fr) = permutation_test(sub001_NS, sub001_count, n_perm);
    sub001_memory = [mean(E_local{8,fr},2); mean(E_local{10,fr},2)];
    pval_sub001_memory(fr) = permutation_test(sub001_NS, sub001_memory, n_perm);

    sub003_NS = [mean(E_local{3,fr},2); mean(E_local{4,fr},2); mean(E_local{5,fr},2); mean(E_local{6,fr},2)];
    sub003_count = mean(E_local{1,fr},2);
    pval_sub003_count(fr) = permutation_test(sub003_NS, sub003_count, n_perm);
    sub003_memory = mean(E_local{2,fr},2);
    pval_sub003_memory(fr) = permutation_test(sub003_NS, sub003_memory, n_perm);
    
    sub010_NS = [mean(E_local{17,fr},2); mean(E_local{18,fr},2); mean(E_local{19,fr},2)];
    sub010_count = [mean(E_local{12,fr},2); mean(E_local{14,fr},2)];
    pval_sub010_count(fr) = permutation_test(sub010_NS, sub010_count, n_perm);
    sub010_memory = [mean(E_local{13,fr},2); mean(E_local{15,fr},2)];
    pval_sub010_memory(fr) = permutation_test(sub010_NS, sub010_memory, n_perm);

    sub029_NS = [mean(E_local{24,fr},2); mean(E_local{25,fr},2)];
    sub029_count = [mean(E_local{20,fr},2); mean(E_local{22,fr},2)];
    pval_sub029_count(fr) = permutation_test(sub029_NS, sub029_count, n_perm);
    sub029_memory = [mean(E_local{21,fr},2); mean(E_local{23,fr},2)];
    pval_sub029_memory(fr) = permutation_test(sub029_NS, sub029_memory, n_perm);   

    sub034_NS = [mean(E_local{30,fr},2); mean(E_local{31,fr},2)];
    sub034_count = [mean(E_local{26,fr},2); mean(E_local{28,fr},2)];
    pval_sub034_count(fr) = permutation_test(sub034_NS, sub034_count, n_perm);
    sub034_memory = [mean(E_local{27,fr},2); mean(E_local{29,fr},2)];
    pval_sub034_memory(fr) = permutation_test(sub034_NS, sub034_memory, n_perm);   
end

[~,~,~,adj_pval_sub001_E_local] = fdr_bh_2015([pval_sub001_count; pval_sub001_memory]);
[~,~,~,adj_pval_sub003_E_local] = fdr_bh_2015([pval_sub003_count; pval_sub003_memory]);
[~,~,~,adj_pval_sub010_E_local] = fdr_bh_2015([pval_sub010_count; pval_sub010_memory]);
[~,~,~,adj_pval_sub029_E_local] = fdr_bh_2015([pval_sub029_count; pval_sub029_memory]);
[~,~,~,adj_pval_sub034_E_local] = fdr_bh_2015([pval_sub034_count; pval_sub034_memory]);

%% Statistical analysis - modularity - within-subject

n_perm = 5000;
pval_sub001_count = zeros(6,1);
pval_sub001_memory = zeros(6,1);
pval_sub003_count = zeros(6,1);
pval_sub003_memory = zeros(6,1);
pval_sub010_count = zeros(6,1);
pval_sub010_memory = zeros(6,1);
pval_sub029_count = zeros(6,1);
pval_sub029_memory = zeros(6,1);
pval_sub034_count = zeros(6,1);
pval_sub034_memory = zeros(6,1);

for fr = 1:6
    sub001_NS = modularity{11,fr};
    sub001_count = [modularity{7,fr}; modularity{9,fr}];
    pval_sub001_count(fr) = permutation_test(sub001_NS, sub001_count, n_perm);
    sub001_memory = [modularity{8,fr}; modularity{10,fr}];
    pval_sub001_memory(fr) = permutation_test(sub001_NS, sub001_memory, n_perm);

    sub003_NS = [modularity{3,fr}; modularity{4,fr}; modularity{5,fr}; modularity{6,fr}];
    sub003_count = modularity{1,fr};
    pval_sub003_count(fr) = permutation_test(sub003_NS, sub003_count, n_perm);
    sub003_memory = modularity{2,fr};
    pval_sub003_memory(fr) = permutation_test(sub003_NS, sub003_memory, n_perm);
    
    sub010_NS = [modularity{17,fr}; modularity{18,fr}; modularity{19,fr}];
    sub010_count = [modularity{12,fr}; modularity{14,fr}];
    pval_sub010_count(fr) = permutation_test(sub010_NS, sub010_count, n_perm);
    sub010_memory = [modularity{13,fr}; modularity{15,fr}];
    pval_sub010_memory(fr) = permutation_test(sub010_NS, sub010_memory, n_perm);

    sub029_NS = [modularity{24,fr}; modularity{25,fr}];
    sub029_count = [modularity{20,fr}; modularity{22,fr}];
    pval_sub029_count(fr) = permutation_test(sub029_NS, sub029_count, n_perm);
    sub029_memory = [modularity{21,fr}; modularity{23,fr}];
    pval_sub029_memory(fr) = permutation_test(sub029_NS, sub029_memory, n_perm);

    sub034_NS = [modularity{30,fr}; modularity{31,fr}];
    sub034_count = [modularity{26,fr}; modularity{28,fr}];
    pval_sub034_count(fr) = permutation_test(sub034_NS, sub034_count, n_perm);
    sub034_memory = [modularity{27,fr}; modularity{29,fr}];
    pval_sub034_memory(fr) = permutation_test(sub034_NS, sub034_memory, n_perm);
end

[~,~,~,adj_pval_sub001_modularity] = fdr_bh_2015([pval_sub001_count; pval_sub001_memory]);
[~,~,~,adj_pval_sub003_modularity] = fdr_bh_2015([pval_sub003_count; pval_sub003_memory]);
[~,~,~,adj_pval_sub010_modularity] = fdr_bh_2015([pval_sub010_count; pval_sub010_memory]);
[~,~,~,adj_pval_sub029_modularity] = fdr_bh_2015([pval_sub029_count; pval_sub029_memory]);
[~,~,~,adj_pval_sub034_modularity] = fdr_bh_2015([pval_sub034_count; pval_sub034_memory]);

%% Statistical analysis - standard deviation of participation coefficient - within-subject

pval_sub001_count = zeros(6,1);
pval_sub001_memory = zeros(6,1);
pval_sub003_count = zeros(6,1);
pval_sub003_memory = zeros(6,1);
pval_sub010_count = zeros(6,1);
pval_sub010_memory = zeros(6,1);
pval_sub029_count = zeros(6,1);
pval_sub029_memory = zeros(6,1);
pval_sub034_count = zeros(6,1);
pval_sub034_memory = zeros(6,1);

for fr = 1:6
    sub001_NS = std(participation{11,fr},0,2);
    sub001_count = [std(participation{7,fr},0,2); std(participation{9,fr},0,2)];
    pval_sub001_count(fr) = permutation_test(sub001_NS, sub001_count, n_perm);
    sub001_memory = [std(participation{8,fr},0,2); std(participation{10,fr},0,2)];
    pval_sub001_memory(fr) = permutation_test(sub001_NS, sub001_memory, n_perm);
    
    sub003_NS = [std(participation{3,fr},0,2); std(participation{4,fr},0,2); std(participation{5,fr},0,2); std(participation{6,fr},0,2)];
    sub003_count = std(participation{1,fr},0,2);
    pval_sub003_count(fr) = permutation_test(sub003_NS, sub003_count, n_perm);
    sub003_memory = std(participation{2,fr},0,2);
    pval_sub003_memory(fr) = permutation_test(sub003_NS, sub003_memory, n_perm);
    
    sub010_NS = [std(participation{17,fr},0,2); std(participation{18,fr},0,2); std(participation{19,fr},0,2)];
    sub010_count = [std(participation{12,fr},0,2); std(participation{14,fr},0,2)];
    pval_sub010_count(fr) = permutation_test(sub010_NS, sub010_count, n_perm);
    sub010_memory = [std(participation{13,fr},0,2); std(participation{15,fr},0,2)];
    pval_sub010_memory(fr) = permutation_test(sub010_NS, sub010_memory, n_perm);

    sub029_NS = [std(participation{24,fr},0,2); std(participation{25,fr},0,2)];
    sub029_count = [std(participation{20,fr},0,2); std(participation{22,fr},0,2)];
    pval_sub029_count(fr) = permutation_test(sub029_NS, sub029_count, n_perm);
    sub029_memory = [std(participation{21,fr},0,2); std(participation{23,fr},0,2)];
    pval_sub029_memory(fr) = permutation_test(sub029_NS, sub029_memory, n_perm);    

    sub034_NS = [std(participation{30,fr},0,2); std(participation{31,fr},0,2)];
    sub034_count = [std(participation{26,fr},0,2); std(participation{28,fr},0,2)];
    pval_sub034_count(fr) = permutation_test(sub034_NS, sub034_count, n_perm);
    sub034_memory = [std(participation{27,fr},0,2); std(participation{29,fr},0,2)];
    pval_sub034_memory(fr) = permutation_test(sub034_NS, sub034_memory, n_perm);   
end

[~,~,~,adj_pval_sub001_participation] = fdr_bh_2015([pval_sub001_count; pval_sub001_memory]);
[~,~,~,adj_pval_sub003_participation] = fdr_bh_2015([pval_sub003_count; pval_sub003_memory]);
[~,~,~,adj_pval_sub010_participation] = fdr_bh_2015([pval_sub010_count; pval_sub010_memory]);
[~,~,~,adj_pval_sub029_participation] = fdr_bh_2015([pval_sub029_count; pval_sub029_memory]);
[~,~,~,adj_pval_sub034_participation] = fdr_bh_2015([pval_sub034_count; pval_sub034_memory]);


%% Statistical analysis - diameter - within-subject

n_perm = 5000;
pval_sub001_count = zeros(6,1);
pval_sub001_memory = zeros(6,1);
pval_sub003_count = zeros(6,1);
pval_sub003_memory = zeros(6,1);
pval_sub010_count = zeros(6,1);
pval_sub010_memory = zeros(6,1);
pval_sub029_count = zeros(6,1);
pval_sub029_memory = zeros(6,1);
pval_sub034_count = zeros(6,1);
pval_sub034_memory = zeros(6,1);

for fr = 1:6
    sub001_NS = diameter{11,fr};
    sub001_count = [diameter{7,fr}; diameter{9,fr}];
    pval_sub001_count(fr) = permutation_test(sub001_NS, sub001_count, n_perm);
    sub001_memory = [diameter{8,fr}; diameter{10,fr}];
    pval_sub001_memory(fr) = permutation_test(sub001_NS, sub001_memory, n_perm);

    sub003_NS = [diameter{3,fr}; diameter{4,fr}; diameter{5,fr}; diameter{6,fr}];
    sub003_count = diameter{1,fr};
    pval_sub003_count(fr) = permutation_test(sub003_NS, sub003_count, n_perm);
    sub003_memory = diameter{2,fr};
    pval_sub003_memory(fr) = permutation_test(sub003_NS, sub003_memory, n_perm);    
    
    sub010_NS = [diameter{17,fr}; diameter{18,fr}; diameter{19,fr}];
    sub010_count = [diameter{12,fr}; diameter{14,fr}];
    pval_sub010_count(fr) = permutation_test(sub010_NS, sub010_count, n_perm);
    sub010_memory = [diameter{13,fr}; diameter{15,fr}];
    pval_sub010_memory(fr) = permutation_test(sub010_NS, sub010_memory, n_perm);

    sub029_NS = [diameter{24,fr}; diameter{25,fr}];
    sub029_count = [diameter{20,fr}; diameter{22,fr}];
    pval_sub029_count(fr) = permutation_test(sub029_NS, sub029_count, n_perm);
    sub029_memory = [diameter{21,fr}; diameter{23,fr}];
    pval_sub029_memory(fr) = permutation_test(sub029_NS, sub029_memory, n_perm);

    sub034_NS = [diameter{24,fr}; diameter{25,fr}];
    sub034_count = [diameter{20,fr}; diameter{22,fr}];
    pval_sub034_count(fr) = permutation_test(sub034_NS, sub034_count, n_perm);
    sub034_memory = [diameter{21,fr}; diameter{23,fr}];
    pval_sub034_memory(fr) = permutation_test(sub034_NS, sub034_memory, n_perm);    
end

[~,~,~,adj_pval_sub001_diameter] = fdr_bh_2015([pval_sub001_count; pval_sub001_memory]);
[~,~,~,adj_pval_sub003_diameter] = fdr_bh_2015([pval_sub003_count; pval_sub003_memory]);
[~,~,~,adj_pval_sub010_diameter] = fdr_bh_2015([pval_sub010_count; pval_sub010_memory]);
[~,~,~,adj_pval_sub029_diameter] = fdr_bh_2015([pval_sub029_count; pval_sub029_memory]);
[~,~,~,adj_pval_sub034_diameter] = fdr_bh_2015([pval_sub034_count; pval_sub034_memory]);

%% Statistical analysis - mean eccentricity - within-subject

n_perm = 5000;

pval_sub001_count = zeros(6,1);
pval_sub001_memory = zeros(6,1);
pval_sub003_count = zeros(6,1);
pval_sub003_memory = zeros(6,1);
pval_sub010_count = zeros(6,1);
pval_sub010_memory = zeros(6,1);
pval_sub029_count = zeros(6,1);
pval_sub029_memory = zeros(6,1);
pval_sub034_count = zeros(6,1);
pval_sub034_memory = zeros(6,1);

for fr = 1:6    
    sub001_NS = mean(eccentricity{11,fr},2);
    sub001_count = [mean(eccentricity{7,fr},2); mean(eccentricity{9,fr},2)];
    pval_sub001_count(fr) = permutation_test(sub001_NS, sub001_count, n_perm);
    sub001_memory = [mean(eccentricity{8,fr},2); mean(eccentricity{10,fr},2)];
    pval_sub001_memory(fr) = permutation_test(sub001_NS, sub001_memory, n_perm);

    sub003_NS = [mean(eccentricity{3,fr},2); mean(eccentricity{4,fr},2); mean(eccentricity{5,fr},2); mean(eccentricity{6,fr},2)];
    sub003_count = mean(eccentricity{1,fr},2);
    pval_sub003_count(fr) = permutation_test(sub003_NS, sub003_count, n_perm);
    sub003_memory = mean(eccentricity{2,fr},2);
    pval_sub003_memory(fr) = permutation_test(sub003_NS, sub003_memory, n_perm);
    
    sub010_NS = [mean(eccentricity{17,fr},2); mean(eccentricity{18,fr},2); mean(eccentricity{19,fr},2)];
    sub010_count = [mean(eccentricity{12,fr},2); mean(eccentricity{14,fr},2)];
    pval_sub010_count(fr) = permutation_test(sub010_NS, sub010_count, n_perm);
    sub010_memory = [mean(eccentricity{13,fr},2); mean(eccentricity{15,fr},2)];
    pval_sub010_memory(fr) = permutation_test(sub010_NS, sub010_memory, n_perm);

    sub029_NS = [mean(eccentricity{24,fr},2); mean(eccentricity{25,fr},2)];
    sub029_count = [mean(eccentricity{20,fr},2); mean(eccentricity{22,fr},2)];
    pval_sub029_count(fr) = permutation_test(sub029_NS, sub029_count, n_perm);
    sub029_memory = [mean(eccentricity{21,fr},2); mean(eccentricity{23,fr},2)];
    pval_sub029_memory(fr) = permutation_test(sub029_NS, sub029_memory, n_perm);  

    sub034_NS = [mean(eccentricity{30,fr},2); mean(eccentricity{31,fr},2)];
    sub034_count = [mean(eccentricity{26,fr},2); mean(eccentricity{28,fr},2)];
    pval_sub034_count(fr) = permutation_test(sub034_NS, sub034_count, n_perm);
    sub034_memory = [mean(eccentricity{27,fr},2); mean(eccentricity{29,fr},2)];
    pval_sub034_memory(fr) = permutation_test(sub034_NS, sub034_memory, n_perm);      
end

[~,~,~,adj_pval_sub001_eccentricity] = fdr_bh_2015([pval_sub001_count; pval_sub001_memory]);
[~,~,~,adj_pval_sub003_eccentricity] = fdr_bh_2015([pval_sub003_count; pval_sub003_memory]);
[~,~,~,adj_pval_sub010_eccentricity] = fdr_bh_2015([pval_sub010_count; pval_sub010_memory]);
[~,~,~,adj_pval_sub029_eccentricity] = fdr_bh_2015([pval_sub029_count; pval_sub029_memory]);
[~,~,~,adj_pval_sub034_eccentricity] = fdr_bh_2015([pval_sub034_count; pval_sub034_memory]);

%% Network properties plotting parameters

subject_labels   = {'sub001', 'sub002', 'sub003', 'sub004', 'sub005'};
condition_labels = {'Counting', 'Memory', 'NS'};

% Column 1: subject index; column 2: condition index
map = [2 1; 2 2; 2 3; 2 3; 2 3; 2 3;
       1 1; 1 2; 1 1; 1 2; 1 3;
       3 1; 3 2; 3 1; 3 2; 3 3; 3 3; 3 3; 3 3;
       4 1; 4 2; 4 1; 4 2; 4 3; 4 3;
       5 1; 5 2; 5 1; 5 2; 5 3; 5 3]; 

%% Plot network properties

subjects   = {};
conditions = {};
all_data   = [];

for i = 1:size(E_global,1)
    d = E_global{i,3};
    all_data = [all_data; d];
    s = subject_labels{map(i,1)};
    c = condition_labels{map(i,2)};
    subjects   = [subjects; repmat({s}, numel(d), 1)];
    conditions = [conditions; repmat({c}, numel(d), 1)];
end

% E_global
figure;
subplot(2,3,1)
hold on

% Explicit mapping from subject label to numeric x-axis position
subjectNames = subject_labels;
conditionNames = condition_labels;
offsets = [-0.2, 0, 0.2];  % Counting, Memory, NS
condOffset = containers.Map(conditionNames, offsets);

% Convert subject and condition to numeric x-positions
x_positions = zeros(size(all_data));
for i = 1:numel(all_data)
    subj_idx = find(strcmp(subjects{i}, subjectNames));
    cond_idx = find(strcmp(conditions{i}, conditionNames));
    x_positions(i) = subj_idx + offsets(cond_idx);
end

% Map conditions to colors (consistent with boxchart defaults)
condColors = containers.Map( ...
    {'Counting', 'Memory', 'NS'}, ...
    {'r', [0.85 0.325 0.098], 'b'});

% Plot each condition separately to preserve color grouping
for i = 1:numel(conditionNames)
    cond = conditionNames{i};
    idx = strcmp(conditions, cond);
    boxchart(x_positions(idx), all_data(idx), ...
        'BoxFaceColor', condColors(cond), 'BoxWidth', 0.15);
end

% Axis formatting
xticks(1:numel(subjectNames))
xticklabels(subjectNames)
xlabel('Subject')
ylabel('Global Efficiency')
title('Global Efficiency (Alpha)')
legend(conditionNames, 'Location', 'northeast')

hold on
offsets = [-0.2, 0, 0.2];  % Left = Counting, Center = Memory, Right = NS
% Significance bars for selected comparisons
y_offset = 0.01;  % vertical space between bars
fontSize = 14;    % size of asterisk labels

% Mapping from condition to offset
condOffset = containers.Map(condition_labels, offsets);

% Significance comparisons: {subject_x, cond1, cond2, stars}
comparisons = {
    1, 'Counting', 'NS', '***';   % sub003: NS vs counting
    1, 'Memory',   'NS', '***';    % sub003: NS vs memory
    4, 'Counting',   'NS', '*';  % sub010: NS vs memory
};

for i = 1:size(comparisons, 1)
    subjX = comparisons{i,1};
    cond1 = comparisons{i,2};
    cond2 = comparisons{i,3};
    stars = comparisons{i,4};

    % x-locations (add offset for condition)
    x1 = subjX + condOffset(cond1);
    x2 = subjX + condOffset(cond2);

    % y-position (choose max y at this subject + some buffer)
    y1 = max(all_data(strcmp(subjects, subject_labels{subjX}))) + (i+1)*y_offset;

    % Plot the line
    h = plot([x1 x2], [y1 y1], 'k', 'LineWidth', 1.2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';  % hide from legend
    % Plot the asterisks
    text(mean([x1 x2]), y1 + y_offset - 0.008, stars, 'HorizontalAlignment', 'center', 'FontSize', fontSize)
end



subjects   = {};
conditions = {};
all_data   = [];

for i = 1:size(E_local,1)
    d = mean(E_local{i,3},2);
    all_data = [all_data; d];
    s = subject_labels{map(i,1)};
    c = condition_labels{map(i,2)};
    subjects   = [subjects; repmat({s}, numel(d), 1)];
    conditions = [conditions; repmat({c}, numel(d), 1)];
end

% Local efficiency

subplot(2,3,2)
hold on

% Explicit mapping from subject label to numeric x-axis position
subjectNames = subject_labels;
conditionNames = condition_labels;
offsets = [-0.2, 0, 0.2];  % Counting, Memory, NS
condOffset = containers.Map(conditionNames, offsets);

% Convert subject and condition to numeric x-positions
x_positions = zeros(size(all_data));
for i = 1:numel(all_data)
    subj_idx = find(strcmp(subjects{i}, subjectNames));
    cond_idx = find(strcmp(conditions{i}, conditionNames));
    x_positions(i) = subj_idx + offsets(cond_idx);
end

% Map conditions to colors (consistent with boxchart defaults)
condColors = containers.Map( ...
    {'Counting', 'Memory', 'NS'}, ...
    {'r', [0.85 0.325 0.098], 'b'});

% Plot each condition separately to preserve color grouping
for i = 1:numel(conditionNames)
    cond = conditionNames{i};
    idx = strcmp(conditions, cond);
    boxchart(x_positions(idx), all_data(idx), ...
        'BoxFaceColor', condColors(cond), 'BoxWidth', 0.15);
end

% Axis formatting
xticks(1:numel(subjectNames))
xticklabels(subjectNames)
xlabel('Subject')
ylabel('Mean Local Efficiency')
title('Mean Local Efficiency (Alpha)')
legend(conditionNames, 'Location', 'northeast')

hold on
offsets = [-0.2, 0, 0.2];  % Left = Counting, Center = Memory, Right = NS
% Significance bars for selected comparisons
y_offset = 0.01;  % vertical space between bars
fontSize = 14;    % size of asterisk labels

% Mapping from condition to offset
condOffset = containers.Map(condition_labels, offsets);

% Significance comparisons: {subject_x, cond1, cond2, stars}
comparisons = {
    1, 'Counting', 'NS', '*';
    1, 'Memory', 'NS', '*';
    3, 'Counting', 'NS', '*';  % sub010: NS vs counting
};

for i = 1:size(comparisons, 1)
    subjX = comparisons{i,1};
    cond1 = comparisons{i,2};
    cond2 = comparisons{i,3};
    stars = comparisons{i,4};

    % x-locations (add offset for condition)
    x1 = subjX + condOffset(cond1);
    x2 = subjX + condOffset(cond2);

    % y-position (choose max y at this subject + some buffer)
    y1 = max(all_data(strcmp(subjects, subject_labels{subjX}))) + (i+1)*y_offset;

    % Plot the line
    h = plot([x1 x2], [y1 y1], 'k', 'LineWidth', 1.2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';  % hide from legend
    % Plot the asterisks
    text(mean([x1 x2]), y1 + y_offset - 0.008, stars, 'HorizontalAlignment', 'center', 'FontSize', fontSize)
end


subjects   = {};
conditions = {};
all_data   = [];

for i = 1:size(modularity,1)
    d = modularity{i,3};
    all_data = [all_data; d];
    s = subject_labels{map(i,1)};
    c = condition_labels{map(i,2)};
    subjects   = [subjects; repmat({s}, numel(d), 1)];
    conditions = [conditions; repmat({c}, numel(d), 1)];
end

% Modularity
subplot(2,3,3)
hold on

% Explicit mapping from subject label to numeric x-axis position
subjectNames = subject_labels;
conditionNames = condition_labels;
offsets = [-0.2, 0, 0.2];  % Counting, Memory, NS
condOffset = containers.Map(conditionNames, offsets);

% Convert subject and condition to numeric x-positions
x_positions = zeros(size(all_data));
for i = 1:numel(all_data)
    subj_idx = find(strcmp(subjects{i}, subjectNames));
    cond_idx = find(strcmp(conditions{i}, conditionNames));
    x_positions(i) = subj_idx + offsets(cond_idx);
end

% Map conditions to colors (consistent with boxchart defaults)
condColors = containers.Map( ...
    {'Counting', 'Memory', 'NS'}, ...
    {'r', [0.85 0.325 0.098], 'b'});

% Plot each condition separately to preserve color grouping
for i = 1:numel(conditionNames)
    cond = conditionNames{i};
    idx = strcmp(conditions, cond);
    boxchart(x_positions(idx), all_data(idx), ...
        'BoxFaceColor', condColors(cond), 'BoxWidth', 0.15);
end

% Axis formatting
xticks(1:numel(subjectNames))
xticklabels(subjectNames)
xlabel('Subject')
ylabel('Modularity')
title('Modularity (Alpha)')
legend(conditionNames, 'Location', 'northeast')

hold on
offsets = [-0.2, 0, 0.2];  % Left = Counting, Center = Memory, Right = NS
% Significance bars for selected comparisons
y_offset = 0.01;  % vertical space between bars
fontSize = 14;    % size of asterisk labels

% Mapping from condition to offset
condOffset = containers.Map(condition_labels, offsets);

% Significance comparisons: {subject_x, cond1, cond2, stars}
comparisons = {
    1, 'Counting', 'NS', '***';
    1, 'Memory', 'NS', '***';
    5, 'Memory',   'NS', '*';  
};

for i = 1:size(comparisons, 1)
    subjX = comparisons{i,1};
    cond1 = comparisons{i,2};
    cond2 = comparisons{i,3};
    stars = comparisons{i,4};

    % x-locations (add offset for condition)
    x1 = subjX + condOffset(cond1);
    x2 = subjX + condOffset(cond2);

    % y-position (choose max y at this subject + some buffer)
    y1 = max(all_data(strcmp(subjects, subject_labels{subjX}))) + (i+1)*y_offset;

    % Plot the line
    h = plot([x1 x2], [y1 y1], 'k', 'LineWidth', 1.2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';  % hide from legend
    % Plot the asterisks
    text(mean([x1 x2]), y1 + y_offset - 0.008, stars, 'HorizontalAlignment', 'center', 'FontSize', fontSize)
end

subjects   = {};
conditions = {};
all_data   = [];

for i = 1:size(participation,1)
    d = std(participation{i,3},0,2);
    all_data = [all_data; d];
    s = subject_labels{map(i,1)};
    c = condition_labels{map(i,2)};
    subjects   = [subjects; repmat({s}, numel(d), 1)];
    conditions = [conditions; repmat({c}, numel(d), 1)];
end

% Standard deviation of participation coefficient
subplot(2,3,4)
hold on

% Explicit mapping from subject label to numeric x-axis position
subjectNames = subject_labels;
conditionNames = condition_labels;
offsets = [-0.2, 0, 0.2];  % Counting, Memory, NS
condOffset = containers.Map(conditionNames, offsets);

% Convert subject and condition to numeric x-positions
x_positions = zeros(size(all_data));
for i = 1:numel(all_data)
    subj_idx = find(strcmp(subjects{i}, subjectNames));
    cond_idx = find(strcmp(conditions{i}, conditionNames));
    x_positions(i) = subj_idx + offsets(cond_idx);
end

% Map conditions to colors (consistent with boxchart defaults)
condColors = containers.Map( ...
    {'Counting', 'Memory', 'NS'}, ...
    {'r', [0.85 0.325 0.098], 'b'});

% Plot each condition separately to preserve color grouping
for i = 1:numel(conditionNames)
    cond = conditionNames{i};
    idx = strcmp(conditions, cond);
    boxchart(x_positions(idx), all_data(idx), ...
        'BoxFaceColor', condColors(cond), 'BoxWidth', 0.15);
end

% Axis formatting
xticks(1:numel(subjectNames))
xticklabels(subjectNames)
xlabel('Subject')
ylabel('Standard Deviation of Participation Coefficient')
title('Standard Deviation of Participation Coefficient (Alpha)')
legend(conditionNames, 'Location', 'northeast')

hold on
offsets = [-0.2, 0, 0.2];  % Left = Counting, Center = Memory, Right = NS
% Significance bars for selected comparisons
y_offset = 0.005;  % vertical space between bars
fontSize = 14;    % size of asterisk labels

% Mapping from condition to offset
condOffset = containers.Map(condition_labels, offsets);

% % Significance comparisons: {subject_x, cond1, cond2, stars}
% comparisons = {
%     2, 'Counting', 'NS', '**';   % sub003: NS vs counting
%     2, 'Memory',   'NS', '***';    % sub003: NS vs memory
%     3, 'Memory',   'NS', '***';  % sub010: NS vs memory
% };
% 
% for i = 1:size(comparisons, 1)
%     subjX = comparisons{i,1};
%     cond1 = comparisons{i,2};
%     cond2 = comparisons{i,3};
%     stars = comparisons{i,4};
% 
%     % x-locations (add offset for condition)
%     x1 = subjX + condOffset(cond1);
%     x2 = subjX + condOffset(cond2);
% 
%     % y-position (choose max y at this subject + some buffer)
%     y1 = max(all_data(strcmp(subjects, subject_labels{subjX}))) + (i+1)*y_offset;
% 
%     % Plot the line
%     h = plot([x1 x2], [y1 y1], 'k', 'LineWidth', 1.2);
%     h.Annotation.LegendInformation.IconDisplayStyle = 'off';  % hide from legend
%     % Plot the asterisks
%     text(mean([x1 x2]), y1 + y_offset - 0.004, stars, 'HorizontalAlignment', 'center', 'FontSize', fontSize)
% end

subjects   = {};
conditions = {};
all_data   = [];

for i = 1:size(eccentricity,1)
    d = mean(eccentricity{i,3},2);
    all_data = [all_data; d];
    s = subject_labels{map(i,1)};
    c = condition_labels{map(i,2)};
    subjects   = [subjects; repmat({s}, numel(d), 1)];
    conditions = [conditions; repmat({c}, numel(d), 1)];
end

% Eccentricity of MST 
subplot(2,3,5)
hold on

% Explicit mapping from subject label to numeric x-axis position
subjectNames = subject_labels;
conditionNames = condition_labels;
offsets = [-0.2, 0, 0.2];  % Counting, Memory, NS
condOffset = containers.Map(conditionNames, offsets);

% Convert subject and condition to numeric x-positions
x_positions = zeros(size(all_data));
for i = 1:numel(all_data)
    subj_idx = find(strcmp(subjects{i}, subjectNames));
    cond_idx = find(strcmp(conditions{i}, conditionNames));
    x_positions(i) = subj_idx + offsets(cond_idx);
end

% Map conditions to colors (consistent with boxchart defaults)
condColors = containers.Map( ...
    {'Counting', 'Memory', 'NS'}, ...
    {'r', [0.85 0.325 0.098], 'b'});

% Plot each condition separately to preserve color grouping
for i = 1:numel(conditionNames)
    cond = conditionNames{i};
    idx = strcmp(conditions, cond);
    boxchart(x_positions(idx), all_data(idx), ...
        'BoxFaceColor', condColors(cond), 'BoxWidth', 0.15);
end

% Axis formatting
xticks(1:numel(subjectNames))
xticklabels(subjectNames)
xlabel('Subject')
ylabel('Mean Eccentricity of MST')
title('Mean Eccentricity of MST (Alpha)')
legend(conditionNames, 'Location', 'northeast')

hold on
offsets = [-0.2, 0, 0.2];  % Left = Counting, Center = Memory, Right = NS
% Significance bars for selected comparisons
y_offset = 0.0005;  % vertical space between bars
fontSize = 14;    % size of asterisk labels

% Mapping from condition to offset
condOffset = containers.Map(condition_labels, offsets);

% % Significance comparisons: {subject_x, cond1, cond2, stars}
% comparisons = {
%     2, 'Counting', 'NS', '*';   % sub003: NS vs counting
%     2, 'Memory', 'NS', '*';   % sub003: NS vs counting
%     3, 'Memory', 'NS', '***';   % sub010: NS vs memory    
% };
% 
% for i = 1:size(comparisons, 1)
%     subjX = comparisons{i,1};
%     cond1 = comparisons{i,2};
%     cond2 = comparisons{i,3};
%     stars = comparisons{i,4};
% 
%     % x-locations (add offset for condition)
%     x1 = subjX + condOffset(cond1);
%     x2 = subjX + condOffset(cond2);
% 
%     % y-position (choose max y at this subject + some buffer)
%     y1 = max(all_data(strcmp(subjects, subject_labels{subjX}))) + (i+1)*y_offset;
% 
%     % Plot the line
%     h = plot([x1 x2], [y1 y1], 'k', 'LineWidth', 1.2);
%     h.Annotation.LegendInformation.IconDisplayStyle = 'off';  % hide from legend
%     % Plot the asterisks
%     text(mean([x1 x2]), y1 + y_offset - 0.0004, stars, 'HorizontalAlignment', 'center', 'FontSize', fontSize)
% end

subjects   = {};
conditions = {};
all_data   = [];

for i = 1:size(diameter,1)
    d = diameter{i,3};
    all_data = [all_data; d];
    s = subject_labels{map(i,1)};
    c = condition_labels{map(i,2)};
    subjects   = [subjects; repmat({s}, numel(d), 1)];
    conditions = [conditions; repmat({c}, numel(d), 1)];
end

% Diameter of MST
subplot(2,3,6)
hold on

% Explicit mapping from subject label to numeric x-axis position
subjectNames = subject_labels;
conditionNames = condition_labels;
offsets = [-0.2, 0, 0.2];  % Counting, Memory, NS
condOffset = containers.Map(conditionNames, offsets);

% Convert subject and condition to numeric x-positions
x_positions = zeros(size(all_data));
for i = 1:numel(all_data)
    subj_idx = find(strcmp(subjects{i}, subjectNames));
    cond_idx = find(strcmp(conditions{i}, conditionNames));
    x_positions(i) = subj_idx + offsets(cond_idx);
end

% Map conditions to colors (consistent with boxchart defaults)
condColors = containers.Map( ...
    {'Counting', 'Memory', 'NS'}, ...
    {'r', [0.85 0.325 0.098], 'b'});

% Plot each condition separately to preserve color grouping
for i = 1:numel(conditionNames)
    cond = conditionNames{i};
    idx = strcmp(conditions, cond);
    boxchart(x_positions(idx), all_data(idx), ...
        'BoxFaceColor', condColors(cond), 'BoxWidth', 0.15);
end

% Axis formatting
xticks(1:numel(subjectNames))
xticklabels(subjectNames)
xlabel('Subject')
ylabel('Diameter of MST')
title('Diameter of MST (Alpha)')
legend(conditionNames, 'Location', 'northeast')

hold on
offsets = [-0.2, 0, 0.2];  % Left = Counting, Center = Memory, Right = NS
% Significance bars for selected comparisons
y_offset = 0.005;  % vertical space between bars
fontSize = 14;    % size of asterisk labels

% Mapping from condition to offset
condOffset = containers.Map(condition_labels, offsets);

% % Significance comparisons: {subject_x, cond1, cond2, stars}
% comparisons = {
%     1, 'Counting', 'NS', '**';   % sub001: NS vs counting
% };
% 
% for i = 1:size(comparisons, 1)
%     subjX = comparisons{i,1};
%     cond1 = comparisons{i,2};
%     cond2 = comparisons{i,3};
%     stars = comparisons{i,4};
% 
%     % x-locations (add offset for condition)
%     x1 = subjX + condOffset(cond1);
%     x2 = subjX + condOffset(cond2);
% 
%     % y-position (choose max y at this subject + some buffer)
%     y1 = max(all_data(strcmp(subjects, subject_labels{subjX}))) + (i+1)*y_offset;
% 
%     % Plot the line
%     h = plot([x1 x2], [y1 y1], 'k', 'LineWidth', 1.2);
%     h.Annotation.LegendInformation.IconDisplayStyle = 'off';  % hide from legend
%     % Plot the asterisks
%     text(mean([x1 x2]), y1 + y_offset - 0.004, stars, 'HorizontalAlignment', 'center', 'FontSize', fontSize)
% end


