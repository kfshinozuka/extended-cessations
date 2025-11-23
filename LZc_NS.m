addpath(genpath('/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Scripts/EntRate-master'))
addpath('/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Clayton/Scripts')
num_reg = 100;
freq_low = [1 4 8 13 30];
freq_high = [4 8 13 30 48];
fs = 200;

%% Combined - Broadband - Epoched

cd("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Results/SourceRecon/Combined/Truncated/CommonFilter")
data_dir = dir('*.mat');
LZc_avg_all = cell(length(data_dir),1);
LZc_reg_all = cell(length(data_dir),1);

ip  = java.net.InetAddress.getLocalHost.getHostAddress().string;
pctconfig('hostname',ip);
parpool;

parfor f = 1:length(data_dir)

    % Load timeseries (vs)
    vs = importdata(data_dir(f).name);
    
    % reshape vs into a 3D matrix
    [p1, Tm]            = size(vs);
    if f < 26 % not sub034 (who has very short data)
        segmentLength   = 3000;   
    else % sub034
        segmentLength   = 300;
    end    
    numFullSegments     = floor(Tm / segmentLength);
    validTimepoints     = numFullSegments * segmentLength;
    vsTrunc             = vs(:, 1:validTimepoints); % Truncate the data to fit whole segments only
    vs_3D               = permute( reshape(vsTrunc, [p1, segmentLength, numFullSegments]), [1 3 2] );
    
    LZc_avg_epoch = zeros(numFullSegments,1);
    LZc_reg_epoch = zeros(numFullSegments,p1);
    for ep = 1:numFullSegments
        vs_2D = squeeze(vs_3D(:,ep,:));
        LZc_reg = zeros(num_reg,1);
        for r = 1:num_reg
            vs_2D_binary_reg = logical(vs_2D(r,:) > mean(vs_2D(r,:)));
            LZc_reg(r) = LZ76(vs_2D_binary_reg)*log2(segmentLength)/segmentLength;
        end
        LZc_avg_epoch(ep) = mean(LZc_reg);
        disp(mean(LZc_reg));
        LZc_reg_epoch(ep,:) = LZc_reg;
        disp(f); 
    end
    LZc_avg_all{f} = LZc_avg_epoch;
    LZc_reg_all{f} = LZc_reg_epoch;
end

%% Combined - Broadband - Epoched - Normalized by phase-shuffling

cd("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Results/SourceRecon/Combined/Truncated/CommonFilter")
data_dir = dir('*.mat');
LZc_avg_all = cell(length(data_dir),1);
LZc_reg_all = cell(length(data_dir),1);

ip  = java.net.InetAddress.getLocalHost.getHostAddress().string;
pctconfig('hostname',ip);
parpool;

parfor f = 1:length(data_dir)

    % Load timeseries (vs)
    vs = importdata(data_dir(f).name);
    
    % reshape vs into a 3D matrix
    [p1, Tm]            = size(vs);
    if f < 26 % not sub034 (who has very short data)
        segmentLength   = 3000;   
    else % sub034
        segmentLength   = 300;
    end    
    numFullSegments     = floor(Tm / segmentLength);
    validTimepoints     = numFullSegments * segmentLength;
    vsTrunc             = vs(:, 1:validTimepoints); % Truncate the data to fit whole segments only
    vs_3D               = permute( reshape(vsTrunc, [p1, segmentLength, numFullSegments]), [1 3 2] );
    
    LZc_avg_epoch = zeros(numFullSegments,1);
    LZc_reg_epoch = zeros(numFullSegments,p1);
    for ep = 1:numFullSegments
        vs_2D = squeeze(vs_3D(:,ep,:));
        LZc_reg = zeros(num_reg,1);
        for r = 1:num_reg
            vs_2D_binary_reg = logical(vs_2D(r,:) > mean(vs_2D(r,:)));
            vs_reg_shuffled = phase_randomize(vs_2D(r,:));
            vs_binary_reg_shuffled = logical(vs_reg_shuffled > mean(vs_reg_shuffled));
            LZc_reg(r) = LZ76(vs_2D_binary_reg)/LZ76(vs_binary_reg_shuffled);
        end
        LZc_avg_epoch(ep) = mean(LZc_reg);
        disp(mean(LZc_reg));
        LZc_reg_epoch(ep,:) = LZc_reg;
        disp(f); 
    end
    LZc_avg_all{f} = LZc_avg_epoch;
    LZc_reg_all{f} = LZc_reg_epoch;
end



%% Statistical analysis - Global LZc

n_perm = 5000;

sub001_NS = LZc_avg_all{11};
sub001_count = [LZc_avg_all{7}; LZc_avg_all{9}];
pval_sub001_count = permutation_test(sub001_NS, sub001_count, n_perm);
sub001_memory = [LZc_avg_all{8}; LZc_avg_all{10}];
pval_sub001_memory = permutation_test(sub001_NS, sub001_memory, n_perm);

sub003_NS = [LZc_avg_all{3}; LZc_avg_all{4}; LZc_avg_all{5}; LZc_avg_all{6}];
sub003_count = LZc_avg_all{1};
pval_sub003_count = permutation_test(sub003_NS, sub003_count, n_perm);
sub003_memory = LZc_avg_all{2};
pval_sub003_memory = permutation_test(sub003_NS, sub003_memory, n_perm);

sub010_NS = [LZc_avg_all{17}; LZc_avg_all{18}; LZc_avg_all{19}];
sub010_count = [LZc_avg_all{12}; LZc_avg_all{14}];
pval_sub010_count = permutation_test(sub010_NS, sub010_count, n_perm);
sub010_memory = [LZc_avg_all{13}; LZc_avg_all{15}];
pval_sub010_memory = permutation_test(sub010_NS, sub010_memory, n_perm);

sub029_NS = [LZc_avg_all{24}; LZc_avg_all{25}];
sub029_count = [LZc_avg_all{20}; LZc_avg_all{22}];
pval_sub029_count = permutation_test(sub029_NS, sub029_count, n_perm);
sub029_memory = [LZc_avg_all{21}; LZc_avg_all{23}];
pval_sub029_memory = permutation_test(sub029_NS, sub029_memory, n_perm);

sub034_NS = [LZc_avg_all{30}; LZc_avg_all{31}];
sub034_count = [LZc_avg_all{26}; LZc_avg_all{28}];
pval_sub034_count = permutation_test(sub034_NS, sub034_count, n_perm);
sub034_memory = [LZc_avg_all{27}; LZc_avg_all{29}];
pval_sub034_memory = permutation_test(sub034_NS, sub034_memory, n_perm);

[~,~,~,adj_pval_sub001] = fdr_bh_2015([pval_sub001_count; pval_sub001_memory]);
[~,~,~,adj_pval_sub003] = fdr_bh_2015([pval_sub003_count; pval_sub003_memory]);
[~,~,~,adj_pval_sub010] = fdr_bh_2015([pval_sub010_count; pval_sub010_memory]);
[~,~,~,adj_pval_sub029] = fdr_bh_2015([pval_sub029_count; pval_sub029_memory]);
[~,~,~,adj_pval_sub034] = fdr_bh_2015([pval_sub034_count; pval_sub034_memory]);

%% Statistical analysis - Regional LZc

n_perm = 5000;
num_reg = 100;
yeo_indices = {[(1:9)'; (51:58)'], [(10:15)'; (59:66)'], [(16:23)'; (67:73)'], [(24:30)'; (74:78)'], [(31:33)'; (79:80)'], [(34:37)'; (81:89)'], [(38:50)'; (90:100)']};

pval_sub001_count_reg = zeros(num_reg,1);
pval_sub001_memory_reg = zeros(num_reg,1);
pval_sub003_count_reg = zeros(num_reg,1);
pval_sub003_memory_reg = zeros(num_reg,1);
pval_sub010_count_reg = zeros(num_reg,1);
pval_sub010_memory_reg = zeros(num_reg,1);
pval_sub029_count_reg = zeros(num_reg,1);
pval_sub029_memory_reg = zeros(num_reg,1);
pval_sub034_count_reg = zeros(num_reg,1);
pval_sub034_memory_reg = zeros(num_reg,1);

for reg = 1:num_reg
    sub001_NS = LZc_reg_all{11}(:,reg);
    sub001_count = [LZc_reg_all{7}(:,reg); LZc_reg_all{9}(:,reg)];
    pval_sub001_count_reg(reg) = permutation_test(sub001_NS, sub001_count, n_perm);
    sub001_memory = [LZc_reg_all{8}(:,reg); LZc_reg_all{10}(:,reg)];
    pval_sub001_memory_reg(reg) = permutation_test(sub001_NS, sub001_memory, n_perm);

    sub003_NS = [LZc_reg_all{3}(:,reg); LZc_reg_all{4}(:,reg); LZc_reg_all{5}(:,reg); LZc_reg_all{6}(:,reg)];
    sub003_count = LZc_reg_all{1}(:,reg);
    pval_sub003_count_reg(reg) = permutation_test(sub003_NS, sub003_count, n_perm);
    sub003_memory = LZc_reg_all{2}(:,reg);
    pval_sub003_memory_reg(reg) = permutation_test(sub003_NS, sub003_memory, n_perm);
    
    sub010_NS = [LZc_reg_all{17}(:,reg); LZc_reg_all{18}(:,reg); LZc_reg_all{19}(:,reg)];
    sub010_count = [LZc_reg_all{12}(:,reg); LZc_reg_all{14}(:,reg)];
    pval_sub010_count_reg(reg) = permutation_test(sub010_NS, sub010_count, n_perm);
    sub010_memory = [LZc_reg_all{13}(:,reg); LZc_reg_all{15}(:,reg)];
    pval_sub010_memory_reg(reg) = permutation_test(sub010_NS, sub010_memory, n_perm);

    sub029_NS = [LZc_reg_all{24}(:,reg); LZc_reg_all{25}(:,reg)];
    sub029_count = [LZc_reg_all{20}(:,reg); LZc_reg_all{22}(:,reg)];
    pval_sub029_count_reg(reg) = permutation_test(sub029_NS, sub029_count, n_perm);
    sub029_memory = [LZc_reg_all{21}(:,reg); LZc_reg_all{23}(:,reg)];
    pval_sub029_memory_reg(reg) = permutation_test(sub029_NS, sub029_memory, n_perm);

    sub034_NS = [LZc_reg_all{30}(:,reg); LZc_reg_all{25}(:,reg)];
    sub034_count = [LZc_reg_all{26}(:,reg); LZc_reg_all{28}(:,reg)];
    pval_sub034_count_reg(reg) = permutation_test(sub034_NS, sub034_count, n_perm);
    sub034_memory = [LZc_reg_all{27}(:,reg); LZc_reg_all{29}(:,reg)];
    pval_sub034_memory_reg(reg) = permutation_test(sub034_NS, sub034_memory, n_perm);
end

[~,~,~,adj_pval_sub001_reg] = fdr_bh_2015([pval_sub001_count_reg; pval_sub001_memory_reg]);
[~,~,~,adj_pval_sub003_reg] = fdr_bh_2015([pval_sub003_count_reg; pval_sub003_memory_reg]);
[~,~,~,adj_pval_sub010_reg] = fdr_bh_2015([pval_sub010_count_reg; pval_sub010_memory_reg]);
[~,~,~,adj_pval_sub029_reg] = fdr_bh_2015([pval_sub029_count_reg; pval_sub029_memory_reg]);
[~,~,~,adj_pval_sub034_reg] = fdr_bh_2015([pval_sub034_count_reg; pval_sub034_memory_reg]);

% Significantly different regions and sign of difference (increase or
% decrease), as well as number of regions per Yeo network that increase or
% decrease

sub001_count_sig_regs = find(adj_pval_sub001_reg(1:100)<0.05);
sub001_count_sig_regs_direction = zeros(length(sub001_count_sig_regs),1);
for reg = 1:length(sub001_count_sig_regs)
    sub001_NS = mean(LZc_reg_all{11}(:,sub001_count_sig_regs(reg)));
    sub001_count = [mean(LZc_reg_all{7}(:,sub001_count_sig_regs(reg))) mean(LZc_reg_all{9}(:,sub001_count_sig_regs(reg)))];
    sub001_count_sig_regs_direction(reg) = sign(mean(sub001_count)-mean(sub001_NS));
end
sub001_count_yeo = zeros(7,1);
sub001_count_yeo_increase = zeros(7,1);
sub001_count_yeo_decrease = zeros(7,1);
for reg = 1:length(sub001_count_sig_regs)
    for y = 1:7
        if ismember(sub001_count_sig_regs(reg), yeo_indices{y})
            sub001_count_yeo(y) = sub001_count_yeo(y) + 1;
            if sub001_count_sig_regs_direction(reg) == 1
                sub001_count_yeo_decrease(y) = sub001_count_yeo_decrease(y) + 1;
            else
                sub001_count_yeo_increase(y) = sub001_count_yeo_increase(y) + 1;
            end
        end
    end
end

sub001_memory_sig_regs = find(adj_pval_sub001_reg(101:200)<0.05);
sub001_memory_sig_regs_direction = zeros(length(sub001_memory_sig_regs),1);
for reg = 1:length(sub001_memory_sig_regs)
    sub001_NS = mean(LZc_reg_all{11}(:,sub001_memory_sig_regs(reg)));
    sub001_memory = [mean(LZc_reg_all{8}(:,sub001_memory_sig_regs(reg))) mean(LZc_reg_all{10}(:,sub001_memory_sig_regs(reg)))];
    sub001_memory_sig_regs_direction(reg) = sign(mean(sub001_memory)-mean(sub001_NS));
end
sub001_memory_yeo = zeros(7,1);
sub001_memory_yeo_increase = zeros(7,1);
sub001_memory_yeo_decrease = zeros(7,1);
for reg = 1:length(sub001_memory_sig_regs)
    for y = 1:7
        if ismember(sub001_memory_sig_regs(reg), yeo_indices{y})
            sub001_memory_yeo(y) = sub001_memory_yeo(y) + 1;
            if sub001_memory_sig_regs_direction(reg) == 1
                sub001_memory_yeo_decrease(y) = sub001_memory_yeo_decrease(y) + 1;
            else
                sub001_memory_yeo_increase(y) = sub001_memory_yeo_increase(y) + 1;
            end
        end
    end
end

sub003_count_sig_regs = find(adj_pval_sub003_reg(1:100)<0.05);
sub003_count_sig_regs_direction = zeros(length(sub003_count_sig_regs),1);
sub003_count_sig_regs_diff = zeros(length(sub003_count_sig_regs),1);
for reg = 1:length(sub003_count_sig_regs)
    sub003_NS = [mean(LZc_reg_all{3}(:,sub003_count_sig_regs(reg))) mean(LZc_reg_all{4}(:,sub003_count_sig_regs(reg))) mean(LZc_reg_all{5}(:,sub003_count_sig_regs(reg))) mean(LZc_reg_all{6}(:,sub003_count_sig_regs(reg)))];
    sub003_count = mean(LZc_reg_all{1}(:,sub003_count_sig_regs(reg)));
    sub003_count_sig_regs_direction(reg) = sign(mean(sub003_count)-mean(sub003_NS));
    sub003_count_sig_regs_diff(reg) = mean(sub003_count)-mean(sub003_NS);
end
sub003_count_yeo = zeros(7,1);
sub003_count_yeo_increase = zeros(7,1);
sub003_count_yeo_decrease = zeros(7,1);
for reg = 1:length(sub003_count_sig_regs)
    for y = 1:7
        if ismember(sub003_count_sig_regs(reg), yeo_indices{y})
            sub003_count_yeo(y) = sub003_count_yeo(y) + 1;
            if sub003_count_sig_regs_direction(reg) == 1
                sub003_count_yeo_decrease(y) = sub003_count_yeo_decrease(y) + 1;
            else
                sub003_count_yeo_increase(y) = sub003_count_yeo_increase(y) + 1;
            end
        end
    end
end

sub003_memory_sig_regs = find(adj_pval_sub003_reg(101:200)<0.05);
sub003_memory_sig_regs_direction = zeros(length(sub003_memory_sig_regs),1);
sub003_memory_sig_regs_diff = zeros(length(sub003_memory_sig_regs),1);
for reg = 1:length(sub003_memory_sig_regs)
    sub003_NS = [mean(LZc_reg_all{3}(:,sub003_memory_sig_regs(reg))) mean(LZc_reg_all{4}(:,sub003_memory_sig_regs(reg))) mean(LZc_reg_all{5}(:,sub003_memory_sig_regs(reg))) mean(LZc_reg_all{6}(:,sub003_memory_sig_regs(reg)))];
    sub003_memory = mean(LZc_reg_all{2}(:,sub003_memory_sig_regs(reg)));
    sub003_memory_sig_regs_direction(reg) = sign(mean(sub003_memory)-mean(sub003_NS));
    sub003_memory_sig_regs_diff(reg) = mean(sub003_memory)-mean(sub003_NS);
end
sub003_memory_yeo = zeros(7,1);
sub003_memory_yeo_increase = zeros(7,1);
sub003_memory_yeo_decrease = zeros(7,1);
for reg = 1:length(sub003_memory_sig_regs)
    for y = 1:7
        if ismember(sub003_memory_sig_regs(reg), yeo_indices{y})
            sub003_memory_yeo(y) = sub003_memory_yeo(y) + 1;
            if sub003_memory_sig_regs_direction(reg) == 1
                sub003_memory_yeo_decrease(y) = sub003_memory_yeo_decrease(y) + 1;
            else
                sub003_memory_yeo_increase(y) = sub003_memory_yeo_increase(y) + 1;
            end
        end
    end
end

sub010_count_sig_regs = find(adj_pval_sub010_reg(1:100)<0.05);
sub010_count_sig_regs_direction = zeros(length(sub010_count_sig_regs),1);
for reg = 1:length(sub010_count_sig_regs)
    sub010_NS = [mean(LZc_reg_all{17}(:,sub010_count_sig_regs(reg))); mean(LZc_reg_all{18}(:,sub010_count_sig_regs(reg))); mean(LZc_reg_all{19}(:,sub010_count_sig_regs(reg)))];
    sub010_count = [mean(LZc_reg_all{12}(:,sub010_count_sig_regs(reg))); mean(LZc_reg_all{14}(:,sub010_count_sig_regs(reg)))];
    sub010_count_sig_regs_direction(reg) = sign(mean(sub010_count)-mean(sub010_NS));
end
sub010_count_yeo = zeros(7,1);
sub010_count_yeo_increase = zeros(7,1);
sub010_count_yeo_decrease = zeros(7,1);
for reg = 1:length(sub010_count_sig_regs)
    for y = 1:7
        if ismember(sub010_count_sig_regs(reg), yeo_indices{y})
            sub010_count_yeo(y) = sub010_count_yeo(y) + 1;
            if sub010_count_sig_regs_direction(reg) == 1
                sub010_count_yeo_decrease(y) = sub010_count_yeo_decrease(y) + 1;
            else
                sub010_count_yeo_increase(y) = sub010_count_yeo_increase(y) + 1;
            end
        end
    end
end

sub010_memory_sig_regs = find(adj_pval_sub010_reg(101:200)<0.05);
sub010_memory_sig_regs_direction = zeros(length(sub010_memory_sig_regs),1);
for reg = 1:length(sub010_memory_sig_regs)
    sub010_NS = [mean(LZc_reg_all{17}(:,sub010_memory_sig_regs(reg))); mean(LZc_reg_all{18}(:,sub010_memory_sig_regs(reg))); mean(LZc_reg_all{19}(:,sub010_memory_sig_regs(reg)))];
    sub010_memory = [mean(LZc_reg_all{13}(:,sub010_memory_sig_regs(reg))); mean(LZc_reg_all{15}(:,sub010_memory_sig_regs(reg)))];
    sub010_memory_sig_regs_direction(reg) = sign(mean(sub010_memory)-mean(sub010_NS));
end
sub010_memory_yeo = zeros(7,1);
sub010_memory_yeo_increase = zeros(7,1);
sub010_memory_yeo_decrease = zeros(7,1);
for reg = 1:length(sub010_memory_sig_regs)
    for y = 1:7
        if ismember(sub010_memory_sig_regs(reg), yeo_indices{y})
            sub010_memory_yeo(y) = sub010_memory_yeo(y) + 1;
            if sub010_memory_sig_regs_direction(reg) == 1
                sub010_memory_yeo_decrease(y) = sub010_memory_yeo_decrease(y) + 1;
            else
                sub010_memory_yeo_increase(y) = sub010_memory_yeo_increase(y) + 1;
            end
        end
    end
end

sub029_count_sig_regs = find(adj_pval_sub029_reg(1:100)<0.05);
sub029_count_sig_regs_direction = zeros(length(sub029_count_sig_regs),1);
for reg = 1:length(sub029_count_sig_regs)
    sub029_NS = [mean(LZc_reg_all{24}(:,sub029_count_sig_regs(reg))); mean(LZc_reg_all{25}(:,sub029_count_sig_regs(reg)))];
    sub029_count = [mean(LZc_reg_all{20}(:,sub029_count_sig_regs(reg))); mean(LZc_reg_all{22}(:,sub029_count_sig_regs(reg)))];
    sub029_count_sig_regs_direction(reg) = sign(mean(sub029_count)-mean(sub029_NS));
end
sub029_count_yeo = zeros(7,1);
sub029_count_yeo_increase = zeros(7,1);
sub029_count_yeo_decrease = zeros(7,1);
for reg = 1:length(sub029_count_sig_regs)
    for y = 1:7
        if ismember(sub029_count_sig_regs(reg), yeo_indices{y})
            sub029_count_yeo(y) = sub029_count_yeo(y) + 1;
            if sub029_count_sig_regs_direction(reg) == 1
                sub029_count_yeo_decrease(y) = sub029_count_yeo_decrease(y) + 1;
            else
                sub029_count_yeo_increase(y) = sub029_count_yeo_increase(y) + 1;
            end
        end
    end
end

sub029_memory_sig_regs = find(adj_pval_sub029_reg(101:200)<0.05);
sub029_memory_sig_regs_direction = zeros(length(sub029_memory_sig_regs),1);
for reg = 1:length(sub029_memory_sig_regs)
    sub029_NS = [mean(LZc_reg_all{24}(:,sub029_memory_sig_regs(reg))); mean(LZc_reg_all{25}(:,sub029_memory_sig_regs(reg)))];
    sub029_memory = [mean(LZc_reg_all{21}(:,sub029_memory_sig_regs(reg))); mean(LZc_reg_all{23}(:,sub029_memory_sig_regs(reg)))];
    sub029_memory_sig_regs_direction(reg) = sign(mean(sub029_memory)-mean(sub029_NS));
end
sub029_memory_yeo = zeros(7,1);
sub029_memory_yeo_increase = zeros(7,1);
sub029_memory_yeo_decrease = zeros(7,1);
for reg = 1:length(sub029_memory_sig_regs)
    for y = 1:7
        if ismember(sub029_memory_sig_regs(reg), yeo_indices{y})
            sub029_memory_yeo(y) = sub029_memory_yeo(y) + 1;
            if sub029_memory_sig_regs_direction(reg) == 1
                sub029_memory_yeo_decrease(y) = sub029_memory_yeo_decrease(y) + 1;
            else
                sub029_memory_yeo_increase(y) = sub029_memory_yeo_increase(y) + 1;
            end
        end
    end
end

sub034_count_sig_regs = find(adj_pval_sub034_reg(1:100)<0.05);
sub034_count_sig_regs_direction = zeros(length(sub034_count_sig_regs),1);
for reg = 1:length(sub034_count_sig_regs)
    sub034_NS = [mean(LZc_reg_all{30}(:,sub034_count_sig_regs(reg))); mean(LZc_reg_all{31}(:,sub034_count_sig_regs(reg)))];
    sub034_count = [mean(LZc_reg_all{26}(:,sub034_count_sig_regs(reg))); mean(LZc_reg_all{28}(:,sub034_count_sig_regs(reg)))];
    sub034_count_sig_regs_direction(reg) = sign(mean(sub034_count)-mean(sub034_NS));
end
sub034_count_yeo = zeros(7,1);
sub034_count_yeo_increase = zeros(7,1);
sub034_count_yeo_decrease = zeros(7,1);
for reg = 1:length(sub034_count_sig_regs)
    for y = 1:7
        if ismember(sub034_count_sig_regs(reg), yeo_indices{y})
            sub034_count_yeo(y) = sub034_count_yeo(y) + 1;
            if sub034_count_sig_regs_direction(reg) == 1
                sub034_count_yeo_decrease(y) = sub034_count_yeo_decrease(y) + 1;
            else
                sub034_count_yeo_increase(y) = sub034_count_yeo_increase(y) + 1;
            end
        end
    end
end

sub034_memory_sig_regs = find(adj_pval_sub034_reg(101:200)<0.05);
sub034_memory_sig_regs_direction = zeros(length(sub034_memory_sig_regs),1);
for reg = 1:length(sub034_memory_sig_regs)
    sub034_NS = [mean(LZc_reg_all{30}(:,sub034_memory_sig_regs(reg))); mean(LZc_reg_all{31}(:,sub034_memory_sig_regs(reg)))];
    sub034_memory = [mean(LZc_reg_all{27}(:,sub034_memory_sig_regs(reg))); mean(LZc_reg_all{29}(:,sub034_memory_sig_regs(reg)))];
    sub034_memory_sig_regs_direction(reg) = sign(mean(sub034_memory)-mean(sub034_NS));
end
sub034_memory_yeo = zeros(7,1);
sub034_memory_yeo_increase = zeros(7,1);
sub034_memory_yeo_decrease = zeros(7,1);
for reg = 1:length(sub034_memory_sig_regs)
    for y = 1:7
        if ismember(sub034_memory_sig_regs(reg), yeo_indices{y})
            sub034_memory_yeo(y) = sub034_memory_yeo(y) + 1;
            if sub034_memory_sig_regs_direction(reg) == 1
                sub034_memory_yeo_decrease(y) = sub034_memory_yeo_decrease(y) + 1;
            else
                sub034_memory_yeo_increase(y) = sub034_memory_yeo_increase(y) + 1;
            end
        end
    end
end


all_sig_regs = {sub001_count_sig_regs, sub001_memory_sig_regs, sub003_count_sig_regs, sub003_memory_sig_regs, sub010_count_sig_regs, sub010_memory_sig_regs, sub029_count_sig_regs, sub029_memory_sig_regs, sub034_count_sig_regs, sub034_memory_sig_regs};
all_sig_regs_direction = {sub001_count_sig_regs_direction, sub001_memory_sig_regs_direction, sub003_count_sig_regs_direction, sub003_memory_sig_regs_direction, sub010_count_sig_regs_direction, sub010_memory_sig_regs_direction, sub029_count_sig_regs_direction, sub029_memory_sig_regs_direction, sub034_count_sig_regs_direction, sub034_memory_sig_regs_direction};
common_sig_regs = all_sig_regs{1};
for i = 2:length(all_sig_regs)
    common_sig_regs = intersect(common_sig_regs, all_sig_regs{i});
end
common_sig_regs_direction = zeros(length(common_sig_regs),length(all_sig_regs));
for i = 1:length(common_sig_regs)
    for j = 1:length(all_sig_regs)
        reg_idx = find(all_sig_regs{j}==common_sig_regs(i));
        common_sig_regs_direction(i,j) = all_sig_regs_direction{j}(reg_idx);
    end
end

all_values = vertcat(all_sig_regs{:});
[unique_vals, ~, idx] = unique(all_values);
counts = accumarray(idx, 1);
[counts_sorted, sort_idx] = sort(counts, 'descend');
most_common_vals = unique_vals(sort_idx);

most_common_vals_direction = zeros(length(most_common_vals),length(all_sig_regs));
for i = 1:length(most_common_vals)
    for j = 1:length(all_sig_regs)
        reg_idx = find(all_sig_regs{j}==most_common_vals(i));
        if isempty(reg_idx)
            most_common_vals_direction(i,j) = 0;
        else
            most_common_vals_direction(i,j) = all_sig_regs_direction{j}(reg_idx);
        end
    end
end

%% Plot regional LZc 

addpath(genpath("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Scripts"))

sub001_count = mean([LZc_reg_all{7}; LZc_reg_all{9}]);
sub001_memory = mean([LZc_reg_all{8}; LZc_reg_all{10}]); 
sub001_NS = mean(LZc_reg_all{11});

sub001_count_diff = sub001_NS-sub001_count;
sub001_count_diff_sig = zeros(num_reg,1);
for r = 1:num_reg
    if pval_sub001_count_reg(r) < 0.05
        sub001_count_diff_sig(r) = sub001_count_diff(r);
    end
end
sub001_memory_diff = sub001_NS-sub001_memory;
sub001_memory_diff_sig = zeros(num_reg,1);
for r = 1:num_reg
    if pval_sub001_memory_reg(r) < 0.05
        sub001_memory_diff_sig(r) = sub001_memory_diff(r);
    end
end

sub003_count = mean(LZc_reg_all{1});
sub003_memory = mean(LZc_reg_all{2});
sub003_NS = mean([LZc_reg_all{3}; LZc_reg_all{4}; LZc_reg_all{5}; LZc_reg_all{6}]);

sub003_count_diff = sub003_NS-sub003_count;
sub003_count_diff_sig = zeros(num_reg,1);
for r = 1:num_reg
    if pval_sub003_count_reg(r) < 0.05
        sub003_count_diff_sig(r) = sub003_count_diff(r);
    end
end
sub003_memory_diff = sub003_NS-sub003_memory;
sub003_memory_diff_sig = zeros(num_reg,1);
for r = 1:num_reg
    if pval_sub003_memory_reg(r) < 0.05
        sub003_memory_diff_sig(r) = sub003_memory_diff(r);
    end
end

sub010_count = mean([LZc_reg_all{12}; LZc_reg_all{14}]);
sub010_memory = mean([LZc_reg_all{13}; LZc_reg_all{15}]); 
sub010_NS = mean([LZc_reg_all{17}; LZc_reg_all{18}; LZc_reg_all{19}]); 

sub010_count_diff = sub010_NS-sub010_count;
sub010_count_diff_sig = zeros(num_reg,1);
for r = 1:num_reg
    if pval_sub010_count_reg(r) < 0.05
        sub010_count_diff_sig(r) = sub010_count_diff(r);
    end
end
sub010_memory_diff = sub010_NS-sub010_memory;
sub010_memory_diff_sig = zeros(num_reg,1);
for r = 1:num_reg
    if pval_sub010_memory_reg(r) < 0.05
        sub010_memory_diff_sig(r) = sub010_memory_diff(r);
    end
end

sub029_count = mean([LZc_reg_all{20}; LZc_reg_all{22}]);
sub029_memory = mean([LZc_reg_all{21}; LZc_reg_all{23}]); 
sub029_NS = mean([LZc_reg_all{24}; LZc_reg_all{25}]); 

sub029_count_diff = sub029_NS-sub029_count;
sub029_count_diff_sig = zeros(num_reg,1);
for r = 1:num_reg
    if pval_sub029_count_reg(r) < 0.05
        sub029_count_diff_sig(r) = sub029_count_diff(r);
    end
end
sub029_memory_diff = sub029_NS-sub029_memory;
sub029_memory_diff_sig = zeros(num_reg,1);
for r = 1:num_reg
    if pval_sub029_memory_reg(r) < 0.05
        sub029_memory_diff_sig(r) = sub029_memory_diff(r);
    end
end

sub034_count = mean([LZc_reg_all{26}; LZc_reg_all{28}]);
sub034_memory = mean([LZc_reg_all{27}; LZc_reg_all{29}]); 
sub034_NS = mean([LZc_reg_all{30}; LZc_reg_all{31}]); 

sub034_count_diff = sub034_NS-sub034_count;
sub034_count_diff_sig = zeros(num_reg,1);
for r = 1:num_reg
    if pval_sub034_count_reg(r) < 0.05
        sub034_count_diff_sig(r) = sub034_count_diff(r);
    end
end
sub034_memory_diff = sub034_NS-sub034_memory;
sub034_memory_diff_sig = zeros(num_reg,1);
for r = 1:num_reg
    if pval_sub034_memory_reg(r) < 0.05
        sub034_memory_diff_sig(r) = sub034_memory_diff(r);
    end
end

min_val = min(min([sub001_count_diff; sub001_memory_diff; sub003_count_diff; sub003_memory_diff; sub010_count_diff; sub010_memory_diff; sub029_count_diff; sub029_memory_diff; sub034_count_diff; sub034_memory_diff]));
max_val = max(max([sub001_count_diff; sub001_memory_diff; sub003_count_diff; sub003_memory_diff; sub010_count_diff; sub010_memory_diff; sub029_count_diff; sub029_memory_diff; sub034_count_diff; sub034_memory_diff]));

%% Plot regional LZc - sub001

sub001_count_surf = parcel_to_surface(sub001_count_diff_sig, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub001_count_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub001_memory_surf = parcel_to_surface(sub001_memory_diff_sig, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub001_memory_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')


%% Plot regional LZc - sub003

sub003_count_surf = parcel_to_surface(sub003_count_diff_sig, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub003_count_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub003_memory_surf = parcel_to_surface(sub003_memory_diff_sig, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub003_memory_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

%% Plot regional LZc - sub010

sub010_count_surf = parcel_to_surface(sub010_count_diff_sig, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub010_count_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub010_memory_surf = parcel_to_surface(sub010_memory_diff_sig, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub010_memory_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

%% Plot regional LZc - sub029

sub029_count_surf = parcel_to_surface(sub029_count_diff_sig, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub029_count_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub029_memory_surf = parcel_to_surface(sub029_memory_diff_sig, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub029_memory_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

%% Plot regional LZc - sub034

sub034_count_surf = parcel_to_surface(sub034_count_diff_sig, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub034_count_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub034_memory_surf = parcel_to_surface(sub034_memory_diff_sig, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub034_memory_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')



%% Plot broadband LZc 

subject_labels   = {'sub001', 'sub003', 'sub010', 'sub029', 'sub034'};
condition_labels = {'Counting', 'Memory', 'NS'};

% Column 1: subject index; column 2: condition index
map = [2 1; 2 2; 2 3; 2 3; 2 3; 2 3;
       1 1; 1 2; 1 1; 1 2; 1 3;
       3 1; 3 2; 3 1; 3 2; 3 3; 3 3; 3 3; 3 3;
       4 1; 4 2; 4 1; 4 2; 4 3; 4 3;
       5 1; 5 2; 5 1; 5 2; 5 3; 5 3];  

subjects   = {};
conditions = {};
all_data   = [];

for i = 1:length(LZc_avg_all)
    d = LZc_avg_all{i};
    all_data = [all_data; d];
    s = subject_labels{map(i,1)};
    c = condition_labels{map(i,2)};
    subjects   = [subjects; repmat({s}, numel(d), 1)];
    conditions = [conditions; repmat({c}, numel(d), 1)];
end

figure;
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
ylabel('Lempel-Ziv Complexity')
title('Lempel-Ziv Complexity (Broadband)')
legend(conditionNames, 'Location', 'northeast')

hold on
offsets = [-0.2, 0, 0.2];  % Left = Counting, Center = Memory, Right = NS
% Significance bars for selected comparisons
y_offset = 0.03;  % vertical space between bars
fontSize = 14;    % size of asterisk labels

% Mapping from condition to offset
condOffset = containers.Map(condition_labels, offsets);

% Significance comparisons: {subject_x, cond1, cond2, stars}
comparisons = {
    % 2, 'Counting', 'NS', '***'; % sub003: NS vs counting
    % 3, 'Counting', 'NS', '**'; % sub010: NS vs counting
    % 3, 'Memory', 'NS', '*'; % sub010: NS vs memory
    % 4, 'Counting', 'NS', '***'; % sub029: NS vs counting
    % 4, 'Memory', 'NS', '***'; % sub029: NS vs memory
    % 5, 'Counting', 'NS', '***'; % sub034: NS vs counting
    % 5, 'Memory', 'NS', '***'; % sub034: NS vs memory
    1, 'Counting', 'NS', '***'; % sub001: NS vs counting
    4, 'Counting', 'NS', '***'; % sub029: NS vs counting
    4, 'Memory', 'NS', '***'; % sub029: NS vs memory
};

% Initialize count of significance bars per subject
subject_sig_count = containers.Map(subject_labels, zeros(1, numel(subject_labels)));

for i = 1:size(comparisons, 1)
    subjX = comparisons{i,1};
    cond1 = comparisons{i,2};
    cond2 = comparisons{i};
    stars = comparisons{i,4};

    % x-locations (add offset for condition)
    x1 = subjX + condOffset(cond1);
    x2 = subjX + condOffset(cond2);

    % Find maximum y for this subject
    subjName = subject_labels{subjX};
    subj_data = all_data(strcmp(subjects, subjName));
    base_y = max(subj_data);

    % How many sig bars already drawn for this subject?
    current_count = subject_sig_count(subjName);
    
    % Compute y-position: base + a small offset per bar
    y1 = base_y + (current_count + 1) * y_offset;

    % Draw horizontal line
    h = plot([x1 x2], [y1 y1], 'k', 'LineWidth', 1.2);
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';

    % Draw asterisks
    text(mean([x1 x2]), y1 + y_offset * 0.2, stars, ...
        'HorizontalAlignment', 'center', 'FontSize', fontSize);

    % Update count for this subject
    subject_sig_count(subjName) = current_count + 1;
end

