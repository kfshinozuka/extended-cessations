%% Phi*

addpath("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Scripts")
addpath("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Scripts/PE")
cd("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Results/SourceRecon/Combined/Truncated/CommonFilter")
data_dir = dir('*.mat');
phi_all = cell(length(data_dir),1);
opt_tau = 100;

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
    
    phi_epoch_avg = zeros(numFullSegments,opt_tau);
    for tau = 1:opt_tau
        for ep = 1:numFullSegments
            vs_2D = squeeze(vs_3D(:,ep,:));  
            phi_epoch = ones(p1);
            for r1 = 1:p1
                for r2 = r1+1:p1
                    phi_epoch(r1,r2) = compute_phi_star_simple(vs_2D,r1,r2,tau);
                    phi_epoch(r2,r1) = phi_epoch(r1,r2);
                end
            end
            phi_epoch_avg(ep,tau) = mean(mean(phi_epoch));
            %disp(phi_epoch_avg(ep,tau));
        end
    end
    disp(f);
    phi_all{f} = phi_epoch_avg;
end

%% Statistics

addpath('/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Scripts')

n_perm = 5000;
opt_tau = 1;

sub001_NS = phi_all{11}(:,opt_tau);
sub001_count = [phi_all{7}(:,opt_tau); phi_all{9}(:,opt_tau)];
pval_sub001_count = permutation_test(sub001_NS, sub001_count, n_perm);
sub001_memory = [phi_all{8}(:,opt_tau); phi_all{10}(:,opt_tau)];
pval_sub001_memory = permutation_test(sub001_NS, sub001_memory, n_perm);

sub003_NS = [phi_all{3}(:,opt_tau); phi_all{4}(:,opt_tau); phi_all{5}(:,opt_tau); phi_all{6}(:,opt_tau)];
sub003_count = phi_all{1}(:,opt_tau);
pval_sub003_count = permutation_test(sub003_NS, sub003_count, n_perm);
sub003_memory = phi_all{2}(:,opt_tau);
pval_sub003_memory = permutation_test(sub003_NS, sub003_memory, n_perm);

sub010_NS = [phi_all{17}(:,opt_tau); phi_all{18}(:,opt_tau); phi_all{19}(:,opt_tau)];
sub010_count = [phi_all{12}(:,opt_tau); phi_all{14}(:,opt_tau)];
pval_sub010_count = permutation_test(sub010_NS, sub010_count, n_perm);
sub010_memory = [phi_all{13}(:,opt_tau); phi_all{15}(:,opt_tau)];
pval_sub010_memory = permutation_test(sub010_NS, sub010_memory, n_perm);

sub029_NS = [phi_all{24}(:,opt_tau); phi_all{25}(:,opt_tau)];
sub029_count = [phi_all{20}(:,opt_tau); phi_all{22}(:,opt_tau)];
pval_sub029_count = permutation_test(sub029_NS, sub029_count, n_perm);
sub029_memory = [phi_all{21}(:,opt_tau); phi_all{23}(:,opt_tau)];
pval_sub029_memory = permutation_test(sub029_NS, sub029_memory, n_perm);

sub034_NS = [phi_all{30}(:,opt_tau); phi_all{31}(:,opt_tau)];
sub034_count = [phi_all{26}(:,opt_tau); phi_all{28}(:,opt_tau)];
pval_sub034_count = permutation_test(sub034_NS, sub034_count, n_perm);
sub034_memory = [phi_all{27}(:,opt_tau); phi_all{29}(:,opt_tau)];
pval_sub034_memory = permutation_test(sub034_NS, sub034_memory, n_perm);

[~,~,~,adj_pval_sub001] = fdr_bh_2015([pval_sub001_count; pval_sub001_memory]);
[~,~,~,adj_pval_sub003] = fdr_bh_2015([pval_sub003_count; pval_sub003_memory]);
[~,~,~,adj_pval_sub010] = fdr_bh_2015([pval_sub010_count; pval_sub010_memory]);
[~,~,~,adj_pval_sub029] = fdr_bh_2015([pval_sub029_count; pval_sub029_memory]);
[~,~,~,adj_pval_sub034] = fdr_bh_2015([pval_sub034_count; pval_sub034_memory]);

%% Plot phi

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

for i = 1:length(phi_all)
    d = phi_all{i}(:,opt_tau);
    all_data = [all_data; d];
    s = subject_labels{map(i,1)};
    c = condition_labels{map(i,2)};
    subjects   = [subjects; repmat({s}, numel(d), 1)];
    conditions = [conditions; repmat({c}, numel(d), 1)];
end

comparisons = {
    1, 'Counting', 'NS', '***';
    1, 'Memory', 'NS', '***';
    2, 'Counting', 'NS', '**';   % sub003: NS vs Counting
    2, 'Memory', 'NS', '**';
    3, 'Counting', 'NS', '***';   % sub010: NS vs Counting
    4, 'Memory', 'NS', '***';
    5, 'Counting', 'NS', '***';
    5, 'Memory',   'NS', '*';   % sub034: NS vs Memory
};



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
ylabel('Phi*')
title('Phi*')
legend(conditionNames, 'Location', 'northeast')

hold on
offsets = [-0.2, 0, 0.2];  % Left = Counting, Center = Memory, Right = NS
% Significance bars for selected comparisons
y_offset = 0.03;  % vertical space between bars
fontSize = 14;    % size of asterisk labels

% Mapping from condition to offset
condOffset = containers.Map(condition_labels, offsets);

% Initialize count of significance bars per subject
subject_sig_count = containers.Map(subject_labels, zeros(1, numel(subject_labels)));

for i = 1:size(comparisons, 1)
    subjX = comparisons{i,1};
    cond1 = comparisons{i,2};
    cond2 = comparisons{i,3};
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

