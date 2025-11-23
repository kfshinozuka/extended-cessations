%% Global coherence - source-space

addpath("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Scripts")
datapath = '/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Results/SourceRecon/Combined/Truncated/CommonFilter/New';
cd(datapath)
d_combined = dir('*.mat');

fs = 200;
winSize = 200;
overlap = 100;
GC = cell(length(d_combined),1);
LE = cell(length(d_combined),1);

ip  = java.net.InetAddress.getLocalHost.getHostAddress().string;
pctconfig('hostname',ip);
parpool;

parfor f = 1:length(d_combined) % loop through datasets

    vs = importdata(d_combined(f).name); % source-reconstructed timeseries
    if f < 26 % not sub034 (who has very short data)
        segmentLength = 3000;
    else % sub034
        segmentLength = 300;
    end
    [GC{f}, LE{f}] = computeGlobalCoherence_epoch(vs,fs,winSize,overlap,segmentLength);
    disp(GC{f})
    fprintf('finished for file %d \n',f)       
    
end 



%% Plot GC - sub001

sub001_GC_counting = [GC{7}(1:51,:) GC{9}(1:51,:)];
sub001_meanGC_counting = mean(sub001_GC_counting,2);
sub001_stdGC_counting = std(sub001_GC_counting,0,2);
sub001_SE_counting = sub001_stdGC_counting / sqrt(size(sub001_GC_counting,2));

sub001_GC_memory = [GC{8}(1:51,:) GC{10}(1:51,:)]; 
sub001_meanGC_memory = mean(sub001_GC_memory,2);
sub001_stdGC_memory = std(sub001_GC_memory,0,2);
sub001_SE_memory = sub001_stdGC_memory / sqrt(size(sub001_GC_memory,2));

sub001_GC_NS = GC{11}(1:51,:);
sub001_meanGC_NS = mean(sub001_GC_NS,2);
sub001_stdGC_NS = std(sub001_GC_NS,0,2);
sub001_SE_NS = sub001_stdGC_NS / sqrt(size(sub001_GC_NS,2));

% Plot the result with standard error shading
figure;
hold on;

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub001_meanGC_counting' + sub001_SE_counting', fliplr(sub001_meanGC_counting' - sub001_SE_counting')], 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p1 = plot(0:50, sub001_meanGC_counting, 'LineWidth', 2, 'Color', 'r');
hold on

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub001_meanGC_memory' + sub001_SE_memory', fliplr(sub001_meanGC_memory' - sub001_SE_memory')], [0.85 0.325 0.098], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p2 = plot((0:50), sub001_meanGC_memory, 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);
hold on

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub001_meanGC_NS' + sub001_SE_NS', fliplr(sub001_meanGC_NS' - sub001_SE_NS')], 'b', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p3 = plot((0:50), sub001_meanGC_NS, 'LineWidth', 2, 'Color', 'b');
hold on

legend([p1,p2,p3],{'Counting','Memory','NS'},'Location','NorthEast')
xlabel('Frequency (Hz)');
ylabel('Global Coherence (%)');
title('sub001 Global Coherence (with SE)');
grid on;
hold off;

%% Plot GC - sub003

sub003_GC_counting = GC{1}(1:51,:);
sub003_meanGC_counting = mean(sub003_GC_counting,2);
sub003_stdGC_counting = std(sub003_GC_counting,0,2);
sub003_SE_counting = sub003_stdGC_counting / sqrt(size(sub003_GC_counting,2));

sub003_GC_memory = GC{2}(1:51,:);
sub003_meanGC_memory = mean(sub003_GC_memory,2);
sub003_stdGC_memory = std(sub003_GC_memory,0,2);
sub003_SE_memory = sub003_stdGC_memory / sqrt(size(sub003_GC_memory,2));

sub003_GC_NS = [GC{3}(1:51,:) GC{4}(1:51,:) GC{5}(1:51,:) GC{6}(1:51,:)];
sub003_meanGC_NS = mean(sub003_GC_NS,2);
sub003_stdGC_NS = std(sub003_GC_NS,0,2);
sub003_SE_NS = sub003_stdGC_NS / sqrt(size(sub003_GC_NS,2));

% Plot the result with standard error shading
figure;
hold on;

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub003_meanGC_counting' + sub003_SE_counting', fliplr(sub003_meanGC_counting' - sub003_SE_counting')], 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p1 = plot(0:50, sub003_meanGC_counting, 'LineWidth', 2, 'Color', 'r');
hold on

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub003_meanGC_memory' + sub003_SE_memory', fliplr(sub003_meanGC_memory' - sub003_SE_memory')], [0.85 0.325 0.098], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p2 = plot((0:50), sub003_meanGC_memory, 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);
hold on

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub003_meanGC_NS' + sub003_SE_NS', fliplr(sub003_meanGC_NS' - sub003_SE_NS')], 'b', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p3 = plot((0:50), sub003_meanGC_NS, 'LineWidth', 2, 'Color', 'b');
hold on

legend([p1,p2,p3],{'Counting','Memory','NS'},'Location','NorthEast')
xlabel('Frequency (Hz)');
ylabel('Global Coherence (%)');
title('sub003 Global Coherence (with SE)');
grid on;
hold off;

%% Plot GC - sub010

sub010_GC_counting = [GC{12}(1:51,:) GC{14}(1:51,:)];
sub010_meanGC_counting = mean(sub010_GC_counting,2);
sub010_stdGC_counting = std(sub010_GC_counting,0,2);
sub010_SE_counting = sub010_stdGC_counting / sqrt(size(sub010_GC_counting,2));

sub010_GC_memory = [GC{13}(1:51,:) GC{15}(1:51,:)]; 
sub010_meanGC_memory = mean(sub010_GC_memory,2);
sub010_stdGC_memory = std(sub010_GC_memory,0,2);
sub010_SE_memory = sub010_stdGC_memory / sqrt(size(sub010_GC_memory,2));

sub010_GC_NS = [GC{17}(1:51,:) GC{18}(1:51,:) GC{19}(1:51,:)]; 
sub010_meanGC_NS = mean(sub010_GC_NS,2);
sub010_stdGC_NS = std(sub010_GC_NS,0,2);
sub010_SE_NS = sub010_stdGC_NS / sqrt(size(sub010_GC_NS,2));

% Plot the result with standard error shading
figure;
hold on;

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub010_meanGC_counting' + sub010_SE_counting', fliplr(sub010_meanGC_counting' - sub010_SE_counting')], 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p1 = plot(0:50, sub010_meanGC_counting, 'LineWidth', 2, 'Color', 'r');
hold on

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub010_meanGC_memory' + sub010_SE_memory', fliplr(sub010_meanGC_memory' - sub010_SE_memory')], [0.85 0.325 0.098], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p2 = plot((0:50), sub010_meanGC_memory, 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);
hold on

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub010_meanGC_NS' + sub010_SE_NS', fliplr(sub010_meanGC_NS' - sub010_SE_NS')], 'b', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p3 = plot((0:50), sub010_meanGC_NS, 'LineWidth', 2, 'Color', 'b');
hold on

legend([p1,p2,p3],{'Counting','Memory','NS'},'Location','NorthEast')
xlabel('Frequency (Hz)');
ylabel('Global Coherence (%)');
title('sub010 Global Coherence (with SE)');
grid on;
hold off;

%% Plot GC - sub029

sub029_GC_counting = [GC{20}(1:51,:) GC{22}(1:51,:)];
sub029_meanGC_counting = mean(sub029_GC_counting,2);
sub029_stdGC_counting = std(sub029_GC_counting,0,2);
sub029_SE_counting = sub029_stdGC_counting / sqrt(size(sub029_GC_counting,2));

sub029_GC_memory = [GC{21}(1:51,:) GC{23}(1:51,:)]; 
sub029_meanGC_memory = mean(sub029_GC_memory,2);
sub029_stdGC_memory = std(sub029_GC_memory,0,2);
sub029_SE_memory = sub029_stdGC_memory / sqrt(size(sub029_GC_memory,2));

sub029_GC_NS = [GC{24}(1:51,:) GC{25}(1:51,:)]; 
sub029_meanGC_NS = mean(sub029_GC_NS,2);
sub029_stdGC_NS = std(sub029_GC_NS,0,2);
sub029_SE_NS = sub029_stdGC_NS / sqrt(size(sub029_GC_NS,2));

% Plot the result with standard error shading
figure;
hold on;

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub029_meanGC_counting' + sub029_SE_counting', fliplr(sub029_meanGC_counting' - sub029_SE_counting')], 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p1 = plot(0:50, sub029_meanGC_counting, 'LineWidth', 2, 'Color', 'r');
hold on

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub029_meanGC_memory' + sub029_SE_memory', fliplr(sub029_meanGC_memory' - sub029_SE_memory')], [0.85 0.325 0.098], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p2 = plot((0:50), sub029_meanGC_memory, 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);
hold on

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub029_meanGC_NS' + sub029_SE_NS', fliplr(sub029_meanGC_NS' - sub029_SE_NS')], 'b', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p3 = plot((0:50), sub029_meanGC_NS, 'LineWidth', 2, 'Color', 'b');
hold on

legend([p1,p2,p3],{'Counting','Memory','NS'},'Location','NorthEast')
xlabel('Frequency (Hz)');
ylabel('Global Coherence (%)');
title('sub029 Global Coherence (with SE)');
grid on;
hold off;

%% Plot GC - sub034

sub034_GC_counting = [GC{26}(1:51,:) GC{28}(1:51,:)];
sub034_meanGC_counting = mean(sub034_GC_counting,2);
sub034_stdGC_counting = std(sub034_GC_counting,0,2);
sub034_SE_counting = sub034_stdGC_counting / sqrt(size(sub034_GC_counting,2));

sub034_GC_memory = [GC{27}(1:51,:) GC{29}(1:51,:)]; 
sub034_meanGC_memory = mean(sub034_GC_memory,2);
sub034_stdGC_memory = std(sub034_GC_memory,0,2);
sub034_SE_memory = sub034_stdGC_memory / sqrt(size(sub034_GC_memory,2));

sub034_GC_NS = [GC{30}(1:51,:) GC{31}(1:51,:)]; 
sub034_meanGC_NS = mean(sub034_GC_NS,2);
sub034_stdGC_NS = std(sub034_GC_NS,0,2);
sub034_SE_NS = sub034_stdGC_NS / sqrt(size(sub034_GC_NS,2));

% Plot the result with standard error shading
figure;
hold on;

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub034_meanGC_counting' + sub034_SE_counting', fliplr(sub034_meanGC_counting' - sub034_SE_counting')], 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p1 = plot(0:50, sub034_meanGC_counting, 'LineWidth', 2, 'Color', 'r');
hold on

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub034_meanGC_memory' + sub034_SE_memory', fliplr(sub034_meanGC_memory' - sub034_SE_memory')], [0.85 0.325 0.098], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p2 = plot((0:50), sub034_meanGC_memory, 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);
hold on

% Shaded area (mean ± SE)
fill([(0:50), fliplr((0:50))], [sub034_meanGC_NS' + sub034_SE_NS', fliplr(sub034_meanGC_NS' - sub034_SE_NS')], 'b', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean GC spectrum
p3 = plot((0:50), sub034_meanGC_NS, 'LineWidth', 2, 'Color', 'b');
hold on

legend([p1,p2,p3],{'Counting','Memory','NS'},'Location','NorthEast')
xlabel('Frequency (Hz)');
ylabel('Global Coherence (%)');
title('sub034 Global Coherence (with SE)');
grid on;
hold off;



%% Statistics - within-subject

n_perm = 5000;
pval_sub001_count = zeros(50,1);
pval_sub001_memory = zeros(50,1);
pval_sub003_count = zeros(50,1);
pval_sub003_memory = zeros(50,1);
pval_sub010_count = zeros(50,1);
pval_sub010_memory = zeros(50,1);
pval_sub029_count = zeros(50,1);
pval_sub029_memory = zeros(50,1);
pval_sub034_count = zeros(50,1);
pval_sub034_memory = zeros(50,1);

for fr = 2:51

    sub001_NS = GC{11}(fr,:)';
    sub001_count = [GC{7}(fr,:)'; GC{9}(fr,:)'];
    pval_sub001_count(fr-1) = permutation_test(sub001_NS, sub001_count, n_perm);
    sub001_memory = [GC{8}(fr,:)'; GC{10}(fr,:)'];
    pval_sub001_memory(fr-1) = permutation_test(sub001_NS, sub001_memory, n_perm);

    sub003_NS = [GC{3}(fr,:)'; GC{4}(fr,:)'; GC{5}(fr,:)'; GC{6}(fr,:)'];
    sub003_count = GC{1}(fr,:)';
    pval_sub003_count(fr-1) = permutation_test(sub003_NS, sub003_count, n_perm);
    sub003_memory = GC{2}(fr,:)';
    pval_sub003_memory(fr-1) = permutation_test(sub003_NS, sub003_memory, n_perm);
    
    sub010_NS = [GC{17}(fr,:)'; GC{18}(fr,:)'; GC{19}(fr,:)'];
    sub010_count = [GC{12}(fr,:)'; GC{14}(fr,:)'];
    pval_sub010_count(fr-1) = permutation_test(sub010_NS, sub010_count, n_perm);
    sub010_memory = [GC{13}(fr,:)'; GC{15}(fr,:)'];
    pval_sub010_memory(fr-1) = permutation_test(sub010_NS, sub010_memory, n_perm);

    sub029_NS = [GC{24}(fr,:)'; GC{25}(fr,:)'];
    sub029_count = [GC{20}(fr,:)'; GC{22}(fr,:)'];
    pval_sub029_count(fr-1) = permutation_test(sub029_NS, sub029_count, n_perm);
    sub029_memory = [GC{21}(fr,:)'; GC{23}(fr,:)'];
    pval_sub029_memory(fr-1) = permutation_test(sub029_NS, sub029_memory, n_perm);

    sub034_NS = [GC{30}(fr,:)'; GC{31}(fr,:)'];
    sub034_count = [GC{26}(fr,:)'; GC{28}(fr,:)'];
    pval_sub034_count(fr-1) = permutation_test(sub034_NS, sub034_count, n_perm);
    sub034_memory = [GC{27}(fr,:)'; GC{29}(fr,:)'];
    pval_sub034_memory(fr-1) = permutation_test(sub034_NS, sub034_memory, n_perm);
end

[~, ~, ~, adj_pval_sub001] = fdr_bh_2015([pval_sub001_count; pval_sub001_memory]);
[~, ~, ~, adj_pval_sub003] = fdr_bh_2015([pval_sub003_count; pval_sub003_memory]);
[~, ~, ~, adj_pval_sub010] = fdr_bh_2015([pval_sub010_count; pval_sub010_memory]);
[~, ~, ~, adj_pval_sub029] = fdr_bh_2015([pval_sub029_count; pval_sub029_memory]);
[~, ~, ~, adj_pval_sub034] = fdr_bh_2015([pval_sub034_count; pval_sub034_memory]);

sub001_NSvscount = zeros(length(find(adj_pval_sub001(1:50)<0.05)),1);
sub001_count_sig = find(adj_pval_sub001(1:50)<0.05);
for i = 1:length(sub001_NSvscount)
    sub001_NS = GC{11}(sub001_count_sig(i)+1,:);
    sub001_count = [GC{7}(sub001_count_sig(i)+1,:) GC{9}(sub001_count_sig(i)+1,:)];
    sub001_NSvscount(i) = mean(sub001_NS)-mean(sub001_count);
end

sub001_NSvsmemory = zeros(length(find(adj_pval_sub001(51:100)<0.05)),1);
sub001_memory_sig = find(adj_pval_sub001(51:100)<0.05);
for i = 1:length(sub001_NSvsmemory)
    sub001_NS = GC{11}(sub001_memory_sig(i)+1,:);
    sub001_memory = [GC{8}(sub001_memory_sig(i)+1,:) GC{10}(sub001_memory_sig(i)+1,:)];
    sub001_NSvsmemory(i) = mean(sub001_NS)-mean(sub001_memory);
end

sub003_NSvscount = zeros(length(find(adj_pval_sub003(1:50)<0.05)),1);
sub003_count_sig = find(adj_pval_sub003(1:50)<0.05);
for i = 1:length(sub003_NSvscount)
    sub003_NS = [GC{3}(sub003_count_sig(i)+1,:) GC{4}(sub003_count_sig(i)+1,:) GC{5}(sub003_count_sig(i)+1,:) GC{6}(sub003_count_sig(i)+1,:)];
    sub003_count = GC{1}(sub003_count_sig(i)+1,:);
    sub003_NSvscount(i) = mean(sub003_NS)-mean(sub003_count);
end

sub003_NSvsmemory = zeros(length(find(adj_pval_sub003(51:100)<0.05)),1);
sub003_memory_sig = find(adj_pval_sub003(51:100)<0.05);
for i = 1:length(sub003_NSvsmemory)
    sub003_NS = [GC{3}(sub003_memory_sig(i)+1,:); GC{4}(sub003_memory_sig(i)+1,:); GC{5}(sub003_memory_sig(i)+1,:); GC{6}(sub003_memory_sig(i)+1,:)];
    sub003_memory = GC{2}(sub003_memory_sig(i)+1,:);
    sub003_NSvsmemory(i) = mean(sub003_NS)-mean(sub003_memory);
end

sub010_NSvscount = zeros(length(find(adj_pval_sub010(1:50)<0.05)),1);
sub010_count_sig = find(adj_pval_sub010(1:50)<0.05);
for i = 1:length(sub010_NSvscount)
    sub010_NS = [GC{17}(sub010_count_sig(i)+1,:) GC{18}(sub010_count_sig(i)+1,:) GC{19}(sub010_count_sig(i)+1,:)];
    sub010_count = [GC{12}(sub010_count_sig(i)+1,:) GC{14}(sub010_count_sig(i)+1,:)];
    sub010_NSvscount(i) = mean(sub010_NS)-mean(sub010_count);
end

sub010_NSvsmemory = zeros(length(find(adj_pval_sub010(51:100)<0.05)),1);
sub010_memory_sig = find(adj_pval_sub010(51:100)<0.05);
for i = 1:length(sub010_NSvsmemory)
    sub010_NS = [GC{17}(sub010_memory_sig(i)+1,:) GC{18}(sub010_memory_sig(i)+1,:) GC{19}(sub010_memory_sig(i)+1,:)];
    sub010_memory = [GC{13}(sub010_memory_sig(i)+1,:) GC{15}(sub010_memory_sig(i)+1,:)];
    sub010_NSvsmemory(i) = mean(sub010_NS)-mean(sub010_memory);
end

sub029_NSvscount = zeros(length(find(adj_pval_sub029(1:50)<0.05)),1);
sub029_count_sig = find(adj_pval_sub029(1:50)<0.05);
for i = 1:length(sub029_NSvscount)
    sub029_NS = [GC{24}(sub029_count_sig(i)+1,:) GC{25}(sub029_count_sig(i)+1,:)];
    sub029_count = [GC{20}(sub029_count_sig(i)+1,:) GC{22}(sub029_count_sig(i)+1,:)];
    sub029_NSvscount(i) = mean(sub029_NS)-mean(sub029_count);
end

sub029_NSvsmemory = zeros(length(find(adj_pval_sub029(51:100)<0.05)),1);
sub029_memory_sig = find(adj_pval_sub029(51:100)<0.05);
for i = 1:length(sub029_NSvsmemory)
    sub029_NS = [GC{24}(sub029_memory_sig(i)+1,:) GC{25}(sub029_memory_sig(i)+1,:)];
    sub029_memory = [GC{21}(sub029_memory_sig(i)+1,:) GC{23}(sub029_memory_sig(i)+1,:)];
    sub029_NSvsmemory(i) = mean(sub029_NS)-mean(sub029_memory);
end

sub034_NSvscount = zeros(length(find(adj_pval_sub034(1:50)<0.05)),1);
sub034_count_sig = find(adj_pval_sub034(1:50)<0.05);
for i = 1:length(sub034_NSvscount)
    sub034_NS = [GC{30}(sub034_count_sig(i)+1,:) GC{31}(sub034_count_sig(i)+1,:)];
    sub034_count = [GC{26}(sub034_count_sig(i)+1,:) GC{28}(sub034_count_sig(i)+1,:)];
    sub034_NSvscount(i) = mean(sub034_NS)-mean(sub034_count);
end

sub034_NSvsmemory = zeros(length(find(adj_pval_sub034(51:100)<0.05)),1);
sub034_memory_sig = find(adj_pval_sub034(51:100)<0.05);
for i = 1:length(sub034_NSvsmemory)
    sub034_NS = [GC{30}(sub034_memory_sig(i)+1,:) GC{31}(sub034_memory_sig(i)+1,:)];
    sub034_memory = [GC{27}(sub034_memory_sig(i)+1,:) GC{29}(sub034_memory_sig(i)+1,:)];
    sub034_NSvsmemory(i) = mean(sub034_NS)-mean(sub034_memory);
end

%% Statistics - network

yeo_indices = {[(1:9)'; (51:58)'], [(10:15)'; (59:66)'], [(16:23)'; (67:73)'], [(24:30)'; (74:78)'], [(31:33)'; (79:80)'], [(34:37)'; (81:89)'], [(38:50)'; (90:100)']};
num_networks = length(yeo_indices);
LE_11Hz_Network = cell(length(d_combined),1);
for f = 1:length(d_combined)
    LE_11Hz_Network{f} = zeros(num_networks,size(LE_11Hz{f},2));
    for net = 1:7
        % Get indices of regions in this network
        idx = yeo_indices{net};
    
        % Take the mean across the 2nd dimension (regions)
        LE_11Hz_Network{f}(net,:) = mean(LE_11Hz{f}(idx,:));
    end
end

pval_sub001_count_yeo = zeros(num_networks,1);
pval_sub001_count_diff_yeo = zeros(num_networks,1);
pval_sub001_memory_yeo = zeros(num_networks,1);
pval_sub001_memory_diff_yeo = zeros(num_networks,1);
pval_sub003_count_yeo = zeros(num_networks,1);
pval_sub003_count_diff_yeo = zeros(num_networks,1);
pval_sub003_memory_yeo = zeros(num_networks,1);
pval_sub003_memory_diff_yeo = zeros(num_networks,1);
pval_sub010_count_yeo = zeros(num_networks,1);
pval_sub010_count_diff_yeo = zeros(num_networks,1);
pval_sub010_memory_yeo = zeros(num_networks,1);
pval_sub010_memory_diff_yeo = zeros(num_networks,1);
pval_sub029_count_yeo = zeros(num_networks,1);
pval_sub029_count_diff_yeo = zeros(num_networks,1);
pval_sub029_memory_yeo = zeros(num_networks,1);
pval_sub029_memory_diff_yeo = zeros(num_networks,1);
pval_sub034_count_yeo = zeros(num_networks,1);
pval_sub034_count_diff_yeo = zeros(num_networks,1);
pval_sub034_memory_yeo = zeros(num_networks,1);
pval_sub034_memory_diff_yeo = zeros(num_networks,1);

for net = 1:num_networks
    sub001_NS = LE_11Hz_Network{11}(net,:)';
    sub001_count = [LE_11Hz_Network{7}(net,:)'; LE_11Hz_Network{9}(net,:)'];
    pval_sub001_count_yeo(net) = permutation_test(sub001_NS, sub001_count, n_perm);
    pval_sub001_count_diff_yeo(net) = sign(mean(sub001_NS)-mean(sub001_count));
    sub001_memory = [LE_11Hz_Network{8}(net,:)'; LE_11Hz_Network{10}(net,:)'];
    pval_sub001_memory_yeo(net) = permutation_test(sub001_NS, sub001_memory, n_perm);
    pval_sub001_memory_diff_yeo(net) = sign(mean(sub001_NS)-mean(sub001_memory));

    sub003_NS = [LE_11Hz_Network{3}(net,:)'; LE_11Hz_Network{4}(net,:)'; LE_11Hz_Network{5}(net,:)'; LE_11Hz_Network{6}(net,:)'];
    sub003_count = LE_11Hz_Network{1}(net,:)';
    pval_sub003_count_yeo(net) = permutation_test(sub003_NS, sub003_count, n_perm);
    pval_sub003_count_diff_yeo(net) = sign(mean(sub003_NS)-mean(sub003_count));
    sub003_memory = LE_11Hz_Network{2}(net,:)';
    pval_sub003_memory_yeo(net) = permutation_test(sub003_NS, sub003_memory, n_perm);
    pval_sub003_memory_diff_yeo(net) = sign(mean(sub003_NS)-mean(sub003_memory));
    
    sub010_NS = [LE_11Hz_Network{17}(net,:)'; LE_11Hz_Network{18}(net,:)'; LE_11Hz_Network{19}(net,:)'];
    sub010_count = [LE_11Hz_Network{12}(net,:)'; LE_11Hz_Network{14}(net,:)'];
    pval_sub010_count_yeo(net) = permutation_test(sub010_NS, sub010_count, n_perm);
    pval_sub010_count_diff_yeo(net) = sign(mean(sub010_NS)-mean(sub010_count));
    sub010_memory = [LE_11Hz_Network{13}(net,:)'; LE_11Hz_Network{15}(net,:)'];
    pval_sub010_memory_yeo(net) = permutation_test(sub010_NS, sub010_memory, n_perm);
    pval_sub010_memory_diff_yeo(net) = sign(mean(sub010_NS)-mean(sub010_memory));

    sub029_NS = [LE_11Hz_Network{24}(net,:)'; LE_11Hz_Network{25}(net,:)'];
    sub029_count = [LE_11Hz_Network{20}(net,:)'; LE_11Hz_Network{22}(net,:)'];
    pval_sub029_count_yeo(net) = permutation_test(sub029_NS, sub029_count, n_perm);
    pval_sub029_count_diff_yeo(net) = sign(mean(sub029_NS)-mean(sub029_count));
    sub029_memory = [LE_11Hz_Network{21}(net,:)'; LE_11Hz_Network{23}(net,:)'];
    pval_sub029_memory_yeo(net) = permutation_test(sub029_NS, sub029_memory, n_perm);
    pval_sub029_memory_diff_yeo(net) = sign(mean(sub029_NS)-mean(sub029_memory));

    sub034_NS = [LE_11Hz_Network{30}(net,:)'; LE_11Hz_Network{31}(net,:)'];
    sub034_count = [LE_11Hz_Network{26}(net,:)'; LE_11Hz_Network{28}(net,:)'];
    pval_sub034_count_yeo(net) = permutation_test(sub034_NS, sub034_count, n_perm);
    pval_sub034_count_diff_yeo(net) = sign(mean(sub034_NS)-mean(sub034_count));
    sub034_memory = [LE_11Hz_Network{27}(net,:)'; LE_11Hz_Network{29}(net,:)'];
    pval_sub034_memory_yeo(net) = permutation_test(sub034_NS, sub034_memory, n_perm);    
    pval_sub034_memory_diff_yeo(net) = sign(mean(sub034_NS)-mean(sub034_memory));

end

[~,~,~,adj_pval_sub001_yeo] = fdr_bh_2015([pval_sub001_count_yeo; pval_sub001_memory_yeo]);
[~,~,~,adj_pval_sub003_yeo] = fdr_bh_2015([pval_sub003_count_yeo; pval_sub003_memory_yeo]);
[~,~,~,adj_pval_sub010_yeo] = fdr_bh_2015([pval_sub010_count_yeo; pval_sub010_memory_yeo]);
[~,~,~,adj_pval_sub029_yeo] = fdr_bh_2015([pval_sub029_count_yeo; pval_sub029_memory_yeo]);
[~,~,~,adj_pval_sub034_yeo] = fdr_bh_2015([pval_sub034_count_yeo; pval_sub034_memory_yeo]);

sum([pval_sub001_count_diff_yeo pval_sub001_memory_diff_yeo pval_sub003_count_diff_yeo pval_sub003_memory_diff_yeo pval_sub010_count_diff_yeo pval_sub010_memory_diff_yeo pval_sub029_count_diff_yeo pval_sub029_memory_diff_yeo pval_sub034_count_diff_yeo pval_sub034_memory_diff_yeo],2)




%% Plot contributions of each source - 11 Hz 

fr_ind = 12;
LE_11Hz_avg = zeros(num_reg,length(LE));
for f = 1:length(LE)
    LE_11Hz_avg(:,f) = mean(abs(squeeze(LE{f}(:,fr_ind,:)).^2),2);
end

addpath(genpath("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Scripts/ENIGMA-2.0.0"))

sub001_count_LE_11Hz_avg = log10(mean(LE_11Hz_avg(:,[7 9]),2));
sub001_memory_LE_11Hz_avg = log10(mean(LE_11Hz_avg(:,[8 10]),2));
sub001_NS_LE_11Hz_avg = log10(LE_11Hz_avg(:,11));

sub001_count_LE_11Hz_avg_diff = sub001_NS_LE_11Hz_avg - sub001_count_LE_11Hz_avg;
sub001_memory_LE_11Hz_avg_diff = sub001_NS_LE_11Hz_avg - sub001_memory_LE_11Hz_avg;

sub003_count_LE_11Hz_avg = log10(LE_11Hz_avg(:,1));
sub003_memory_LE_11Hz_avg = log10(LE_11Hz_avg(:,2));
sub003_NS_LE_11Hz_avg = log10(mean(LE_11Hz_avg(:,[3 4 5 6]),2));

sub003_count_LE_11Hz_avg_diff = sub003_NS_LE_11Hz_avg - sub003_count_LE_11Hz_avg;
sub003_memory_LE_11Hz_avg_diff = sub003_NS_LE_11Hz_avg - sub003_memory_LE_11Hz_avg;

sub010_count_LE_11Hz_avg = log10(mean(LE_11Hz_avg(:,[12 14]),2));
sub010_memory_LE_11Hz_avg = log10(mean(LE_11Hz_avg(:,[13 15]),2));
sub010_NS_LE_11Hz_avg = log10(mean(LE_11Hz_avg(:,[17 18 19]),2));

sub010_count_LE_11Hz_avg_diff = sub010_NS_LE_11Hz_avg - sub010_count_LE_11Hz_avg;
sub010_memory_LE_11Hz_avg_diff = sub010_NS_LE_11Hz_avg - sub010_memory_LE_11Hz_avg;

sub029_count_LE_11Hz_avg = log10(mean(LE_11Hz_avg(:,[20 22]),2));
sub029_memory_LE_11Hz_avg = log10(mean(LE_11Hz_avg(:,[21 23]),2));
sub029_NS_LE_11Hz_avg = log10(mean(LE_11Hz_avg(:,[24 25]),2));

sub029_count_LE_11Hz_avg_diff = sub029_NS_LE_11Hz_avg - sub029_count_LE_11Hz_avg;
sub029_memory_LE_11Hz_avg_diff = sub029_NS_LE_11Hz_avg - sub029_memory_LE_11Hz_avg;

sub034_count_LE_11Hz_avg = log10(mean(LE_11Hz_avg(:,[26 28]),2));
sub034_memory_LE_11Hz_avg = log10(mean(LE_11Hz_avg(:,[27 29]),2));
sub034_NS_LE_11Hz_avg = log10(mean(LE_11Hz_avg(:,[30 31]),2));

sub034_count_LE_11Hz_avg_diff = sub034_NS_LE_11Hz_avg - sub034_count_LE_11Hz_avg;
sub034_memory_LE_11Hz_avg_diff = sub034_NS_LE_11Hz_avg - sub034_memory_LE_11Hz_avg;

min_val = min([sub001_count_LE_11Hz_avg_diff; sub001_memory_LE_11Hz_avg_diff; sub003_count_LE_11Hz_avg_diff; sub003_memory_LE_11Hz_avg_diff; sub010_count_LE_11Hz_avg_diff; sub010_memory_LE_11Hz_avg_diff; sub029_count_LE_11Hz_avg_diff; sub029_memory_LE_11Hz_avg_diff; sub034_count_LE_11Hz_avg_diff; sub034_memory_LE_11Hz_avg_diff]);
max_val = max([sub001_count_LE_11Hz_avg_diff; sub001_memory_LE_11Hz_avg_diff; sub003_count_LE_11Hz_avg_diff; sub003_memory_LE_11Hz_avg_diff; sub010_count_LE_11Hz_avg_diff; sub010_memory_LE_11Hz_avg_diff; sub029_count_LE_11Hz_avg_diff; sub029_memory_LE_11Hz_avg_diff; sub034_count_LE_11Hz_avg_diff; sub034_memory_LE_11Hz_avg_diff]);

%% Plot contributions of each source - 11 Hz - sub001 

sub001_count_surf = parcel_to_surface(sub001_count_LE_11Hz_avg_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub001_count_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub001_memory_surf = parcel_to_surface(sub001_memory_LE_11Hz_avg_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub001_memory_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')


%% Plot contributions of each source - 11 Hz - sub003

sub003_count_surf = parcel_to_surface(sub003_count_LE_11Hz_avg_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub003_count_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub003_memory_surf = parcel_to_surface(sub003_memory_LE_11Hz_avg_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub003_memory_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

%% Plot contributions of each source - 11 Hz - sub010

sub010_count_surf = parcel_to_surface(sub010_count_LE_11Hz_avg_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub010_count_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub010_memory_surf = parcel_to_surface(sub010_memory_LE_11Hz_avg_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub010_memory_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

%% Plot contributions of each source - 11 Hz - sub029

sub029_count_surf = parcel_to_surface(sub029_count_LE_11Hz_avg_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub029_count_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub029_memory_surf = parcel_to_surface(sub029_memory_LE_11Hz_avg_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub029_memory_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

%% Plot contributions of each source - 11 Hz - sub034

sub034_count_surf = parcel_to_surface(sub034_count_LE_11Hz_avg_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub034_count_surf, 'surface_name', 'fsa5', 'color_range', ...
              [min_val -min_val], 'cmap', 'RdBu_r')

sub034_memory_surf = parcel_to_surface(sub034_memory_LE_11Hz_avg_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub034_memory_surf, 'surface_name', 'fsa5', 'color_range', ...
              [min_val -min_val], 'cmap', 'RdBu_r')

%% Statistics - contributions of each source - 11 Hz

fr_ind = 12;
LE_11Hz = cell(31,1);
for f = 1:length(LE_11Hz)
    LE_11Hz{f} = abs(squeeze(LE{f}(:,fr_ind,:))).^2;
end

n_perm = 5000;
num_reg = 100;

pval_sub001_count_LE = zeros(num_reg,1);
pval_sub001_memory_LE = zeros(num_reg,1);
pval_sub003_count_LE = zeros(num_reg,1);
pval_sub003_memory_LE = zeros(num_reg,1);
pval_sub010_count_LE = zeros(num_reg,1);
pval_sub010_memory_LE = zeros(num_reg,1);
pval_sub029_count_LE = zeros(num_reg,1);
pval_sub029_memory_LE = zeros(num_reg,1);
pval_sub034_count_LE = zeros(num_reg,1);
pval_sub034_memory_LE = zeros(num_reg,1);

for reg = 1:num_reg

    sub001_NS = LE_11Hz{11}(reg,:)';
    sub001_count = [LE_11Hz{7}(reg,:)'; LE_11Hz{9}(reg,:)'];
    pval_sub001_count_LE(reg) = permutation_test(sub001_NS, sub001_count, n_perm);
    sub001_memory = [LE_11Hz{8}(reg,:)'; LE_11Hz{10}(reg,:)'];
    pval_sub001_memory_LE(reg) = permutation_test(sub001_NS, sub001_memory, n_perm);

    sub003_NS = [LE_11Hz{3}(reg,:)'; LE_11Hz{4}(reg,:)'; LE_11Hz{5}(reg,:)'; LE_11Hz{6}(reg,:)'];
    sub003_count = LE_11Hz{1}(reg,:)';
    pval_sub003_count_LE(reg) = permutation_test(sub003_NS, sub003_count, n_perm);
    sub003_memory = LE_11Hz{2}(reg,:)';
    pval_sub003_memory_LE(reg) = permutation_test(sub003_NS, sub003_memory, n_perm);
    
    sub010_NS = [LE_11Hz{17}(reg,:)'; LE_11Hz{18}(reg,:)'; LE_11Hz{19}(reg,:)'];
    sub010_count = [LE_11Hz{12}(reg,:)'; LE_11Hz{14}(reg,:)'];
    pval_sub010_count_LE(reg) = permutation_test(sub010_NS, sub010_count, n_perm);
    sub010_memory = [LE_11Hz{13}(reg,:)'; LE_11Hz{15}(reg,:)'];
    pval_sub010_memory_LE(reg) = permutation_test(sub010_NS, sub010_memory, n_perm);

    sub029_NS = [LE_11Hz{24}(reg,:)'; LE_11Hz{25}(reg,:)'];
    sub029_count = [LE_11Hz{20}(reg,:)'; LE_11Hz{22}(reg,:)'];
    pval_sub029_count_LE(reg) = permutation_test(sub029_NS, sub029_count, n_perm);
    sub029_memory = [LE_11Hz{21}(reg,:)'; LE_11Hz{23}(reg,:)'];
    pval_sub029_memory_LE(reg) = permutation_test(sub029_NS, sub029_memory, n_perm);

    sub034_NS = [LE_11Hz{30}(reg,:)'; LE_11Hz{31}(reg,:)'];
    sub034_count = [LE_11Hz{26}(reg,:)'; LE_11Hz{28}(reg,:)'];
    pval_sub034_count_LE(reg) = permutation_test(sub034_NS, sub034_count, n_perm);
    sub034_memory = [LE_11Hz{27}(reg,:)'; LE_11Hz{29}(reg,:)'];
    pval_sub034_memory_LE(reg) = permutation_test(sub034_NS, sub034_memory, n_perm);
end

[~,~,~,adj_pval_sub001_LE] = fdr_bh_2015([pval_sub001_count_LE; pval_sub001_memory_LE]);
[~,~,~,adj_pval_sub003_LE] = fdr_bh_2015([pval_sub003_count_LE; pval_sub003_memory_LE]);
[~,~,~,adj_pval_sub010_LE] = fdr_bh_2015([pval_sub010_count_LE; pval_sub010_memory_LE]);
[~,~,~,adj_pval_sub029_LE] = fdr_bh_2015([pval_sub029_count_LE; pval_sub029_memory_LE]);
[~,~,~,adj_pval_sub034_LE] = fdr_bh_2015([pval_sub034_count_LE; pval_sub034_memory_LE]);

% Significantly different regions and the direction of change (increase or
% decrease), as well as direction of change of regions per Yeo network

sub001_count_LE_sig_regs = find(adj_pval_sub001_LE(1:100)<0.05);
sub001_count_LE_sig_regs_direction = zeros(length(sub001_count_LE_sig_regs),1);
for reg = 1:length(sub001_count_LE_sig_regs)
    sub001_NS = LE_11Hz{11}(sub001_count_LE_sig_regs(reg),:)';
    sub001_count = [LE_11Hz{7}(sub001_count_LE_sig_regs(reg),:)'; LE_11Hz{9}(sub001_count_LE_sig_regs(reg),:)'];
    sub001_count_LE_sig_regs_direction(reg) = sign(mean(sub001_count)-mean(sub001_NS));
end
sub001_count_LE_yeo = zeros(7,1);
sub001_count_LE_yeo_increase = zeros(7,1);
sub001_count_LE_yeo_decrease = zeros(7,1);
for reg = 1:length(sub001_count_LE_sig_regs)
    for y = 1:7
        if ismember(sub001_count_LE_sig_regs(reg), yeo_indices{y})
            sub001_count_LE_yeo(y) = sub001_count_LE_yeo(y) + 1;
            if sub001_count_LE_sig_regs_direction(reg) == 1
                sub001_count_LE_yeo_decrease(y) = sub001_count_LE_yeo_decrease(y) + 1;
            else
                sub001_count_LE_yeo_increase(y) = sub001_count_LE_yeo_increase(y) + 1;
            end
        end
    end
end

sub001_memory_LE_sig_regs = find(adj_pval_sub001_LE(101:200)<0.05);
sub001_memory_LE_sig_regs_direction = zeros(length(sub001_memory_LE_sig_regs),1);
for reg = 1:length(sub001_memory_LE_sig_regs)
    sub001_NS = LE_11Hz{11}(sub001_memory_LE_sig_regs(reg),:)';
    sub001_memory = [LE_11Hz{8}(sub001_memory_LE_sig_regs(reg),:)'; LE_11Hz{10}(sub001_memory_LE_sig_regs(reg),:)'];
    sub001_memory_LE_sig_regs_direction(reg) = sign(mean(sub001_memory)-mean(sub001_NS));
end
sub001_memory_LE_yeo = zeros(7,1);
sub001_memory_LE_yeo_increase = zeros(7,1);
sub001_memory_LE_yeo_decrease = zeros(7,1);
for reg = 1:length(sub001_memory_LE_sig_regs)
    for y = 1:7
        if ismember(sub001_memory_LE_sig_regs(reg), yeo_indices{y})
            sub001_memory_LE_yeo(y) = sub001_memory_LE_yeo(y) + 1;
            if sub001_memory_LE_sig_regs_direction(reg) == 1
                sub001_memory_LE_yeo_decrease(y) = sub001_memory_LE_yeo_decrease(y) + 1;
            else
                sub001_memory_LE_yeo_increase(y) = sub001_memory_LE_yeo_increase(y) + 1;
            end
        end
    end
end

sub003_count_LE_sig_regs = find(adj_pval_sub003_LE(1:100)<0.05);
sub003_count_LE_sig_regs_direction = zeros(length(sub003_count_LE_sig_regs),1);
sub003_count_LE_sig_regs_diff = zeros(length(sub003_count_LE_sig_regs),1);
for reg = 1:length(sub003_count_LE_sig_regs)
    sub003_NS = [LE_11Hz{3}(sub003_count_LE_sig_regs(reg),:)'; LE_11Hz{4}(sub003_count_LE_sig_regs(reg),:)'; LE_11Hz{5}(sub003_count_LE_sig_regs(reg),:)'; LE_11Hz{6}(sub003_count_LE_sig_regs(reg),:)'];
    sub003_count = LE_11Hz{1}(sub003_count_LE_sig_regs(reg),:)';
    sub003_count_LE_sig_regs_direction(reg) = sign(mean(sub003_count)-mean(sub003_NS));
    sub003_count_LE_sig_regs_diff(reg) = mean(sub003_count)-mean(sub003_NS);
end
sub003_count_LE_yeo = zeros(7,1);
sub003_count_LE_yeo_increase = zeros(7,1);
sub003_count_LE_yeo_decrease = zeros(7,1);
for reg = 1:length(sub003_count_LE_sig_regs)
    for y = 1:7
        if ismember(sub003_count_LE_sig_regs(reg), yeo_indices{y})
            sub003_count_LE_yeo(y) = sub003_count_LE_yeo(y) + 1;
            if sub003_count_LE_sig_regs_direction(reg) == 1
                sub003_count_LE_yeo_decrease(y) = sub003_count_LE_yeo_decrease(y) + 1;
            else
                sub003_count_LE_yeo_increase(y) = sub003_count_LE_yeo_increase(y) + 1;
            end
        end
    end
end

sub003_memory_LE_sig_regs = find(adj_pval_sub003_LE(101:200)<0.05);
sub003_memory_LE_sig_regs_direction = zeros(length(sub003_memory_LE_sig_regs),1);
sub003_memory_LE_sig_regs_diff = zeros(length(sub003_memory_LE_sig_regs),1);
for reg = 1:length(sub003_memory_LE_sig_regs)
    sub003_NS = [LE_11Hz{3}(sub003_memory_LE_sig_regs(reg),:)'; LE_11Hz{4}(sub003_memory_LE_sig_regs(reg),:)'; LE_11Hz{5}(sub003_memory_LE_sig_regs(reg),:)'; LE_11Hz{6}(sub003_memory_LE_sig_regs(reg),:)'];
    sub003_memory = LE_11Hz{2}(sub003_memory_LE_sig_regs(reg),:)';
    sub003_memory_LE_sig_regs_direction(reg) = sign(mean(sub003_memory)-mean(sub003_NS));
    sub003_memory_LE_sig_regs_diff(reg) = mean(sub003_memory)-mean(sub003_NS);
end
sub003_memory_LE_yeo = zeros(7,1);
sub003_memory_LE_yeo_increase = zeros(7,1);
sub003_memory_LE_yeo_decrease = zeros(7,1);
for reg = 1:length(sub003_memory_LE_sig_regs)
    for y = 1:7
        if ismember(sub003_memory_LE_sig_regs(reg), yeo_indices{y})
            sub003_memory_LE_yeo(y) = sub003_memory_LE_yeo(y) + 1;
            if sub003_memory_LE_sig_regs_direction(reg) == 1
                sub003_memory_LE_yeo_decrease(y) = sub003_memory_LE_yeo_decrease(y) + 1;
            else
                sub003_memory_LE_yeo_increase(y) = sub003_memory_LE_yeo_increase(y) + 1;
            end
        end
    end
end

sub010_count_LE_sig_regs = find(adj_pval_sub010_LE(1:100)<0.05);
sub010_count_LE_sig_regs_direction = zeros(length(sub010_count_LE_sig_regs),1);
for reg = 1:length(sub010_count_LE_sig_regs)
    sub010_NS = [LE_11Hz{17}(sub010_count_LE_sig_regs(reg),:)'; LE_11Hz{18}(sub010_count_LE_sig_regs(reg),:)'; LE_11Hz{19}(sub010_count_LE_sig_regs(reg),:)'];
    sub010_count = [LE_11Hz{12}(sub010_count_LE_sig_regs(reg),:)'; LE_11Hz{14}(sub010_count_LE_sig_regs(reg),:)'];
    sub010_count_LE_sig_regs_direction(reg) = sign(mean(sub010_count)-mean(sub010_NS));
end
sub010_count_LE_yeo = zeros(7,1);
sub010_count_LE_yeo_increase = zeros(7,1);
sub010_count_LE_yeo_decrease = zeros(7,1);
for reg = 1:length(sub010_count_LE_sig_regs)
    for y = 1:7
        if ismember(sub010_count_LE_sig_regs(reg), yeo_indices{y})
            sub010_count_LE_yeo(y) = sub010_count_LE_yeo(y) + 1;
            if sub010_count_LE_sig_regs_direction(reg) == 1
                sub010_count_LE_yeo_decrease(y) = sub010_count_LE_yeo_decrease(y) + 1;
            else
                sub010_count_LE_yeo_increase(y) = sub010_count_LE_yeo_increase(y) + 1;
            end
        end
    end
end

sub010_memory_LE_sig_regs = find(adj_pval_sub010_LE(101:200)<0.05);
sub010_memory_LE_sig_regs_direction = zeros(length(sub010_memory_LE_sig_regs),1);
for reg = 1:length(sub010_memory_LE_sig_regs)
    sub010_NS = [LE_11Hz{17}(sub010_memory_LE_sig_regs(reg),:)'; LE_11Hz{18}(sub010_memory_LE_sig_regs(reg),:)'; LE_11Hz{19}(sub010_memory_LE_sig_regs(reg),:)'];
    sub010_memory = [LE_11Hz{13}(sub010_memory_LE_sig_regs(reg),:)'; LE_11Hz{15}(sub010_memory_LE_sig_regs(reg),:)'];
    sub010_memory_LE_sig_regs_direction(reg) = sign(mean(sub010_memory)-mean(sub010_NS));
end
sub010_memory_LE_yeo = zeros(7,1);
sub010_memory_LE_yeo_increase = zeros(7,1);
sub010_memory_LE_yeo_decrease = zeros(7,1);
for reg = 1:length(sub010_memory_LE_sig_regs)
    for y = 1:7
        if ismember(sub010_memory_LE_sig_regs(reg), yeo_indices{y})
            sub010_memory_LE_yeo(y) = sub010_memory_LE_yeo(y) + 1;
            if sub010_memory_LE_sig_regs_direction(reg) == 1
                sub010_memory_LE_yeo_decrease(y) = sub010_memory_LE_yeo_decrease(y) + 1;
            else
                sub010_memory_LE_yeo_increase(y) = sub010_memory_LE_yeo_increase(y) + 1;
            end
        end
    end
end

sub029_count_LE_sig_regs = find(adj_pval_sub029_LE(1:100)<0.05);
sub029_count_LE_sig_regs_direction = zeros(length(sub029_count_LE_sig_regs),1);
for reg = 1:length(sub029_count_LE_sig_regs)
    sub029_NS = [LE_11Hz{24}(sub029_count_LE_sig_regs(reg),:)'; LE_11Hz{25}(sub029_count_LE_sig_regs(reg),:)'];
    sub029_count = [LE_11Hz{20}(sub029_count_LE_sig_regs(reg),:)'; LE_11Hz{22}(sub029_count_LE_sig_regs(reg),:)'];
    sub029_count_LE_sig_regs_direction(reg) = sign(mean(sub029_count)-mean(sub029_NS));
end
sub029_count_LE_yeo = zeros(7,1);
sub029_count_LE_yeo_increase = zeros(7,1);
sub029_count_LE_yeo_decrease = zeros(7,1);
for reg = 1:length(sub029_count_LE_sig_regs)
    for y = 1:7
        if ismember(sub029_count_LE_sig_regs(reg), yeo_indices{y})
            sub029_count_LE_yeo(y) = sub029_count_LE_yeo(y) + 1;
            if sub029_count_LE_sig_regs_direction(reg) == 1
                sub029_count_LE_yeo_decrease(y) = sub029_count_LE_yeo_decrease(y) + 1;
            else
                sub029_count_LE_yeo_increase(y) = sub029_count_LE_yeo_increase(y) + 1;
            end
        end
    end
end

sub029_memory_LE_sig_regs = find(adj_pval_sub029_LE(101:200)<0.05);
sub029_memory_LE_sig_regs_direction = zeros(length(sub029_memory_LE_sig_regs),1);
for reg = 1:length(sub029_memory_LE_sig_regs)
    sub029_NS = [LE_11Hz{24}(sub029_memory_LE_sig_regs(reg),:)'; LE_11Hz{25}(sub029_memory_LE_sig_regs(reg),:)'];
    sub029_memory = [LE_11Hz{21}(sub029_memory_LE_sig_regs(reg),:)'; LE_11Hz{23}(sub029_memory_LE_sig_regs(reg),:)'];
    sub029_memory_LE_sig_regs_direction(reg) = sign(mean(sub029_memory)-mean(sub029_NS));
end
sub029_memory_LE_yeo = zeros(7,1);
sub029_memory_LE_yeo_increase = zeros(7,1);
sub029_memory_LE_yeo_decrease = zeros(7,1);
for reg = 1:length(sub029_memory_LE_sig_regs)
    for y = 1:7
        if ismember(sub029_memory_LE_sig_regs(reg), yeo_indices{y})
            sub029_memory_LE_yeo(y) = sub029_memory_LE_yeo(y) + 1;
            if sub029_memory_LE_sig_regs_direction(reg) == 1
                sub029_memory_LE_yeo_decrease(y) = sub029_memory_LE_yeo_decrease(y) + 1;
            else
                sub029_memory_LE_yeo_increase(y) = sub029_memory_LE_yeo_increase(y) + 1;
            end
        end
    end
end

sub034_count_LE_sig_regs = find(adj_pval_sub034_LE(1:100)<0.05);
sub034_count_LE_sig_regs_direction = zeros(length(sub034_count_LE_sig_regs),1);
for reg = 1:length(sub034_count_LE_sig_regs)
    sub034_NS = [LE_11Hz{30}(sub034_count_LE_sig_regs(reg),:)'; LE_11Hz{31}(sub034_count_LE_sig_regs(reg),:)'];
    sub034_count = [LE_11Hz{26}(sub034_count_LE_sig_regs(reg),:)'; LE_11Hz{28}(sub034_count_LE_sig_regs(reg),:)'];
    sub034_count_LE_sig_regs_direction(reg) = sign(mean(sub034_count)-mean(sub034_NS));
end
sub034_count_LE_yeo = zeros(7,1);
sub034_count_LE_yeo_increase = zeros(7,1);
sub034_count_LE_yeo_decrease = zeros(7,1);
for reg = 1:length(sub034_count_LE_sig_regs)
    for y = 1:7
        if ismember(sub034_count_LE_sig_regs(reg), yeo_indices{y})
            sub034_count_LE_yeo(y) = sub034_count_LE_yeo(y) + 1;
            if sub034_count_LE_sig_regs_direction(reg) == 1
                sub034_count_LE_yeo_decrease(y) = sub034_count_LE_yeo_decrease(y) + 1;
            else
                sub034_count_LE_yeo_increase(y) = sub034_count_LE_yeo_increase(y) + 1;
            end
        end
    end
end

sub034_memory_LE_sig_regs = find(adj_pval_sub034_LE(101:200)<0.05);
sub034_memory_LE_sig_regs_direction = zeros(length(sub034_memory_LE_sig_regs),1);
for reg = 1:length(sub034_memory_LE_sig_regs)
    sub034_NS = [LE_11Hz{30}(sub034_memory_LE_sig_regs(reg),:)'; LE_11Hz{31}(sub034_memory_LE_sig_regs(reg),:)'];
    sub034_memory = [LE_11Hz{27}(sub034_memory_LE_sig_regs(reg),:)'; LE_11Hz{29}(sub034_memory_LE_sig_regs(reg),:)'];
    sub034_memory_LE_sig_regs_direction(reg) = sign(mean(sub034_memory)-mean(sub034_NS));
end
sub034_memory_LE_yeo = zeros(7,1);
sub034_memory_LE_yeo_increase = zeros(7,1);
sub034_memory_LE_yeo_decrease = zeros(7,1);
for reg = 1:length(sub034_memory_LE_sig_regs)
    for y = 1:7
        if ismember(sub034_memory_LE_sig_regs(reg), yeo_indices{y})
            sub034_memory_LE_yeo(y) = sub034_memory_LE_yeo(y) + 1;
            if sub034_memory_LE_sig_regs_direction(reg) == 1
                sub034_memory_LE_yeo_decrease(y) = sub034_memory_LE_yeo_decrease(y) + 1;
            else
                sub034_memory_LE_yeo_increase(y) = sub034_memory_LE_yeo_increase(y) + 1;
            end
        end
    end
end

all_LE_sig_regs = {sub001_count_LE_sig_regs, sub001_memory_LE_sig_regs, sub003_count_LE_sig_regs, sub003_memory_LE_sig_regs, sub010_count_LE_sig_regs, sub010_memory_LE_sig_regs, sub029_count_LE_sig_regs, sub029_memory_LE_sig_regs, sub034_count_LE_sig_regs, sub034_memory_LE_sig_regs};
all_LE_sig_regs_direction = {sub001_count_LE_sig_regs_direction, sub001_memory_LE_sig_regs_direction, sub003_count_LE_sig_regs_direction, sub003_memory_LE_sig_regs_direction, sub010_count_LE_sig_regs_direction, sub010_memory_LE_sig_regs_direction, sub029_count_LE_sig_regs_direction, sub029_memory_LE_sig_regs_direction, sub034_count_LE_sig_regs_direction, sub034_memory_LE_sig_regs_direction};
common_LE_sig_regs = all_LE_sig_regs{1};
for i = 2:length(all_LE_sig_regs)
    common_LE_sig_regs = intersect(common_LE_sig_regs, all_LE_sig_regs{i});
end
common_LE_sig_regs_direction = zeros(length(common_LE_sig_regs),length(all_LE_sig_regs));
for i = 1:length(common_LE_sig_regs)
    for j = 1:length(all_LE_sig_regs)
        reg_idx = find(all_LE_sig_regs{j}==common_LE_sig_regs(i));
        common_LE_sig_regs_direction(i,j) = all_LE_sig_regs_direction{j}(reg_idx);
    end
end

all_LE_values = vertcat(all_LE_sig_regs{:});
[unique_vals, ~, idx] = unique(all_LE_values);
counts = accumarray(idx, 1);
[counts_sorted, sort_idx] = sort(counts, 'descend');
most_common_LE_vals = unique_vals(sort_idx);

most_common_LE_vals_direction = zeros(length(most_common_LE_vals),length(all_LE_sig_regs));
for i = 1:length(most_common_LE_vals)
    for j = 1:length(all_LE_sig_regs)
        reg_idx = find(all_LE_sig_regs{j}==most_common_LE_vals(i));
        if isempty(reg_idx)
            most_common_LE_vals_direction(i,j) = 0;
        else
            most_common_LE_vals_direction(i,j) = all_LE_sig_regs_direction{j}(reg_idx);
        end
    end
end


