%% Compute PSD

cd("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Results/SourceRecon/Combined/Truncated/CommonFilter/New")
d = dir('*.mat');
fs = 200; % sampling rate 
windowLength = 200;
nfft = windowLength;
overlap = 100;
num_reg = 100;
num_freqs = floor(nfft/2)+1;

power = cell(length(d),1);
powerPerSource = cell(length(d),1);

for f = 1:length(d)

    % Load timeseries 
    vs = importdata(d(f).name);
    
    % Reshape vs into a 3D matrix
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

    power{f} = zeros(numFullSegments, num_freqs);
    powerPerSource{f} = zeros(numFullSegments, num_reg, num_freqs);

    for ep = 1:numFullSegments
        vs_2D = squeeze(vs_3D(:,ep,:));
        for reg = 1:num_reg
            [Pxx, freqs] = pwelch(vs_2D(reg,:), windowLength, overlap, nfft, fs);
            powerPerSource{f}(ep,reg,:) = Pxx;  % Store the power spectrum for this source
        end
        power{f}(ep,:) = mean(squeeze(powerPerSource{f}(ep,:,:)));
    end
    f
end

%% Save each cell individually

cd("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Results/PSD/SourceSpace/ForPython")
for f = 1:length(power)
    pow = power{f};
    save(['power_' num2str(f) '.mat'], 'pow')
end

%% Plot - EEG - sub001

sub001_counting_avg = mean([power{7}; power{9}])';
sub001_counting_std =  std([power{7}; power{9}],0,1);
sub001_counting_SE = sub001_counting_std/sqrt(size([power{7}; power{9}],1));
sub001_memory_avg = mean([power{8}; power{10}])';
sub001_memory_std = std([power{8}; power{10}],0,1);
sub001_memory_SE = sub001_memory_std/sqrt(size([power{8}; power{10}],1));
sub001_NS_avg = mean(power{11})';
sub001_NS_std = std(power{11},0,1);
sub001_NS_SE = sub001_NS_std/sqrt(size(power{11},1));

% Plot the result with standard error shading
figure;
hold on;

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub001_counting_avg(1:50)' + sub001_counting_SE(1:50), fliplr(sub001_counting_avg(1:50)' - sub001_counting_SE(1:50))], 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p1 = plot((0:49)', sub001_counting_avg(1:50), 'LineWidth', 2, 'Color', 'r');
hold on

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub001_memory_avg(1:50)' + sub001_memory_SE(1:50), fliplr(sub001_memory_avg(1:50)' - sub001_memory_SE(1:50))], [0.85 0.325 0.098], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p2 = plot((0:49)', sub001_memory_avg(1:50), 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);
hold on

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub001_NS_avg(1:50)' + sub001_NS_SE(1:50), fliplr(sub001_NS_avg(1:50)' - sub001_NS_SE(1:50))], 'b', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p3 = plot((0:49)', sub001_NS_avg(1:50), 'LineWidth', 2, 'Color', 'b');

legend([p1,p2,p3],{'Counting','Memory','NS'},'Location','NorthEast')
xlabel('Frequency (Hz)');
ylabel('Power');
title('sub001 EEG Power Spectrum (with SE)');
grid on;
hold off;

%% Plot - combined EEG+MEG - sub003

sub003_counting_avg = mean(power{1})';
sub003_counting_std = std(power{1},0,1);
sub003_counting_SE = sub003_counting_std/sqrt(size(power{1},1));
sub003_memory_avg = mean(power{2})';
sub003_memory_std = std(power{2},0,1);
sub003_memory_SE = sub003_memory_std/sqrt(size(power{2},1));
sub003_NS_avg = mean([power{3}; power{4}; power{5}; power{6}])';
sub003_NS_std = std([power{3}; power{4}; power{5}; power{6}],0,1);
sub003_NS_SE = sub003_NS_std/sqrt(size([power{3}; power{4}; power{5}; power{6}],1));

% Plot the result with standard error shading
figure;
hold on;

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub003_counting_avg(1:50)' + sub003_counting_SE(1:50), fliplr(sub003_counting_avg(1:50)' - sub003_counting_SE(1:50))], 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p1 = plot((0:49)', sub003_counting_avg(1:50), 'LineWidth', 2, 'Color', 'r');
hold on

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub003_memory_avg(1:50)' + sub003_memory_SE(1:50), fliplr(sub003_memory_avg(1:50)' - sub003_memory_SE(1:50))], [0.85 0.325 0.098], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p2 = plot((0:49)', sub003_memory_avg(1:50), 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);
hold on

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub003_NS_avg(1:50)' + sub003_NS_SE(1:50), fliplr(sub003_NS_avg(1:50)' - sub003_NS_SE(1:50))], 'b', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p3 = plot((0:49)', sub003_NS_avg(1:50), 'LineWidth', 2, 'Color', 'b');

legend([p1,p2,p3],{'Counting','Memory','NS'},'Location','NorthEast')
xlabel('Frequency (Hz)');
ylabel('Power');
title('sub003 Power Spectrum (with SE)');
grid on;
hold off;

%% Plot - combined EEG+MEG - sub010

sub010_counting_avg = mean([power{12}; power{14}])';
sub010_counting_std =  std([power{12}; power{14}],0,1);
sub010_counting_SE = sub010_counting_std/sqrt(size([power{12}; power{14}],1));
sub010_memory_avg = mean([power{13}; power{15}])';
sub010_memory_std = std([power{13}; power{15}],0,1);
sub010_memory_SE = sub010_memory_std/sqrt(size([power{13}; power{15}],1));
sub010_NS_avg = mean([power{17}; power{18}; power{19}])';
sub010_NS_std = std([power{17}; power{18}; power{19}],0,1);
sub010_NS_SE = sub010_NS_std/sqrt(size([power{17}; power{18}; power{19}],1));

% Plot the result with standard error shading
figure;
hold on;

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub010_counting_avg(1:50)' + sub010_counting_SE(1:50), fliplr(sub010_counting_avg(1:50)' - sub010_counting_SE(1:50))], 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p1 = plot((0:49)', sub010_counting_avg(1:50), 'LineWidth', 2, 'Color', 'r');
hold on

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub010_memory_avg(1:50)' + sub010_memory_SE(1:50), fliplr(sub010_memory_avg(1:50)' - sub010_memory_SE(1:50))], [0.85 0.325 0.098], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p2 = plot((0:49)', sub010_memory_avg(1:50), 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);
hold on

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub010_NS_avg(1:50)' + sub010_NS_SE(1:50), fliplr(sub010_NS_avg(1:50)' - sub010_NS_SE(1:50))], 'b', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p3 = plot((0:49)', sub010_NS_avg(1:50), 'LineWidth', 2, 'Color', 'b');

legend([p1,p2,p3],{'Counting','Memory','NS'},'Location','NorthEast')
xlabel('Frequency (Hz)');
ylabel('Power');
title('sub010 Power Spectrum (with SE)');
grid on;
hold off;

%% Plot - MEG - sub029

sub029_counting_avg = mean([power{20}; power{22}])';
sub029_counting_std =  std([power{20}; power{22}],0,1);
sub029_counting_SE = sub029_counting_std/sqrt(size([power{20}; power{22}],1));
sub029_memory_avg = mean([power{21}; power{23}])';
sub029_memory_std = std([power{21}; power{23}],0,1);
sub029_memory_SE = sub029_memory_std/sqrt(size([power{21}; power{23}],1));
sub029_NS_avg = mean([power{24}; power{25}])';
sub029_NS_std = std([power{24}; power{25}],0,1);
sub029_NS_SE = sub029_NS_std/sqrt(size([power{24}; power{25}],1));

% Plot the result with standard error shading
figure;
hold on;

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub029_counting_avg(1:50)' + sub029_counting_SE(1:50), fliplr(sub029_counting_avg(1:50)' - sub029_counting_SE(1:50))], 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p1 = plot((0:49)', sub029_counting_avg(1:50), 'LineWidth', 2, 'Color', 'r');
hold on

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub029_memory_avg(1:50)' + sub029_memory_SE(1:50), fliplr(sub029_memory_avg(1:50)' - sub029_memory_SE(1:50))], [0.85 0.325 0.098], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p2 = plot((0:49)', sub029_memory_avg(1:50), 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);
hold on

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub029_NS_avg(1:50)' + sub029_NS_SE(1:50), fliplr(sub029_NS_avg(1:50)' - sub029_NS_SE(1:50))], 'b', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p3 = plot((0:49)', sub029_NS_avg(1:50), 'LineWidth', 2, 'Color', 'b');

legend([p1,p2,p3],{'Counting','Memory','NS'},'Location','NorthEast')
xlabel('Frequency (Hz)');
ylabel('Power');
title('sub029 Power Spectrum (with SE)');
grid on;
hold off;

%% Plot - EEG+MEG - sub034

sub034_counting_avg = mean([power{26}; power{28}])';
sub034_counting_std =  std([power{26}; power{28}],0,1);
sub034_counting_SE = sub034_counting_std/sqrt(size([power{26}; power{28}],1));
sub034_memory_avg = mean([power{27}; power{29}])';
sub034_memory_std = std([power{27}; power{29}],0,1);
sub034_memory_SE = sub034_memory_std/sqrt(size([power{27}; power{29}],1));
sub034_NS_avg = mean([power{30}; power{31}])';
sub034_NS_std = std([power{30}; power{31}],0,1);
sub034_NS_SE = sub034_NS_std/sqrt(size([power{30}; power{31}],1));

% Plot the result with standard error shading
figure;
hold on;

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub034_counting_avg(1:50)' + sub034_counting_SE(1:50), fliplr(sub034_counting_avg(1:50)' - sub034_counting_SE(1:50))], 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p1 = plot((0:49)', sub034_counting_avg(1:50), 'LineWidth', 2, 'Color', 'r');
hold on

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub034_memory_avg(1:50)' + sub034_memory_SE(1:50), fliplr(sub034_memory_avg(1:50)' - sub034_memory_SE(1:50))], [0.85 0.325 0.098], ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p2 = plot((0:49)', sub034_memory_avg(1:50), 'LineWidth', 2, 'Color', [0.85 0.325 0.098]);
hold on

% Shaded area (mean ± SE)
fill([(0:49), fliplr((0:49))], [sub034_NS_avg(1:50)' + sub034_NS_SE(1:50), fliplr(sub034_NS_avg(1:50)' - sub034_NS_SE(1:50))], 'b', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.2);
hold on 
% Plot the mean power spectrum
p3 = plot((0:49)', sub034_NS_avg(1:50), 'LineWidth', 2, 'Color', 'b');

legend([p1,p2,p3],{'Counting','Memory','NS'},'Location','NorthEast')
xlabel('Frequency (Hz)');
ylabel('Power');
title('sub034 Power Spectrum (with SE)');
grid on;
hold off;

%% Topography of alpha power - values

sub001_counting_alpha = log10(mean([squeeze(mean(powerPerSource{7}(:,:,9:14),3)); squeeze(mean(powerPerSource{9}(:,:,9:14),3))]));
sub001_memory_alpha = log10(mean([squeeze(mean(powerPerSource{8}(:,:,9:14),3)); squeeze(mean(powerPerSource{10}(:,:,9:14),3))]));
sub001_NS_alpha = log10(mean(squeeze(mean(powerPerSource{11}(:,:,9:14),3))));

sub001_counting_alpha_diff = sub001_NS_alpha-sub001_counting_alpha;
sub001_memory_alpha_diff = sub001_NS_alpha-sub001_memory_alpha;

sub003_counting_alpha = log10(mean(squeeze(mean(powerPerSource{1}(:,:,9:14),3))));
sub003_memory_alpha = log10(mean(squeeze(mean(powerPerSource{2}(:,:,9:14),3))));
sub003_NS_alpha = log10(mean([squeeze(mean(powerPerSource{3}(:,:,9:14),3)); squeeze(mean(powerPerSource{4}(:,:,9:14),3)); squeeze(mean(powerPerSource{5}(:,:,9:14),3)); squeeze(mean(powerPerSource{6}(:,:,9:14),3))]));

sub003_counting_alpha_diff = sub003_NS_alpha-sub003_counting_alpha;
sub003_memory_alpha_diff = sub003_NS_alpha-sub003_memory_alpha;

sub010_counting_alpha = log10(mean([squeeze(mean(powerPerSource{12}(:,:,9:14),3)); squeeze(mean(powerPerSource{14}(:,:,9:14),3))]));
sub010_memory_alpha = log10(mean([squeeze(mean(powerPerSource{13}(:,:,9:14),3)); squeeze(mean(powerPerSource{15}(:,:,9:14),3))]));
sub010_NS_alpha = log10(mean([squeeze(mean(powerPerSource{17}(:,:,9:14),3)); squeeze(mean(powerPerSource{18}(:,:,9:14),3)); squeeze(mean(powerPerSource{19}(:,:,9:14),3))]));

sub010_counting_alpha_diff = sub010_NS_alpha-sub010_counting_alpha;
sub010_memory_alpha_diff = sub010_NS_alpha-sub010_memory_alpha;

sub029_counting_alpha = log10(mean([squeeze(mean(powerPerSource{20}(:,:,9:14),3)); squeeze(mean(powerPerSource{22}(:,:,9:14),3))]));
sub029_memory_alpha = log10(mean([squeeze(mean(powerPerSource{21}(:,:,9:14),3)); squeeze(mean(powerPerSource{23}(:,:,9:14),3))]));
sub029_NS_alpha = log10(mean([squeeze(mean(powerPerSource{24}(:,:,9:14),3)); squeeze(mean(powerPerSource{25}(:,:,9:14),3))]));

sub029_counting_alpha_diff = sub029_NS_alpha-sub029_counting_alpha;
sub029_memory_alpha_diff = sub029_NS_alpha-sub029_memory_alpha;

sub034_counting_alpha = log10(mean([squeeze(mean(powerPerSource{26}(:,:,9:14),3)); squeeze(mean(powerPerSource{28}(:,:,9:14),3))]));
sub034_memory_alpha = log10(mean([squeeze(mean(powerPerSource{27}(:,:,9:14),3)); squeeze(mean(powerPerSource{29}(:,:,9:14),3))]));
sub034_NS_alpha = log10(mean([squeeze(mean(powerPerSource{30}(:,:,9:14),3)); squeeze(mean(powerPerSource{31}(:,:,9:14),3))]));

sub034_counting_alpha_diff = sub034_NS_alpha-sub034_counting_alpha;
sub034_memory_alpha_diff = sub034_NS_alpha-sub034_memory_alpha;

min_val = min([sub001_counting_alpha_diff sub001_memory_alpha_diff sub003_counting_alpha_diff sub003_memory_alpha_diff sub010_counting_alpha_diff sub010_memory_alpha_diff sub029_counting_alpha_diff sub029_memory_alpha_diff sub034_counting_alpha_diff sub034_memory_alpha_diff]);
max_val = max([sub001_counting_alpha_diff sub001_memory_alpha_diff sub003_counting_alpha_diff sub003_memory_alpha_diff sub010_counting_alpha_diff sub010_memory_alpha_diff sub029_counting_alpha_diff sub029_memory_alpha_diff sub034_counting_alpha_diff sub034_memory_alpha_diff]);

%% Topography of alpha power - sub001

sub001_counting_alpha_surf = parcel_to_surface(sub001_counting_alpha_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub001_counting_alpha_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub001_memory_alpha_surf = parcel_to_surface(sub001_memory_alpha_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub001_memory_alpha_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

%% Topography of alpha power - sub003

sub003_counting_alpha_surf = parcel_to_surface(sub003_counting_alpha_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure, 
plot_cortical(sub003_counting_alpha_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub003_memory_alpha_surf = parcel_to_surface(sub003_memory_alpha_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub003_memory_alpha_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

%% Topography of alpha power - sub010

sub010_counting_alpha_surf = parcel_to_surface(sub010_counting_alpha_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub010_counting_alpha_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')
% hold on

sub010_memory_alpha_surf = parcel_to_surface(sub010_memory_alpha_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub010_memory_alpha_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')
% hold on

%% Topography of alpha power - sub029

sub029_counting_alpha_surf = parcel_to_surface(sub029_counting_alpha_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub029_counting_alpha_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub029_memory_alpha_surf = parcel_to_surface(sub029_memory_alpha_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub029_memory_alpha_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

%% Topography of alpha power - sub034

sub034_counting_alpha_surf = parcel_to_surface(sub034_counting_alpha_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub034_counting_alpha_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

sub034_memory_alpha_surf = parcel_to_surface(sub034_memory_alpha_diff, 'schaefer_100_fsa5');
% Project the results on the surface brain
f = figure,
plot_cortical(sub034_memory_alpha_surf, 'surface_name', 'fsa5', 'color_range', ...
              [-max_val max_val], 'cmap', 'RdBu_r')

%% Statistics - source-space - individual frequencies - within-subject 

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
    
    sub001_NS = power{11}(:,fr);
    sub001_count = [power{7}(:,fr); power{9}(:,fr)];
    pval_sub001_count(fr-1) = permutation_test(sub001_NS, sub001_count, n_perm);
    sub001_memory = [power{8}(:,fr); power{10}(:,fr)];
    pval_sub001_memory(fr-1) = permutation_test(sub001_NS, sub001_memory, n_perm);

    sub003_NS = [power{3}(:,fr); power{4}(:,fr); power{5}(:,fr); power{6}(:,fr)];
    sub003_count = power{1}(:,fr);
    pval_sub003_count(fr-1) = permutation_test(sub003_NS, sub003_count, n_perm);
    sub003_memory = power{2}(:,fr);
    pval_sub003_memory(fr-1) = permutation_test(sub003_NS, sub003_memory, n_perm);
    
    sub010_NS = [power{17}(:,fr); power{18}(:,fr); power{19}(:,fr)];
    sub010_count = [power{12}(:,fr); power{14}(:,fr)];
    pval_sub010_count(fr-1) = permutation_test(sub010_NS, sub010_count, n_perm);
    sub010_memory = [power{13}(:,fr); power{15}(:,fr)];
    pval_sub010_memory(fr-1) = permutation_test(sub010_NS, sub010_memory, n_perm);
    
    sub029_NS = [power{24}(:,fr); power{25}(:,fr)];
    sub029_count = [power{20}(:,fr); power{22}(:,fr)];
    pval_sub029_count(fr-1) = permutation_test(sub029_NS, sub029_count, n_perm);
    sub029_memory = [power{21}(:,fr); power{23}(:,fr)];
    pval_sub029_memory(fr-1) = permutation_test(sub029_NS, sub029_memory, n_perm);

    sub034_NS = [power{30}(:,fr); power{31}(:,fr)];
    sub034_count = [power{26}(:,fr); power{28}(:,fr)];
    pval_sub034_count(fr-1) = permutation_test(sub034_NS, sub034_count, n_perm);
    sub034_memory = [power{27}(:,fr); power{29}(:,fr)];
    pval_sub034_memory(fr-1) = permutation_test(sub034_NS, sub034_memory, n_perm);
end

[~,~,~,adj_pval_sub001] = fdr_bh_2015([pval_sub001_count; pval_sub001_memory]);
[~,~,~,adj_pval_sub003] = fdr_bh_2015([pval_sub003_count; pval_sub003_memory]);
[~,~,~,adj_pval_sub010] = fdr_bh_2015([pval_sub010_count; pval_sub010_memory]);
[~,~,~,adj_pval_sub029] = fdr_bh_2015([pval_sub029_count; pval_sub029_memory]);
[~,~,~,adj_pval_sub034] = fdr_bh_2015([pval_sub034_count; pval_sub034_memory]);

sub003_NSvscount = zeros(length(find(adj_pval_sub003(1:50)<0.05)),1);
sub003_count_sig = find(adj_pval_sub003(1:50)<0.05);
for i = 1:length(sub003_NSvscount)
    sub003_NS = [power{3}(:,sub003_count_sig(i)+1); power{4}(:,sub003_count_sig(i)+1); power{5}(:,sub003_count_sig(i)+1); power{6}(:,sub003_count_sig(i)+1)];
    sub003_count = power{1}(:,sub003_count_sig(i)+1);
    sub003_NSvscount(i) = mean(sub003_NS)-mean(sub003_count);
end

sub003_NSvsmemory = zeros(length(find(adj_pval_sub003(51:100)<0.05)),1);
sub003_memory_sig = find(adj_pval_sub003(51:100)<0.05);
for i = 1:length(sub003_NSvsmemory)
    sub003_NS = [power{3}(:,sub003_memory_sig(i)+1); power{4}(:,sub003_memory_sig(i)+1); power{5}(:,sub003_memory_sig(i)+1); power{6}(:,sub003_memory_sig(i)+1)];
    sub003_memory = power{2}(:,sub003_memory_sig(i)+1);
    sub003_NSvsmemory(i) = mean(sub003_NS)-mean(sub003_memory);
end

sub010_NSvscount = zeros(length(find(adj_pval_sub010(1:50)<0.05)),1);
sub010_count_sig = find(adj_pval_sub010(1:50)<0.05);
for i = 1:length(sub010_NSvscount)
    sub010_NS = [power{17}(:,sub010_count_sig(i)+1); power{18}(:,sub010_count_sig(i)+1); power{19}(:,sub010_count_sig(i)+1)];
    sub010_count = [power{12}(:,sub010_count_sig(i)+1); power{14}(:,sub010_count_sig(i)+1)];
    sub010_NSvscount(i) = mean(sub010_NS)-mean(sub010_count);
end

sub010_NSvsmemory = zeros(length(find(adj_pval_sub010(51:100)<0.05)),1);
sub010_memory_sig = find(adj_pval_sub010(51:100)<0.05);
for i = 1:length(sub010_NSvsmemory)
    sub010_NS = [power{17}(:,sub010_memory_sig(i)+1); power{18}(:,sub010_memory_sig(i)+1); power{19}(:,sub010_memory_sig(i)+1)];
    sub010_memory = [power{13}(:,sub010_memory_sig(i)+1); power{15}(:,sub010_memory_sig(i)+1)];
    sub010_NSvsmemory(i) = mean(sub010_NS)-mean(sub010_memory);
end

sub029_NSvscount = zeros(length(find(adj_pval_sub029(1:50)<0.05)),1);
sub029_count_sig = find(adj_pval_sub029(1:50)<0.05);
for i = 1:length(sub029_NSvscount)
    sub029_NS = [power{24}(:,sub029_count_sig(i)+1); power{25}(:,sub029_count_sig(i)+1)];
    sub029_count = [power{20}(:,sub029_count_sig(i)+1); power{22}(:,sub029_count_sig(i)+1)];
    sub029_NSvscount(i) = mean(sub029_NS)-mean(sub029_count);
end

sub029_NSvsmemory = zeros(length(find(adj_pval_sub029(51:100)<0.05)),1);
sub029_memory_sig = find(adj_pval_sub029(51:100)<0.05);
for i = 1:length(sub029_NSvsmemory)
    sub029_NS = [power{24}(:,sub029_memory_sig(i)+1); power{25}(:,sub029_memory_sig(i)+1)];
    sub029_memor = [power{21}(:,sub029_memory_sig(i)+1); power{23}(:,sub029_memory_sig(i)+1)];
    sub029_NSvsmemory(i) = mean(sub029_NS)-mean(sub029_memory);
end

sub034_NSvscount = zeros(length(find(adj_pval_sub034(1:50)<0.05)),1);
sub034_count_sig = find(adj_pval_sub034(1:50)<0.05);
for i = 1:length(sub034_NSvscount)
    sub034_NS = [power{30}(:,sub034_count_sig(i)+1); power{31}(:,sub034_count_sig(i)+1)];
    sub034_count = [power{26}(:,sub034_count_sig(i)+1); power{28}(:,sub034_count_sig(i)+1)];
    sub034_NSvscount(i) = mean(sub034_NS)-mean(sub034_count);
end

sub034_NSvsmemory = zeros(length(find(adj_pval_sub034(51:100)<0.05)),1);
sub034_memory_sig = find(adj_pval_sub034(51:100)<0.05);
for i = 1:length(sub034_NSvsmemory)
    sub034_NS = [power{30}(:,sub034_memory_sig(i)+1); power{31}(:,sub034_memory_sig(i)+1)];
    sub034_memory = [power{27}(:,sub034_memory_sig(i)+1); power{29}(:,sub034_memory_sig(i)+1)];
    sub034_NSvsmemory(i) = mean(sub034_NS)-mean(sub034_memory);
end


%% Statistics - source-space - regional alpha - within-subject

n_perm = 5000;
num_reg = 100;

pval_sub001_count_alpha = zeros(num_reg,1);
pval_sub001_memory_alpha = zeros(num_reg,1);
pval_sub003_count_alpha = zeros(num_reg,1);
pval_sub003_memory_alpha = zeros(num_reg,1);
pval_sub010_count_alpha = zeros(num_reg,1);
pval_sub010_memory_alpha = zeros(num_reg,1);
pval_sub029_count_alpha = zeros(num_reg,1);
pval_sub029_memory_alpha = zeros(num_reg,1);
pval_sub034_count_alpha = zeros(num_reg,1);
pval_sub034_memory_alpha = zeros(num_reg,1);

for reg = 1:num_reg
    sub001_NS = squeeze(mean(powerPerSource{11}(:,reg,9:14),3));
    sub001_count = [squeeze(mean(powerPerSource{7}(:,reg,9:14),3)); squeeze(mean(powerPerSource{9}(:,reg,9:14),3))];
    pval_sub001_count_alpha(reg) = permutation_test(sub001_NS, sub001_count, n_perm);
    sub001_memory = [squeeze(mean(powerPerSource{8}(:,reg,9:14),3)); squeeze(mean(powerPerSource{10}(:,reg,9:14),3))];
    pval_sub001_memory_alpha(reg) = permutation_test(sub001_NS, sub001_memory, n_perm);

    sub003_NS = [squeeze(mean(powerPerSource{3}(:,reg,9:14),3)); squeeze(mean(powerPerSource{4}(:,reg,9:14),3)); squeeze(mean(powerPerSource{5}(:,reg,9:14),3)); squeeze(mean(powerPerSource{6}(:,reg,9:14),3))];
    sub003_count = squeeze(mean(powerPerSource{1}(:,reg,9:14),3));
    pval_sub003_count_alpha(reg) = permutation_test(sub003_NS, sub003_count, n_perm);
    sub003_memory = squeeze(mean(powerPerSource{2}(:,reg,9:14),3));
    pval_sub003_memory_alpha(reg) = permutation_test(sub003_NS, sub003_memory, n_perm);
    
    sub010_NS = [squeeze(mean(powerPerSource{17}(:,reg,9:14),3)); squeeze(mean(powerPerSource{18}(:,reg,9:14),3)); squeeze(mean(powerPerSource{19}(:,reg,9:14),3))];
    sub010_count = [squeeze(mean(powerPerSource{12}(:,reg,9:14),3)); squeeze(mean(powerPerSource{14}(:,reg,9:14),3))];
    pval_sub010_count_alpha(reg) = permutation_test(sub010_NS, sub010_count, n_perm);
    sub010_memory = [squeeze(mean(powerPerSource{13}(:,reg,9:14),3)); squeeze(mean(powerPerSource{15}(:,reg,9:14),3))];
    pval_sub010_memory_alpha(reg) = permutation_test(sub010_NS, sub010_memory, n_perm);

    sub029_NS = [squeeze(mean(powerPerSource{24}(:,reg,9:14),3)); squeeze(mean(powerPerSource{25}(:,reg,9:14),3))];
    sub029_count = [squeeze(mean(powerPerSource{20}(:,reg,9:14),3)); squeeze(mean(powerPerSource{22}(:,reg,9:14),3))];
    pval_sub029_count_alpha(reg) = permutation_test(sub029_NS, sub029_count, n_perm);
    sub029_memory = [squeeze(mean(powerPerSource{21}(:,reg,9:14),3)); squeeze(mean(powerPerSource{23}(:,reg,9:14),3))];
    pval_sub029_memory_alpha(reg) = permutation_test(sub029_NS, sub029_memory, n_perm);

    sub034_NS = [squeeze(mean(powerPerSource{30}(:,reg,9:14),3)); squeeze(mean(powerPerSource{31}(:,reg,9:14),3))];
    sub034_count = [squeeze(mean(powerPerSource{26}(:,reg,9:14),3)); squeeze(mean(powerPerSource{28}(:,reg,9:14),3))];
    pval_sub034_count_alpha(reg) = permutation_test(sub034_NS, sub034_count, n_perm);
    sub034_memory = [squeeze(mean(powerPerSource{27}(:,reg,9:14),3)); squeeze(mean(powerPerSource{29}(:,reg,9:14),3))];
    pval_sub034_memory_alpha(reg) = permutation_test(sub034_NS, sub034_memory, n_perm);    
end

[~,~,~,adj_pval_sub001_alpha] = fdr_bh_2015([pval_sub001_count_alpha; pval_sub001_memory_alpha]);
[~,~,~,adj_pval_sub003_alpha] = fdr_bh_2015([pval_sub003_count_alpha; pval_sub003_memory_alpha]);
[~,~,~,adj_pval_sub010_alpha] = fdr_bh_2015([pval_sub010_count_alpha; pval_sub010_memory_alpha]);
[~,~,~,adj_pval_sub029_alpha] = fdr_bh_2015([pval_sub029_count_alpha; pval_sub029_memory_alpha]);
[~,~,~,adj_pval_sub034_alpha] = fdr_bh_2015([pval_sub034_count_alpha; pval_sub034_memory_alpha]);

yeo_indices = {[(1:9)'; (51:58)'], [(10:15)'; (59:66)'], [(16:23)'; (67:73)'], [(24:30)'; (74:78)'], [(31:33)'; (79:80)'], [(34:37)'; (81:89)'], [(38:50)'; (90:100)']};

sub001_count_alpha_sig_regs = find(adj_pval_sub001_alpha(1:100)<0.05);
sub001_count_alpha_sig_regs_direction = zeros(length(sub001_count_alpha_sig_regs),1);
for reg = 1:length(sub001_count_alpha_sig_regs)
    sub001_NS = squeeze(mean(powerPerSource{11}(:,sub001_count_alpha_sig_regs(reg),9:14),3));
    sub001_count = [squeeze(mean(powerPerSource{7}(:,sub001_count_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{9}(:,sub001_count_alpha_sig_regs(reg),9:14),3))];
    sub001_count_alpha_sig_regs_direction(reg) = sign(mean(sub001_count)-mean(sub001_NS));
end
sub001_count_alpha_yeo = zeros(7,1);
sub001_count_alpha_yeo_increase = zeros(7,1);
sub001_count_alpha_yeo_decrease = zeros(7,1);
for reg = 1:length(sub001_count_alpha_sig_regs)
    for y = 1:7
        if ismember(sub001_count_alpha_sig_regs(reg), yeo_indices{y})
            sub001_count_alpha_yeo(y) = sub001_count_alpha_yeo(y) + 1;
            if sub001_count_alpha_sig_regs_direction(reg) == 1
                sub001_count_alpha_yeo_decrease(y) = sub001_count_alpha_yeo_decrease(y) + 1;
            else
                sub001_count_alpha_yeo_increase(y) = sub001_count_alpha_yeo_increase(y) + 1;
            end
        end
    end
end

sub001_memory_alpha_sig_regs = find(adj_pval_sub001_alpha(101:200)<0.05);
sub001_memory_alpha_sig_regs_direction = zeros(length(sub001_memory_alpha_sig_regs),1);
for reg = 1:length(sub001_memory_alpha_sig_regs)
    sub001_NS = squeeze(mean(powerPerSource{11}(:,sub001_memory_alpha_sig_regs(reg),9:14),3));
    sub001_memory = [squeeze(mean(powerPerSource{8}(:,sub001_memory_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{10}(:,sub001_memory_alpha_sig_regs(reg),9:14),3))];
    sub001_memory_alpha_sig_regs_direction(reg) = sign(mean(sub001_memory)-mean(sub001_NS));
end
sub001_memory_alpha_yeo = zeros(7,1);
sub001_memory_alpha_yeo_increase = zeros(7,1);
sub001_memory_alpha_yeo_decrease = zeros(7,1);
for reg = 1:length(sub001_memory_alpha_sig_regs)
    for y = 1:7
        if ismember(sub001_memory_alpha_sig_regs(reg), yeo_indices{y})
            sub001_memory_alpha_yeo(y) = sub001_memory_alpha_yeo(y) + 1;
            if sub001_memory_alpha_sig_regs_direction(reg) == 1
                sub001_memory_alpha_yeo_decrease(y) = sub001_memory_alpha_yeo_decrease(y) + 1;
            else
                sub001_memory_alpha_yeo_increase(y) = sub001_memory_alpha_yeo_increase(y) + 1;
            end
        end
    end
end

sub003_count_alpha_sig_regs = find(adj_pval_sub003_alpha(1:100)<0.05);
sub003_count_alpha_sig_regs_direction = zeros(length(sub003_count_alpha_sig_regs),1);
sub003_count_alpha_sig_regs_diff = zeros(length(sub003_count_alpha_sig_regs),1);
for reg = 1:length(sub003_count_alpha_sig_regs)
    sub003_NS = [squeeze(mean(powerPerSource{3}(:,sub003_count_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{4}(:,sub003_count_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{5}(:,sub003_count_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{6}(:,sub003_count_alpha_sig_regs(reg),9:14),3))];
    sub003_count = squeeze(mean(powerPerSource{1}(:,sub003_count_alpha_sig_regs(reg),9:14),3));
    sub003_count_alpha_sig_regs_direction(reg) = sign(mean(sub003_count)-mean(sub003_NS));
    sub003_count_alpha_sig_regs_diff(reg) = mean(sub003_count)-mean(sub003_NS);
end
sub003_count_alpha_yeo = zeros(7,1);
sub003_count_alpha_yeo_increase = zeros(7,1);
sub003_count_alpha_yeo_decrease = zeros(7,1);
for reg = 1:length(sub003_count_alpha_sig_regs)
    for y = 1:7
        if ismember(sub003_count_alpha_sig_regs(reg), yeo_indices{y})
            sub003_count_alpha_yeo(y) = sub003_count_alpha_yeo(y) + 1;
            if sub003_count_alpha_sig_regs_direction(reg) == 1
                sub003_count_alpha_yeo_decrease(y) = sub003_count_alpha_yeo_decrease(y) + 1;
            else
                sub003_count_alpha_yeo_increase(y) = sub003_count_alpha_yeo_increase(y) + 1;
            end
        end
    end
end

sub003_memory_alpha_sig_regs = find(adj_pval_sub003_alpha(101:200)<0.05);
sub003_memory_alpha_sig_regs_direction = zeros(length(sub003_memory_alpha_sig_regs),1);
sub003_memory_alpha_sig_regs_diff = zeros(length(sub003_memory_alpha_sig_regs),1);
for reg = 1:length(sub003_memory_alpha_sig_regs)
    sub003_NS = [squeeze(mean(powerPerSource{3}(:,sub003_memory_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{4}(:,sub003_memory_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{5}(:,sub003_memory_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{6}(:,sub003_memory_alpha_sig_regs(reg),9:14),3))];
    sub003_memory = squeeze(mean(powerPerSource{2}(:,sub003_memory_alpha_sig_regs(reg),9:14),3));
    sub003_memory_alpha_sig_regs_direction(reg) = sign(mean(sub003_memory)-mean(sub003_NS));
    sub003_memory_alpha_sig_regs_diff(reg) = mean(sub003_memory)-mean(sub003_NS);
end
sub003_memory_alpha_yeo = zeros(7,1);
sub003_memory_alpha_yeo_increase = zeros(7,1);
sub003_memory_alpha_yeo_decrease = zeros(7,1);
for reg = 1:length(sub003_memory_alpha_sig_regs)
    for y = 1:7
        if ismember(sub003_memory_alpha_sig_regs(reg), yeo_indices{y})
            sub003_memory_alpha_yeo(y) = sub003_memory_alpha_yeo(y) + 1;
            if sub003_memory_alpha_sig_regs_direction(reg) == 1
                sub003_memory_alpha_yeo_decrease(y) = sub003_memory_alpha_yeo_decrease(y) + 1;
            else
                sub003_memory_alpha_yeo_increase(y) = sub003_memory_alpha_yeo_increase(y) + 1;
            end
        end
    end
end

sub010_count_alpha_sig_regs = find(adj_pval_sub010_alpha(1:100)<0.05);
sub010_count_alpha_sig_regs_direction = zeros(length(sub010_count_alpha_sig_regs),1);
for reg = 1:length(sub010_count_alpha_sig_regs)
    sub010_NS = [squeeze(mean(powerPerSource{17}(:,sub010_count_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{18}(:,sub010_count_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{19}(:,sub010_count_alpha_sig_regs(reg),9:14),3))];
    sub010_count = [squeeze(mean(powerPerSource{12}(:,sub010_count_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{14}(:,sub010_count_alpha_sig_regs(reg),9:14),3))];
    sub010_count_alpha_sig_regs_direction(reg) = sign(mean(sub010_count)-mean(sub010_NS));
end
sub010_count_alpha_yeo = zeros(7,1);
sub010_count_alpha_yeo_increase = zeros(7,1);
sub010_count_alpha_yeo_decrease = zeros(7,1);
for reg = 1:length(sub010_count_alpha_sig_regs)
    for y = 1:7
        if ismember(sub010_count_alpha_sig_regs(reg), yeo_indices{y})
            sub010_count_alpha_yeo(y) = sub010_count_alpha_yeo(y) + 1;
            if sub010_count_alpha_sig_regs_direction(reg) == 1
                sub010_count_alpha_yeo_decrease(y) = sub010_count_alpha_yeo_decrease(y) + 1;
            else
                sub010_count_alpha_yeo_increase(y) = sub010_count_alpha_yeo_increase(y) + 1;
            end
        end
    end
end

sub010_memory_alpha_sig_regs = find(adj_pval_sub010_alpha(101:200)<0.05);
sub010_memory_alpha_sig_regs_direction = zeros(length(sub010_memory_alpha_sig_regs),1);
for reg = 1:length(sub010_memory_alpha_sig_regs)
    sub010_NS = [squeeze(mean(powerPerSource{17}(:,sub010_memory_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{18}(:,sub010_memory_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{19}(:,sub010_memory_alpha_sig_regs(reg),9:14),3))];
    sub010_memory = [squeeze(mean(powerPerSource{13}(:,sub010_memory_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{15}(:,sub010_memory_alpha_sig_regs(reg),9:14),3))];
    sub010_memory_alpha_sig_regs_direction(reg) = sign(mean(sub010_memory)-mean(sub010_NS));
end
sub010_memory_alpha_yeo = zeros(7,1);
sub010_memory_alpha_yeo_increase = zeros(7,1);
sub010_memory_alpha_yeo_decrease = zeros(7,1);
for reg = 1:length(sub010_memory_alpha_sig_regs)
    for y = 1:7
        if ismember(sub010_memory_alpha_sig_regs(reg), yeo_indices{y})
            sub010_memory_alpha_yeo(y) = sub010_memory_alpha_yeo(y) + 1;
            if sub010_memory_alpha_sig_regs_direction(reg) == 1
                sub010_memory_alpha_yeo_decrease(y) = sub010_memory_alpha_yeo_decrease(y) + 1;
            else
                sub010_memory_alpha_yeo_increase(y) = sub010_memory_alpha_yeo_increase(y) + 1;
            end
        end
    end
end

sub029_count_alpha_sig_regs = find(adj_pval_sub029_alpha(1:100)<0.05);
sub029_count_alpha_sig_regs_direction = zeros(length(sub029_count_alpha_sig_regs),1);
for reg = 1:length(sub029_count_alpha_sig_regs)
    sub029_NS = [squeeze(mean(powerPerSource{24}(:,sub029_count_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{25}(:,sub029_count_alpha_sig_regs(reg),9:14),3))];
    sub029_count = [squeeze(mean(powerPerSource{20}(:,sub029_count_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{22}(:,sub029_count_alpha_sig_regs(reg),9:14),3))];
    sub029_count_alpha_sig_regs_direction(reg) = sign(mean(sub029_count)-mean(sub029_NS));
end
sub029_count_alpha_yeo = zeros(7,1);
sub029_count_alpha_yeo_increase = zeros(7,1);
sub029_count_alpha_yeo_decrease = zeros(7,1);
for reg = 1:length(sub029_count_alpha_sig_regs)
    for y = 1:7
        if ismember(sub029_count_alpha_sig_regs(reg), yeo_indices{y})
            sub029_count_alpha_yeo(y) = sub029_count_alpha_yeo(y) + 1;
            if sub029_count_alpha_sig_regs_direction(reg) == 1
                sub029_count_alpha_yeo_decrease(y) = sub029_count_alpha_yeo_decrease(y) + 1;
            else
                sub029_count_alpha_yeo_increase(y) = sub029_count_alpha_yeo_increase(y) + 1;
            end
        end
    end
end

sub029_memory_alpha_sig_regs = find(adj_pval_sub029_alpha(101:200)<0.05);
sub029_memory_alpha_sig_regs_direction = zeros(length(sub029_memory_alpha_sig_regs),1);
for reg = 1:length(sub029_memory_alpha_sig_regs)
    sub029_NS = [squeeze(mean(powerPerSource{24}(:,sub029_memory_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{25}(:,sub029_memory_alpha_sig_regs(reg),9:14),3))];
    sub029_memory = [squeeze(mean(powerPerSource{21}(:,sub029_memory_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{23}(:,sub029_memory_alpha_sig_regs(reg),9:14),3))];
    sub029_memory_alpha_sig_regs_direction(reg) = sign(mean(sub029_memory)-mean(sub029_NS));
end
sub029_memory_alpha_yeo = zeros(7,1);
sub029_memory_alpha_yeo_increase = zeros(7,1);
sub029_memory_alpha_yeo_decrease = zeros(7,1);
for reg = 1:length(sub029_memory_alpha_sig_regs)
    for y = 1:7
        if ismember(sub029_memory_alpha_sig_regs(reg), yeo_indices{y})
            sub029_memory_alpha_yeo(y) = sub029_memory_alpha_yeo(y) + 1;
            if sub029_memory_alpha_sig_regs_direction(reg) == 1
                sub029_memory_alpha_yeo_decrease(y) = sub029_memory_alpha_yeo_decrease(y) + 1;
            else
                sub029_memory_alpha_yeo_increase(y) = sub029_memory_alpha_yeo_increase(y) + 1;
            end
        end
    end
end

sub034_count_alpha_sig_regs = find(adj_pval_sub034_alpha(1:100)<0.05);
sub034_count_alpha_sig_regs_direction = zeros(length(sub034_count_alpha_sig_regs),1);
for reg = 1:length(sub034_count_alpha_sig_regs)
    sub034_NS = [squeeze(mean(powerPerSource{30}(:,sub034_count_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{31}(:,sub034_count_alpha_sig_regs(reg),9:14),3))];
    sub034_count = [squeeze(mean(powerPerSource{26}(:,sub034_count_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{28}(:,sub034_count_alpha_sig_regs(reg),9:14),3))];
    sub034_count_alpha_sig_regs_direction(reg) = sign(mean(sub034_count)-mean(sub034_NS));
end
sub034_count_alpha_yeo = zeros(7,1);
sub034_count_alpha_yeo_increase = zeros(7,1);
sub034_count_alpha_yeo_decrease = zeros(7,1);
for reg = 1:length(sub034_count_alpha_sig_regs)
    for y = 1:7
        if ismember(sub034_count_alpha_sig_regs(reg), yeo_indices{y})
            sub034_count_alpha_yeo(y) = sub034_count_alpha_yeo(y) + 1;
            if sub034_count_alpha_sig_regs_direction(reg) == 1
                sub034_count_alpha_yeo_decrease(y) = sub034_count_alpha_yeo_decrease(y) + 1;
            else
                sub034_count_alpha_yeo_increase(y) = sub034_count_alpha_yeo_increase(y) + 1;
            end
        end
    end
end

sub034_memory_alpha_sig_regs = find(adj_pval_sub034_alpha(101:200)<0.05);
sub034_memory_alpha_sig_regs_direction = zeros(length(sub034_memory_alpha_sig_regs),1);
for reg = 1:length(sub034_memory_alpha_sig_regs)
    sub034_NS = [squeeze(mean(powerPerSource{30}(:,sub034_memory_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{31}(:,sub034_memory_alpha_sig_regs(reg),9:14),3))];
    sub034_memory = [squeeze(mean(powerPerSource{27}(:,sub034_memory_alpha_sig_regs(reg),9:14),3)); squeeze(mean(powerPerSource{29}(:,sub034_memory_alpha_sig_regs(reg),9:14),3))];
    sub034_memory_alpha_sig_regs_direction(reg) = sign(mean(sub034_memory)-mean(sub034_NS));
end
sub034_memory_alpha_yeo = zeros(7,1);
sub034_memory_alpha_yeo_increase = zeros(7,1);
sub034_memory_alpha_yeo_decrease = zeros(7,1);
for reg = 1:length(sub034_memory_alpha_sig_regs)
    for y = 1:7
        if ismember(sub034_memory_alpha_sig_regs(reg), yeo_indices{y})
            sub034_memory_alpha_yeo(y) = sub034_memory_alpha_yeo(y) + 1;
            if sub034_memory_alpha_sig_regs_direction(reg) == 1
                sub034_memory_alpha_yeo_decrease(y) = sub034_memory_alpha_yeo_decrease(y) + 1;
            else
                sub034_memory_alpha_yeo_increase(y) = sub034_memory_alpha_yeo_increase(y) + 1;
            end
        end
    end
end

all_alpha_sig_regs = {sub001_count_alpha_sig_regs, sub001_memory_alpha_sig_regs, sub003_count_alpha_sig_regs, sub003_memory_alpha_sig_regs, sub010_count_alpha_sig_regs, sub010_memory_alpha_sig_regs, sub029_count_alpha_sig_regs, sub029_memory_alpha_sig_regs, sub034_count_alpha_sig_regs, sub034_memory_alpha_sig_regs};
all_alpha_sig_regs_direction = {sub001_count_alpha_sig_regs_direction, sub001_memory_alpha_sig_regs_direction, sub003_count_alpha_sig_regs_direction, sub003_memory_alpha_sig_regs_direction, sub010_count_alpha_sig_regs_direction, sub010_memory_alpha_sig_regs_direction, sub029_count_alpha_sig_regs_direction, sub029_memory_alpha_sig_regs_direction, sub034_count_alpha_sig_regs_direction, sub034_memory_alpha_sig_regs_direction};
common_alpha_sig_regs = all_alpha_sig_regs{1};
for i = 2:length(all_alpha_sig_regs)
    common_alpha_sig_regs = intersect(common_alpha_sig_regs, all_alpha_sig_regs{i});
end
common_alpha_sig_regs_direction = zeros(length(common_alpha_sig_regs),length(all_alpha_sig_regs));
for i = 1:length(common_alpha_sig_regs)
    for j = 1:length(all_alpha_sig_regs)
        reg_idx = find(all_alpha_sig_regs{j}==common_alpha_sig_regs(i));
        common_alpha_sig_regs_direction(i,j) = all_alpha_sig_regs_direction{j}(reg_idx);
    end
end

all_alpha_values = vertcat(all_alpha_sig_regs{:});
[unique_alpha_vals, ~, idx] = unique(all_alpha_values);
alpha_counts = accumarray(idx, 1);
[alpha_counts_sorted, alpha_sort_idx] = sort(alpha_counts, 'descend');
most_common_alpha_vals = unique_alpha_vals(alpha_sort_idx);

most_common_alpha_vals_direction = zeros(length(most_common_alpha_vals),length(all_alpha_sig_regs));
for i = 1:length(most_common_alpha_vals)
    for j = 1:length(all_alpha_sig_regs)
        reg_idx = find(all_alpha_sig_regs{j}==most_common_alpha_vals(i));
        if isempty(reg_idx)
            most_common_alpha_vals_direction(i,j) = 0;
        else
            most_common_alpha_vals_direction(i,j) = all_alpha_sig_regs_direction{j}(reg_idx);
        end
    end
end




