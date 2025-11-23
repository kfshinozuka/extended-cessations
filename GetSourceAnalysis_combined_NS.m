% Compute spatial filters jointly for EEG + MEG data, using beamformer
% eeg_data_preprocessed, meg_data_preprocessed - preprocessed data, result of ft_preprocessing
% sourcemodel - positions at which the virtual sensors will be placed
% eeg_headmodel, meg_headmodel - volume conduction (aka head) model of the subject
% source - spatial filters

function source = GetSourceAnalysis_combined_alt_NS(eeg_data_preprocessed, meg_data_preprocessed, sourcemodel, eeg_headmodel, meg_headmodel)

% Compute channel-wise standard deviation 
eeg_std = std(cell2mat(cellfun(@(x) x, eeg_data_preprocessed.trial, 'UniformOutput', false)), 0, 2);
meg_std = std(cell2mat(cellfun(@(x) x, meg_data_preprocessed.trial, 'UniformOutput', false)), 0, 2);

mean_eeg_std = mean(eeg_std);
mean_meg_std = mean(meg_std);

scaling_factor = mean_meg_std / mean_eeg_std;

% Scale EEG trials before merging
for t = 1:length(eeg_data_preprocessed.trial)
    eeg_data_preprocessed.trial{t} = eeg_data_preprocessed.trial{t} * scaling_factor;
end

eeg_data_concat = concatenate_trials(eeg_data_preprocessed);
meg_data_concat = concatenate_trials(meg_data_preprocessed);

% Combine EEG and MEG data
data = ft_appenddata([], eeg_data_concat, meg_data_concat);
data.grad = meg_data_concat.grad;
data.elec = eeg_data_concat.elec;

% Calculate means and covariances for source analysis
cfg                   = [];
cfg.covariance        = 'yes';
cfg.covariancewindow  = 'all';
cfg.channel           = {'EEG', 'MEG'};
cfg.keeptrial         = 'yes';
tlock                 = ft_timelockanalysis(cfg, data);
% tlock.grad            = ft_convert_units(tlock.grad, 'mm');

% Prepare EEG leadfield
cfg                 = [];
cfg.elec            = tlock.elec;  % EEG electrode positions
cfg.headmodel       = eeg_headmodel; % Use both headmodels
cfg.channel         = 'EEG';
cfg.sourcemodel     = sourcemodel;  % Inverse-warped template
cfg.normalize       = 'yes';
cfg.normalizeparam  = 1;
lf_eeg              = ft_prepare_leadfield(cfg);

% Prepare MEG leadfield
cfg                 = [];
cfg.grad            = tlock.grad;  % EEG electrode positions
cfg.headmodel       = meg_headmodel; % Use both headmodels
cfg.reducerank      = 2; % Typically needed for MEG, not for EEG
cfg.channel         = 'MEG';
cfg.sourcemodel     = sourcemodel;  % Inverse-warped template
cfg.normalize       = 'yes';
cfg.normalizeparam  = 1;
lf_meg              = ft_prepare_leadfield(cfg);

lf = lf_meg;
for i = 1:length(lf_meg.leadfield)
    lf.leadfield{i} = cat(1, lf_eeg.leadfield{i}, lf_meg.leadfield{i});
end
lf.label = [lf_eeg.label; lf_meg.label];

% Run source analysis on the subject-specific inverse-warped template positions
cfg                    = [];
cfg.method             = 'lcmv';
cfg.sourcemodel        = lf;
cfg.headmodel          = {eeg_headmodel, meg_headmodel};
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.lambda        = '5%';
cfg.lcmv.fixedori      = 'yes';
source  = ft_sourceanalysis(cfg, tlock);
