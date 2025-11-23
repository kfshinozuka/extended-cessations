% Compute spatial filters for MEG data, using beamformer
% meg_data - preprocessed data, result of ft_preprocessing 
% sourcemodel - positions at which the virtual sensors will be placed
% headmodel - volume conduction (aka head) model of the subject
% source - spatial filters

function [ source, meg_data ] = GetSourceAnalysis_meg_NS(meg_data, sourcemodel, headmodel)

addpath("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Scripts")

meg_data = concatenate_trials(meg_data);

% Calculate means and covariances for source analysis
cfg                   = [];
cfg.covariance        = 'yes';
cfg.covariancewindow  = 'all';
cfg.channel           = 'MEG';
cfg.keeptrial         = 'yes';    
tlock                 = ft_timelockanalysis(cfg, meg_data);
% tlock.elec.label      = new_labels; % update labels

% Prepare MEG leadfield
cfg                 = [];
cfg.grad            = tlock.grad;  % EEG electrode positions
cfg.headmodel       = headmodel; % Use both headmodels
cfg.reducerank      = 2; % Typically needed for MEG, not for EEG
cfg.channel         = 'MEG';
cfg.sourcemodel     = sourcemodel;  % Inverse-warped template
cfg.normalize       = 'yes';
cfg.normalizeparam  = 1;
lf                  = ft_prepare_leadfield(cfg);

% Run source analysis on the subject-specific inverse-warped template positions
cfg                    = [];
cfg.method             = 'lcmv';
cfg.sourcemodel        = lf;
cfg.headmodel          = headmodel;
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.lambda        = '5%';
cfg.lcmv.fixedori      = 'yes';
source  = ft_sourceanalysis(cfg, tlock);

