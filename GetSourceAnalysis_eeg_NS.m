% Compute spatial filters for EEG data, using beamformer
% eeg_data - preprocessed data, result of ft_preprocessing 
% sourcemodel - positions at which the virtual sensors will be placed
% headmodel - volume conduction (aka head) model of the subject
% subid - subject ID
% source - spatial filters

function [ source, eeg_data ] = GetSourceAnalysis_eeg_NS(eeg_data, sourcemodel, headmodel, subid)

if strcmp(subid, 'sub001')
    elec = ft_read_sens('/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Data/Original/sub001/digitizer_data/sub_001_day1.sfp');   
    labels = str2double(eeg_data.label);
    eeg_data.label = elec.label(labels+3);
end

% Calculate means and covariances for source analysis
cfg                   = [];
cfg.covariance        = 'yes';
cfg.covariancewindow  = 'all';
cfg.method            = 'ridge';
if ~strcmp(subid, 'sub001')
    cfg.channel       = 'EEG';
end
cfg.keeptrial         = 'yes';    
tlock                 = ft_timelockanalysis(cfg, eeg_data);

% Prepare leadfield model
cfg                 = [];
if strcmp(subid, 'sub001')
    cfg.elec = elec;
else
    cfg.elec = tlock.elec;
end
cfg.headmodel       = headmodel;
cfg.reducerank      = 2;
if ~strcmp(subid, 'sub001')
    cfg.channel     = 'EEG';
end
cfg.sourcemodel     = sourcemodel;  % Inverse-warped template
cfg.normalize       = 'yes';
cfg.normalizeparam  = 1;
[lf]                = ft_prepare_leadfield(cfg);

% Run source analysis on the subject-specific inverse-warped template positions
cfg                    = [];
cfg.method             = 'lcmv';
cfg.sourcemodel        = lf;
cfg.headmodel          = headmodel;
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.lambda        = '20%';
cfg.lcmv.fixedori      = 'yes';
source  = ft_sourceanalysis(cfg, tlock);

