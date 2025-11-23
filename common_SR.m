% Compute joint EEG+MEG source reconstruction for sub003, sub010, sub034;
% MEG source reconstruction for sub029; EEG source reconstruction for
% sub001

raw_meeg_datapath = '/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Data/Original';
preprocessed_meeg_datapath = '/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Data/Preprocessed';
mri_datapath = '/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Data/MRI';
mri_realigned_datapath = '/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Data/MRI/Realigned';
elec_realigned_datapath = '/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Data/Elec';
headmodel_datapath = '/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Data/Headmodel';
resultspath = '/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Results/SourceRecon/Combined/Truncated/CommonFilter';
addpath('/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Oxford/Research/Rebirth/Sheaf/Scripts')
addpath("/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Scripts")

load('/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Data/schaefer100_centroids.mat')
positions = schaefer100_centroids;

d_mri = dir(fullfile(mri_datapath, '*T1w.nii'));

freqs = [1 4; 4 8; 8 13; 13 20; 20 30; 30 48];

for f_mri = 1:length(d_mri) % Loop through MRIs of the five subjects

    cd(mri_datapath)
    subid = d_mri(f_mri).name(1:6);
    mri = ft_read_mri(d_mri(f_mri).name);
    if strcmp(subid, 'sub001')
        mri.coordsys = 'ras';
    else
        mri.coordsys = 'acpc';
    end
    if strcmp(subid, 'sub034')
        mri_unbias = ft_volumebiascorrect([],mri);
    end

    if strcmp(subid, 'sub001') % EEG data

        % Transform MRI to head space (1st pass)
        cfg = [];
        cfg.method = 'interactive';
        cfg.coordsys = 'ctf';
        mri_realigned = ft_volumerealign(cfg, mri);  
        % load(fullfile(mri_realigned_datapath, [subid '_mri_realigned.mat']))

        % Source model in head space
        cfg = [];
        cfg.sourcemodel.pos = positions;
        template = ft_prepare_sourcemodel(cfg);
        sourcemodel = SourcemodelFromTemplate(template, mri_realigned);

        % Head model in head space
        % EEG
        cfg = [];
        cfg.output = {'brain', 'skull', 'scalp'};
        segmentedmri = ft_volumesegment(cfg, mri_realigned);
        % load(fullfile(mri_datapath, 'Segmented', [subid '_segmentedmri.mat']))        
        cfg = [];
        cfg.method = 'concentricspheres';
        headmodel = ft_prepare_headmodel(cfg, segmentedmri);
        % load(fullfile(headmodel_datapath, [subid '_eeg_headmodel.mat']))

        % Electrode realignment
        elec = ft_read_sens('/Users/kennethshinozuka/Macintosh HD - Data/Users/kshinozuka/Documents/Nirodha Samapatti/Data/Original/sub001/digitizer_data/sub_001_day1.sfp');
        elec = ft_convert_units(elec, 'mm');
        cfg = [];
        cfg.tissue      = 'scalp';
        scalp = ft_prepare_mesh(cfg, segmentedmri);
        cfg            = [];
        cfg.method     = 'interactive';
        cfg.elec       = elec;
        cfg.headshape  = scalp;
        elec_realigned = ft_electroderealign(cfg);
        % load(fullfile(elec_realigned_datapath, [subid '_elec_realigned.mat']))

        % Check alignment between electrodes, head model, and source model
        figure;
        ft_plot_sens(elec_realigned, 'style', 'ob');
        ft_plot_headmodel(headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
        ft_plot_mesh(sourcemodel.pos, 'vertexcolor', 'r');

        cd(fullfile(preprocessed_meeg_datapath, subid, 'Truncated'))
        d_eeg = dir('*auto.mat');

        % Compute common spatial filters on all of sub001's EEG data (all
        % conditions and all runs of each condition)
        eeg_data = {};
        for i = 1:numel(d_eeg)
            load(d_eeg(i).name)
            eeg_data{end+1} = concatenate_trials(eeg_data_preprocessed);
        end
        % Align time and restrict all datasets to shared channels (in same order)
        eeg_data_all = trim_and_align_data(eeg_data); 
        eeg_data_all.elec = elec_realigned;
        source_all = GetSourceAnalysis_eeg_NS(eeg_data_all, sourcemodel, headmodel, subid);
        for f_eeg = 1:length(eeg_data_all.trial)
            data = [];
            data.trial = {eeg_data_all.trial{f_eeg}};
            vs = GetVirtualSensors(data, source_all);
            save(fullfile(resultspath, [d_eeg(f_eeg).name(1:end-44) 'preprocessed_noln_nobt_auto_combined_simple_SR_realigned_ER_common.mat']), 'vs')
        end

    elseif strcmp(subid, 'sub029') % MEG data

        % Transform MRI to head space (1st pass - manual alignment)
        cfg = [];
        cfg.method = 'interactive';
        cfg.coordsys = 'neuromag';
        mri_realigned1 = ft_volumerealign(cfg, mri);

        cd(fullfile(raw_meeg_datapath, subid, 'meg_eeg', 'control_new'))
        d_raw_meeg = dir('*.fif');
        dataset = d_raw_meeg(1).name;
        shape = ft_read_headshape(dataset);

        % Transform MRI to head space (2nd pass - alignment using headshape)
        cfg = [];
        cfg.method = 'headshape';
        cfg.headshape = shape; 
        mri_realigned2 = ft_volumerealign(cfg, mri_realigned1);
        % load(fullfile(mri_datapath, 'Realigned', [subid '_mri_realigned2.mat']))

        % Source model in head space
        cfg = [];
        cfg.sourcemodel.pos = positions;
        template = ft_prepare_sourcemodel(cfg);
        sourcemodel = SourcemodelFromTemplate(template, mri_realigned2);

        % Head model in head space
        % MEG
        cfg = [];
        cfg.output = {'brain'};
        segmentedmri_meg = ft_volumesegment(cfg, mri_realigned2);
        cfg = [];
        cfg.method = 'singleshell';
        headmodel = ft_prepare_headmodel(cfg, segmentedmri_meg);
        load(fullfile(headmodel_datapath, [subid '_meg_headmodel.mat']))

        % Check alignment between head model, source model, and sensor
        % positions
        % MEG
        figure;
        grad = ft_read_sens(dataset, 'senstype', 'meg');
        grad = ft_convert_units(grad, 'mm');
        ft_plot_sens(grad, 'style', 'ob');
        ft_plot_headmodel(meg_headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
        ft_plot_mesh(sourcemodel.pos, 'vertexcolor', 'r');

        % Compute common spatial filters on all of sub029's MEG data (all
        % conditions and all runs of each condition)
        cd(fullfile(preprocessed_meeg_datapath, subid, 'Truncated'))
        d_meg = dir('*meg*auto.mat');
        meg_data = {};
        for i = 1:numel(d_meg)
            load(d_meg(i).name)
            meg_data{end+1} = concatenate_trials(meg_data_preprocessed);
        end
        % Align time and restrict all datasets to shared channels (in same order)
        meg_data_all = trim_and_align_data(meg_data);   
        meg_data_all.grad = meg_data{1}.grad;
        source_all = GetSourceAnalysis_meg_NS(meg_data_all, sourcemodel, headmodel);
        for f_meg = 1:length(meg_data_all.trial)
            data = [];
            data.trial = {meg_data_all.trial{f_meg}};
            vs = GetVirtualSensors(data, source_all);
            save(fullfile(resultspath, [d_meg(f_meg).name(1:end-4) 'preprocessed_noln_nobt_auto_SR_realigned_common.mat']), 'vs')
        end

    else % EEG/MEG data

        % Transform MRI to head space (1st pass - manual alignment)
        cfg = [];
        cfg.method = 'interactive';
        cfg.coordsys = 'neuromag';
        if strcmp(subid, 'sub034')
            mri_realigned1 = ft_volumerealign(cfg, mri_unbias);
        else
            mri_realigned1 = ft_volumerealign(cfg, mri);  
        end

        cd(fullfile(raw_meeg_datapath, subid, 'meg_eeg', 'control'))
        d_raw_meeg = dir('*.fif');
        dataset = d_raw_meeg(1).name;
        shape = ft_read_headshape(dataset);

        % Transform MRI to head space (2nd pass - alignment using headshape)
        cfg = [];
        cfg.method = 'headshape';
        cfg.headshape = shape; 
        mri_realigned2 = ft_volumerealign(cfg, mri_realigned1);
        % load(fullfile(mri_datapath, 'Realigned', [subid '_mri_realigned2.mat']))

        % Source model in head space
        cfg = [];
        cfg.sourcemodel.pos = positions;
        template = ft_prepare_sourcemodel(cfg);
        sourcemodel = SourcemodelFromTemplate(template, mri_realigned2);

        % Head model in head space
        % EEG
        cfg = [];
        cfg.output = {'brain', 'skull', 'scalp'};
        if strcmp(subid, 'sub034')
            cfg.spmversion = 'spm12';      % use SPM12 explicitly
            cfg.scalpsmooth = 'no';
            cfg.scalpthreshold = 0.1;
        end
        segmentedmri_eeg = ft_volumesegment(cfg, mri_realigned2);
        % load(fullfile(mri_datapath, 'Segmented', [subid '_segmentedmri_eeg.mat']))
        cfg = [];
        cfg.method = 'concentricspheres';
        eeg_headmodel = ft_prepare_headmodel(cfg, segmentedmri_eeg);
        % load(fullfile(headmodel_datapath, [subid '_eeg_headmodel.mat']))
        % MEG
        cfg = [];
        cfg.output = {'brain'};
        if strcmp(subid, 'sub034')
            cfg.spmversion = 'spm12';      % use SPM12 explicitly
        end
        segmentedmri_meg = ft_volumesegment(cfg, mri_realigned2);
        cfg = [];
        cfg.method = 'singleshell';
        meg_headmodel = ft_prepare_headmodel(cfg, segmentedmri_meg);
        % load(fullfile(headmodel_datapath, [subid '_meg_headmodel.mat']))

        % Electrode realignment
        elec = ft_read_sens(dataset, 'senstype', 'eeg');
        elec = ft_convert_units(elec, 'mm');
        cfg = [];
        cfg.tissue      = 'scalp';
        scalp          = ft_prepare_mesh(cfg, segmentedmri_eeg);
        cfg            = [];
        cfg.method     = 'interactive';
        cfg.elec       = elec;
        cfg.headshape  = scalp;
        elec_realigned = ft_electroderealign(cfg);
        % load(fullfile(elec_realigned_datapath, [subid '_elec_realigned.mat']))

        % Check alignment between head model, source model, and sensor
        % positions
        % EEG
        figure;
        ft_plot_sens(elec_realigned, 'style', 'ob');
        ft_plot_headmodel(eeg_headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
        ft_plot_mesh(sourcemodel.pos, 'vertexcolor', 'r');
        % MEG
        figure;
        grad = ft_read_sens(dataset, 'senstype', 'meg');
        grad = ft_convert_units(grad, 'mm');
        ft_plot_sens(grad, 'style', 'ob');
        ft_plot_headmodel(meg_headmodel, 'facealpha', 0.5, 'edgecolor', 'none');
        ft_plot_mesh(sourcemodel.pos, 'vertexcolor', 'r');

        % Compute common spatial filters on all of this subjec'ts EEG and MEG data (all
        % conditions and all runs of each condition)
        cd(fullfile(preprocessed_meeg_datapath, subid, 'Truncated'))
        d_eeg = dir('*truncatedraw_eeg_preprocessed_noln_nobt_auto_combined.mat');
        d_meg = dir('*truncatedraw_meg_preprocessed_noln_nobt_auto_combined.mat');
        eeg_data = {};
        for i = 1:numel(d_eeg)
            load(d_eeg(i).name)
            eeg_data{end+1} = concatenate_trials(eeg_data_preprocessed);
        end
        meg_data = {};
        for i = 1:numel(d_meg)
            load(d_meg(i).name)
            meg_data{end+1} = concatenate_trials(meg_data_preprocessed);
        end
        % Align time and restrict all datasets to shared channels (in same order)
        eeg_data_all = trim_and_align_data(eeg_data);
        eeg_data_all.elec = elec_realigned;
        meg_data_all = trim_and_align_data(meg_data);
        meg_data_all.grad = meg_data{1}.grad;
        source_all = GetSourceAnalysis_combined_alt_NS(eeg_data_all, meg_data_all, sourcemodel, eeg_headmodel, meg_headmodel);
        for f_eeg = 1:length(eeg_data_all.trial)
            data = [];
            data.trial = {[eeg_data_all.trial{f_eeg}; meg_data_all.trial{f_eeg}]};
            vs = GetVirtualSensors(data, source_all);
            save(fullfile(resultspath, [d_eeg(f_eeg).name(1:end-44) 'preprocessed_noln_nobt_auto_combined_simple_SR_realigned_ER_common.mat']), 'vs')
        end
    end
end
