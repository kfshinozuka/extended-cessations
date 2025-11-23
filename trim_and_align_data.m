% Aligns time and restricts all datasets to shared channels (in same order)
% data_list - cell array containing data from each subject
% data_all - FieldTrip struct with data from each subject

function data_all = trim_and_align_data(data_list)

    % Find shortest trial length
    min_len = min(cellfun(@(x) length(x.time{1}), data_list));

    % Define common time vector
    fs = unique(cellfun(@(x) x.fsample, data_list));
    if numel(fs) > 1
        error('Sampling rates are inconsistent.');
    end
    dt = 1 / fs(1);
    common_time = 0 + (0:min_len-1) * dt;

    % Find common channel labels
    all_labels = cellfun(@(x) x.label, data_list, 'UniformOutput', false);
    common_labels = all_labels{1};
    for i = 2:length(all_labels)
        common_labels = intersect(common_labels, all_labels{i}, 'stable');
    end
    if isempty(common_labels)
        error('No common channels across datasets.');
    end

    % Select and align each dataset
    for i = 1:length(data_list)
        cfg = [];
        cfg.channel = common_labels;
        data_list{i} = ft_selectdata(cfg, data_list{i});

        % Trim trial and align time
        data_list{i}.trial{1} = data_list{i}.trial{1}(:, 1:min_len);
        data_list{i}.time{1}  = common_time;
    end

    % Append safely
    data_all = ft_appenddata([], data_list{:});
end
