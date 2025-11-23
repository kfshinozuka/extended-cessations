% Concatenate all the trials in FieldTrip data

function data = concatenate_trials(data)

if length(data.trial) > 1
    total_timepoints = 0;
    trl = data.trial;
    for t = 1:length(trl)
        total_timepoints = total_timepoints + size(trl{t},2);
    end
    num_chans = length(data.label);
    data_concat = zeros(num_chans, total_timepoints);
    time_concat = zeros(1, total_timepoints);
    start_timepoint = 0;
    for t = 1:length(trl)
        end_timepoint = start_timepoint + size(trl{t},2);
        data_concat(:, start_timepoint+1:end_timepoint) = trl{t}; 
        time_concat(1, start_timepoint+1:end_timepoint) = data.time{t};
        start_timepoint = end_timepoint;
    end
    data.trial = {data_concat};
    data.time = {time_concat};
end