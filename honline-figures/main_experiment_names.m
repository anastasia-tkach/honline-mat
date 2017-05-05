%close all;
clear;

start_offset = 1;
half_window_size = 1;
folder_name = 'sridhar2';
folder_name = 'tompson_no_wrist';
folder_name = 'synthetic_02';
sequence_names = {'teaser_short'};

ylim_time_max = 5.0;
ylim_time_min = 2.3;
display_time = false;
display_statistics = true;
display_variances = false;
weighted = false;
full = false;
markers = false;
synthetic = true;

if (weighted)
    xlim_max = 3.0;
    xlim_min = 1.0;
end
if (full)
    xlim_max = 3.0;
    xlim_min = 5.2;
end
if (markers)
    xlim_max = 30;
    xlim_min = 0;
end
if (synthetic)
    xlim_max = inf;
    xlim_min = 0;
end

estimation_types = {'KALMAN_DIAGONAL_EVALUATION_0', 'KALMAN_STANDARD_EVALUATION_0', 'KALMAN_EXTENDED_EVALUATION_10', ...
    'ONLINE_EVALUATION_0', 'ONLINE_CALIBRATION_0', 'BATCH_OFFLINE_EVALUATION_10', 'BATCH_ONLINE_EVALUATION_10', 'ORIGINAL_10'};

%estimation_types = {'KALMAN_STANDARD_EVALUATION_0', 'KALMAN_EXTENDED_EVALUATION_10', 'ORIGINAL_10', 'ONLINE_CALIBRATION_0', 'BATCH_OFFLINE_EVALUATION_10', 'BATCH_ONLINE_EVALUATION_10'};

%estimation_types = {'marker_based_metrics_full_iters_15_time_5_wrist'};
%estimation_types = {'KALMAN_STANDARD_CALIBRATION_0', 'KALMAN_STANDARD_EVALUATION_0', 'KALMAN_EXTENDED_EVALUATION_10', 'BATCH_OFFLINE_EVALUATION_0', 'BATCH_ONLINE_EVALUATION_0'};
estimation_types = {'ONLINE_0'};

%% Processing
listing = dir(['E:\Data\honline-results\', folder_name]);
series_names = {};
for i = 3:length(listing)
    split = strsplit(listing(i).name,'_');
    series_names{end + 1} = [split{end - 1}, '_'];
end
series_names = unique(series_names);

experiment_names = {};

for s = 1:length(sequence_names)
    for i = 1:length(series_names)
        for j = 1:length(estimation_types)
            weight = 0;
            
            if (weighted)
                experiment_names{end + 1} = [estimation_types{j}, '_', sequence_names{s}, '_', series_names{i}, '_weighted'];
            end
            if (markers)
                 experiment_names{end + 1} = [estimation_types{j}, '_', sequence_names{s}, '_', series_names{i}, '_markers'];
            end
            if (full)
                experiment_names{end + 1} = [estimation_types{j}, '_', sequence_names{s}, '_', series_names{i}];
            end
            if (synthetic)
                experiment_names{end + 1} = [estimation_types{j}, '_', sequence_names{s}, '_', series_names{i}, '_solutions'];
            end
            
        end
    end
end

experiment_names = sort(experiment_names);
data_path = ['E:\Data\honline-results\', folder_name, '\'];
legend_names = cell(length(experiment_names), 1);
for i = 1:length(experiment_names)
    legend_names{i} = strrep(experiment_names{i}, '_', ' ');
end

%% Compute marker based metrics
if (markers)   
    main_marker_based_metrics;
end

%% Compute online metrics
if (full || weighted)
    main_online_metrics;
end

%% Compute synthetic metrics
if (synthetic)
    main_synthetic_metrics;
end