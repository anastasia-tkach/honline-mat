%close all;
clear;

start_offset = 1;
half_window_size = 1;

weighted = false;
full = true;
markers = false;
synthetic = false;
display_time = false;
display_statistics = false;
display_variances = false;

%% Synthetic data
folder_name = 'synthetic_anastasia_easy'; sequence_names = {'anastasia_easy'};
estimation_types = {'ONLINE_0.050000', 'ONLINE_0.075000', 'ONLINE_0.100000', 'ONLINE_0.125000', 'ONLINE_0.150000', 'ONLINE_0.175000', 'ONLINE_0.200000', 'ONLINE_0.225000', 'ONLINE_0.250000', 'ONLINE_0.275000', 'ONLINE_0.300000', 'ONLINE_0.325000', 'ONLINE_0.350000', 'ONLINE_0.375000', 'ONLINE_0.400000'};
synthetic = true;

%% Real data - Tompson
%folder_name = 'tompson_FINAL'; sequence_names = {'tompson'};
%estimation_types = {'KALMAN_STANDARD_EVALUATION_0.000000', 'KALMAN_EXTENDED_EVALUATION_0.000000', 'ONLINE_CALIBRATION_0.000000', 'BATCH_OFFLINE_EVALUATION_0.000000', 'BATCH_ONLINE_EVALUATION_0.000000', 'TEMPLATE_0.000000'};
%markers = true;

%% Real data - teaser
%folder_name = 'teaser_final_30_perturbed_02'; sequence_names = {'teaser'};
%estimation_types = {'BATCH_ONLINE_1_EVALUATION_0.200000', 'KALMAN_STANDARD_EVALUATION_0.200000', 'KALMAN_EXTENDED_EVALUATION_0.200000', 'ONLINE_CALIBRATION_0.200000', 'TEMPLATE_0.200000', 'TAYLOR', 'SHARP', 'HTRACK'};
%estimation_types = {'BATCH_ONLINE_1_EVALUATION_0.200000', 'KALMAN_STANDARD_EVALUATION_0.200000', 'KALMAN_EXTENDED_EVALUATION_0.200000', 'ONLINE_CALIBRATION_0.200000', 'TEMPLATE_0.200000', 'ORIGINAL_0.200000'};
%estimation_types = {'KALMAN_STANDARD_EVALUATION_0.000000', 'KALMAN_EXTENDED_EVALUATION_0.000000', 'ONLINE_CALIBRATION_0.000000', 'BATCH_OFFLINE_EVALUATION_0.000000', 'BATCH_ONLINE_EVALUATION_0.000000', 'TEMPLATE_0.000000'};
% full = true;
% weighted = false;

display_statistics = true;

ylim_time_max = 5.0;
ylim_time_min = 2.3;

if (weighted)
    xlim_max = 3.0;
    xlim_min = 1.0;
end
if (full)
    xlim_min = 3.0;
    xlim_max = 8.0;
end
if (markers)
    xlim_max = 30;
    xlim_min = 0;
end
if (synthetic)
    xlim_max = inf;
    xlim_min = 0;
end



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
for i = 1:length(experiment_names), legend_names{i} = strrep(experiment_names{i}, '_', ' '); end
legend_type_names = cell(length(estimation_types), 1);
for i = 1:length(estimation_types), legend_type_names{i} = strrep(estimation_types{i}, '_', ' '); end
colors = {[0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840], [0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290 0.6940 0.1250],  [0.77, 0.77, 0.77]};

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