%clear; close all; clc;
start_offset = 50;
line_width = 1.5;
display_title = true;
num_markers = 25;

dataset_type = 'tkach';
display_time_metrics = false;

if strcmp(dataset_type, 'tkach')
    data_path = 'E:\Data\honline-results\marker_positions_metrics\tkach_dataset\';
    num_frames = 1570;
    
     experiment_names = {...%'template',...
         'batch_unperturbed', ...%'batch_shape_prior_20_perturbed', 'batch_shape_prior_0_perturbed', ...%'batch_shape_prior_20_strongly_perturbed'...        
         ...%'kalman_shape_prior_20',  'kalman_shape_prior_0', ...
         ...%'kalman_shape_prior_20_R_10000', 'kalman_shape_prior_20_R_1000', 'kalman_shape_prior_20_R_100', 'kalman_shape_prior_20_R_10', 'kalman_shape_prior_20_R_1', ...
         ...%'kalman_shape_prior_20_R_1000_Q_0.1', 'kalman_shape_prior_20_R_1000_Q_0.01', 'kalman_shape_prior_20_R_1000_Q_0.05', 'kalman_shape_prior_20_R_1000_Q_0.001', 'kalman_shape_prior_20_R_1000', ...         ... %'kalman_shape_prior_20_R_100_covariance', ...         
         ...%'kalman_two_stages', 
         ...%'kalman_shape_prior_20_R_1000_Q_0.0001_I_5', 'kalman_shape_prior_20_R_1000_Q_0.001_I_15', 'kalman_shape_prior_20_R_1000_Q_0.001_I_20', 'kalman_shape_prior_20_R_1000_Q_0.001_I_30', 'kalman_shape_prior_20_R_1000_Q_0.001_I_30-2', 'kalman_shape_prior_20_R_1000_Q_0.001_I_60', 'kalman_shape_prior_20_R_1000_Q_0.001_I_120',...
         ...%'kalman_shape_prior_20_R_1000_Q_0.001_I_30', 'kalman_shape_prior_20_R_1000_Q_0.001_I_30_attract_5', 'kalman_shape_prior_20_R_1000_Q_0.001_I_30_attract_10', 'kalman_shape_prior_20_R_1000_Q_0.001_I_30_attract_100', 'kalman_shape_prior_20_R_1000_Q_0.001_I_30_attract_300', 'kalman_shape_prior_20_R_1000_Q_0.001_I_30_attract_1000',  'kalman_shape_prior_20_R_1000_Q_0.001_I_30_attract_10000'...
         'kalman_shape_prior_0_R_1000_Q_0.001_I_30_attract_0', 'kalman_shape_prior_0_R_1000_Q_0.001_I_30_attract_0_perturbed', 'kalman_shape_prior_20_R_1000_Q_0.001_I_30', 'kalman_shape_prior_0_R_1000_Q_0.001_I_30_attract_100_perturbed', ...
         ...%'independent_fingers_shape_prior_20', 'independent_fingers_shape_prior_0', ...'independent_fingers_shape_prior_100', ...
         ...%'marker_based_metrics',...
         };
end
if strcmp(dataset_type, 'tompson')
    data_path = 'E:\Data\honline-results\marker_positions_metrics\tompson_dataset\';
    num_frames = 2300;
    
    experiment_names = {'uniform_scaling', ...
        'batch_unperturbed', 'batch_shape_prior_20_perturbed', 'batch_shape_prior_0_perturbed', ...
        'kalman_shape_prior_20', 'kalman_shape_prior_0', ...%'kalman_shape_prior_0_R_1000_Q_0.001_I_5_attract_100', 'kalman_shape_prior_20_R_1000_Q_0.001_I_5_attract_100', ...
        'independent_fingers_shape_prior_20', 'independent_fingers_shape_prior_0', ...% 'independent_fingers_shape_prior_100'
        };
end

%% Rigid reweight


legend_names = cell(length(experiment_names), 1);
for i = 1:length(experiment_names)
    legend_names{i} = strrep(experiment_names{i}, '_', ' ');
end

%% Data Hmodel
errors = cell(length(experiment_names), 1);
for i = 1:length(experiment_names)
    fileID = fopen([data_path, experiment_names{i}, '.txt'], 'r');
    
    error = fscanf(fileID, '%f');
    N = length(error)/num_markers;
    error = reshape(error, num_markers, N)';
    error = error(start_offset:min(N, num_frames), :);
    
    errors{i} = mean(error, 2);
end

%% Find limits
min_error = 0; max_error = 0;
for i = 1:length(errors)
    current_max_error = mean(errors{i}) + 3 * std(errors{i});
    if (current_max_error > max_error)
        max_error = current_max_error;
    end
end

%% Plot data metric
figure_size = [0.3, 0.3, 0.3, 0.35];
figure_borders = [0.05 0.08 0.93 0.84];

if (display_time_metrics)
    figure('units', 'normalized', 'outerposition', figure_size); hold on;
    for i = 1:length(experiment_names)
        plot(1:length(errors{i}), errors{i}(:, 1), 'lineWidth', 1);
    end
    ylim([0, max_error]);
    legend(legend_names);
    if display_title, title('average data-model distance'); end
    xlabel('frame number');
    ylabel('metric');
    set(gca,'position', figure_borders, 'units','normalized');
end


%% Statistics
num_bins = 100;
thresholds = linspace(min_error, max_error, num_bins);
figure('units', 'normalized', 'outerposition', figure_size); hold on;

for i = 1:length(experiment_names)
    statistics = zeros(length(num_bins), 1);
    for j = 1:length(thresholds)
        statistics(j) = numel(find(errors{i} < thresholds(j))) / numel(errors{i});
    end
    plot(thresholds, statistics, 'lineWidth', line_width);
end
%if strcmp(dataset_type, 'tompson'), xlim([0, 30]); end
legend(legend_names, 'Location','southeast');
xlabel('error threshold');
ylabel('% frames with error < threshold');
set(gca,'fontsize', 13);
title('average markers distance');



