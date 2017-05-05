%close all;
clear;

start_offset = 700;
half_window_size = 1;
line_width = 1;
figure_size = [0.3, 0.3, 0.3, 0.35];
figure_borders = [0.05 0.08 0.93 0.84];

ylim_time_max = 5.0;
ylim_time_min = 2.3;
display_time = false;
display_statistics = true;
display_variances = false;
weighted = false;
markers = true;
estimation_types = {};

%experiment_names = {'standard_shape_0', 'standard_shape_20', 'standard_shape_40', 'extended_shape_40', 'extended_shape_20', 'extended_shape_0'};
experiment_names = {'extended_shape_40', 'extended_shape_0_full', 'top_15000', 'top_5000', 'limits', 'limits_shape_10'};

data_path = 'C:\Users\tkach\Desktop\';

legend_names = cell(length(experiment_names), 1);
for i = 1:length(experiment_names)
    legend_names{i} = strrep(experiment_names{i}, '_', ' ');
end

%% Data Hmodel
errors1 = cell(length(experiment_names), 1);
errors2 = cell(length(experiment_names), 1);
for i = 1:length(experiment_names)
    disp([data_path, experiment_names{i}, '.txt']);
    fileID = fopen([data_path, experiment_names{i}, '.txt'], 'r');
    
    error = fscanf(fileID, '%f');
    N = length(error)/2;
    error = reshape(error, 2, N)';
    error = error(start_offset:N, :);
    
    errors1{i} = error(:, 1);
    
    errors2{i} = error(:, 2);
    
    if mean(errors2{i}) < 0.1
        errors2{i} = errors2{i} * 1000;
    end
    
    errors1{i} = sliding_window_averaging(errors1{i}, half_window_size);
    errors2{i} = sliding_window_averaging(errors2{i}, half_window_size);
    
    errors1{i} = errors1{i}(half_window_size + 1:end - half_window_size - 1, :);
    errors2{i} = errors2{i}(half_window_size + 1:end - half_window_size - 1, :);
end


%%{
%% Plot data metric
figure('units', 'normalized', 'outerposition', figure_size); hold on;
for i = 1:length(experiment_names)
    plot(1:length(errors1{i}), errors1{i}(:, 1), 'lineWidth', 1);
    %plot(1:1000, errors1{i}(end - 1000 + 1:end, 1), 'lineWidth', 1);
end
legend(legend_names);

xlabel('frame number');
ylabel('metric');
set(gca,'position', figure_borders, 'units','normalized');


%% Statistics

num_bins = 100;
min_error = 3.3;
max_error = 5.1;
errors = errors1;
thresholds = linspace(min_error, max_error, num_bins);
figure('units', 'normalized', 'outerposition', figure_size); hold on;

for i = 1:length(experiment_names)
    statistics = zeros(length(num_bins), 1);
    for j = 1:length(thresholds)
        statistics(j) = numel(find(errors{i} < thresholds(j))) / numel(errors{i});
    end
    plot(thresholds, statistics, 'lineWidth', line_width);
end
xlim([min_error, max_error]);
legend(legend_names, 'Location','southeast');

set(gca,'position', figure_borders, 'units','normalized');
