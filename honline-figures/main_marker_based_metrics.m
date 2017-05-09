line_width = 1.5;
display_title = true;
num_markers = 36;

dataset_type = 'tkach';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data_path = 'E:/Data/sensor-sequences/sridhar/';
% experiment_names = {'marker_based_metrics_nof','marker_based_metrics'};
% % data_path = 'E:\OneDrive\EPFL\Code\honline-mat\honline-figures\htrack_results_on_dexter\';
% % experiment_names = {'errors_all'};
% num_markers = 6;
% xlim_max = 30;
% xlim_min = 0;
% 
% subsequence_limits = [1, 425, 940, 1280, 1672, 2113, 2647, 3154];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data Hmodel
errors = cell(length(experiment_names), 1);

%% Markers analysis
%{
for k = 1:num_markers / 7
    errors = cell(length(experiment_names), 1);
    start = 7 * (k - 1) + 1;
    finish = 7 * k;
    for i = start:finish
        display([data_path, experiment_names{1}, '.txt']);
        fileID = fopen([data_path, experiment_names{1} , '.txt'], 'r');
        
        error = fscanf(fileID, '%f');
        N = length(error)/num_markers;
        error = reshape(error, num_markers, N)';
        error = error(start_offset:N, :);
        
        errors{i} = error(:, i); legend_names{i - start + 1} = num2str(i);
    end
%}


%% Experiments

for i = 1:length(experiment_names)
    display([data_path, experiment_names{i}, '.txt']);
    fileID = fopen([data_path, experiment_names{i} , '.txt'], 'r');
    
    error = fscanf(fileID, '%f');
    N = length(error)/num_markers;
    error = reshape(error, num_markers, N)';
    error = error(start_offset:N, :);    
   
    %errors{i} = mean(error, 2);
    %errors{i} = max(error, [], 2);
    errors{i} = max(error(:, [1:27]), [], 2);
    %errors{i} = mean(error(:, [1:27]), 2);
    
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
figure_size = [0.25, 0.25, 0.5, 0.6];
figure_borders = [0.05 0.08 0.93 0.84];

if (display_time)
    figure('units', 'normalized', 'outerposition', figure_size); hold on;
    for i = 1:length(experiment_names)
        plot((1:length(errors{i})), errors{i}(:, 1), 'lineWidth', 1);
    end
    ylim([0, 60]);
    legend(legend_names);
    if display_title, title('average data-model distance'); end
    xlabel('frame number');
    ylabel('metric');
    set(gca,'position', figure_borders, 'units','normalized');
end

%% Statistics
num_bins = 100;
thresholds = linspace(min_error, min(xlim_max, max_error), num_bins);
thresholds = linspace(0,30, 200);
figure('units', 'normalized', 'outerposition', figure_size); hold on;

for i = 1:length(errors)
    statistics = zeros(length(num_bins), 1);
    for j = 1:length(thresholds)
        statistics(j) = numel(find(errors{i} < thresholds(j))) / numel(errors{i});
    end
    plot(thresholds, statistics, 'lineWidth', 3);
end
%if strcmp(dataset_type, 'tompson'), xlim([0, 30]); end
legend(legend_names, 'Location','southeast');
xlabel('error threshold');
ylabel('% frames with error < threshold');
set(gca,'fontsize', 13);
title('average markers distance');
xlim([xlim_min, xlim_max]); ylim([0, 1]);
%end
return

%% Split and reweight
figure('units', 'normalized', 'outerposition', figure_size); hold on;
subsequence_statistics = cell(length(subsequence_limits) - 1, 1);
for i = 1:length(experiment_names)
    for k = 1:length(subsequence_limits) - 1
        subsequence_errors = errors{i}(subsequence_limits(k):subsequence_limits(k + 1));
        subsequence_statistics{k} = zeros(length(num_bins), 1);
        for j = 1:length(thresholds)
            subsequence_statistics{k}(j) = numel(find(subsequence_errors < thresholds(j))) / numel(subsequence_errors);
        end
        plot(thresholds,  subsequence_statistics{k}, 'lineWidth', line_width);
    end
end


weights = zeros(length(subsequence_limits) - 1, 1);
for k = 1:length(weights)
    weights(k) = (subsequence_limits(k + 1) - subsequence_limits(k)) / (subsequence_limits(end) - subsequence_limits(1));
end
weights = 1/7 * ones(7, 1);

weighted_statistics = zeros(length(thresholds), 1);
for i = 1:length(thresholds)
    for k = 1:length(subsequence_statistics)
        weighted_statistics(i) =  weighted_statistics(i) + weights(k) * subsequence_statistics{k}(i);
    end
end
plot(thresholds,  weighted_statistics, 'lineWidth', 3);
set(gca,'fontsize', 13);
legend({'1', '2', '3', '4', '5', '6', '7'}, 'Location','southeast');
xlim([xlim_min, xlim_max]); ylim([0, 1]);