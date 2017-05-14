%close all;
clear;

start_offset = 1;
half_window_size = 1;

weighted = false;
full = false;
markers = false;
synthetic = false;
display_time = false;


%% Real data - Tompson
folder_name = 'tompson_FINAL'; sequence_names = {'tompson'};
estimation_types = {'TEMPLATE_0.000000', 'KALMAN_STANDARD_EVALUATION_0.000000','BATCH_OFFLINE_EVALUATION_0.000000', ...
    'KALMAN_EXTENDED_EVALUATION_0.000000', 'ONLINE_CALIBRATION_0.000000', ...
    'BATCH_ONLINE_EVALUATION_0.000000', };
markers = true;


display_statistics = true;

xlim_max = 30;
xlim_min = 0;

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
            experiment_names{end + 1} = [estimation_types{j}, '_', sequence_names{s}, '_', series_names{i}, '_markers'];
        end
    end
end

data_path = ['E:\Data\honline-results\', folder_name, '\'];
legend_names = cell(length(experiment_names), 1);
for i = 1:length(experiment_names), legend_names{i} = strrep(experiment_names{i}, '_', ' '); end
legend_type_names = cell(length(estimation_types), 1);
for i = 1:length(estimation_types), legend_type_names{i} = strrep(estimation_types{i}, '_', ' '); end

%% Compute marker based metrics
display_title = true;
num_markers = 36;
errors = cell(length(experiment_names), 1);

%% Experiments
for i = 1:length(experiment_names)
    display([data_path, experiment_names{i}, '.txt']);
    fileID = fopen([data_path, experiment_names{i} , '.txt'], 'r');
    
    error = fscanf(fileID, '%f');
    N = length(error)/num_markers;
    error = reshape(error, num_markers, N)';
    error = error(start_offset:N, :);
    
    %errors{i} = max(error(:, [1:27]), [], 2);
    errors{i} = mean(error(:, [1:27]), 2);
end


%% Plot
figure_size = [0.25, 0.25, 0.5, 0.6];
figure_borders = [0.05 0.08 0.93 0.84];

dark_red = [178, 68, 117]/255;
light_red = [217, 154, 143]/255;
dark_green = [61, 131, 119]/255;
light_green = [144, 194, 171]/230;
grey = [0.7, 0.7, 0.7];
black = [0.3, 0.3, 0.3];
colors = {black, light_red,  grey, dark_red, dark_green, light_green, };

num_bins = 200;
thresholds = linspace(0,30, num_bins);
figure('units', 'normalized', 'outerposition', figure_size); hold on;

for i = 1:length(errors)
    statistics = zeros(length(num_bins), 1);
    for j = 1:length(thresholds)
        statistics(j) = numel(find(errors{i} < thresholds(j))) / numel(errors{i});
    end
    plot(thresholds, statistics, 'lineWidth', 1.5, 'color', colors{i});
end
legend( {'Scaled Template', 'Kalman Standard','Batch Offline', 'Kalman Extended', 'Online', 'Batch Online' });
xlabel('marker error threshold \epsilon');
ylabel('% frames with max marker distance < \epsilon');
set(gca,'fontsize', 13);
xlim([xlim_min, xlim_max]); ylim([0, 1]);
set(gca,'fontsize', 15, 'fontname', 'Cambria');
box on; set(gca,'linewidth', 1.5);
set(gca, 'xcolor', [0.2, 0.2, 0.2]);  set(gca, 'ycolor', [0.2, 0.2, 0.2]);


