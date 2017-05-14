%close all;
clear;

start_offset = 1;
half_window_size = 1;



%% Real data - teaser
folder_name = 'teaser_FINAL'; sequence_names = {'teaser'};
estimation_types = {'SHARP', 'HTRACK', 'TAYLOR', ...
    'TEMPLATE_0.200000', 'KALMAN_STANDARD_EVALUATION_0.200000',  'KALMAN_EXTENDED_EVALUATION_0.200000', 'ONLINE_CALIBRATION_0.200000'};

display_time = false;
display_statistics = true;
xlim_min = 3.0;
xlim_max = 8.0;

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
            experiment_names{end + 1} = [estimation_types{j}, '_', sequence_names{s}, '_', series_names{i}];
        end
    end
end

data_path = ['E:\Data\honline-results\', folder_name, '\'];
legend_names = cell(length(experiment_names), 1);
for i = 1:length(experiment_names), legend_names{i} = strrep(experiment_names{i}, '_', ' '); end
legend_type_names = cell(length(estimation_types), 1);
for i = 1:length(estimation_types), legend_type_names{i} = strrep(estimation_types{i}, '_', ' '); end


%% Compute online metrics
line_width = 1.5;
figure_size = [0.25, 0.25, 0.5, 0.6];
display_title = true;

legend_names = cell(length(experiment_names), 1);
for i = 1:length(experiment_names)
    legend_names{i} = strrep(experiment_names{i}, '_', ' ');
end

%% Data Hmodel
errors = cell(length(experiment_names), 1);
sequence_length = inf;
for i = 1:length(experiment_names)
    fileID = fopen([data_path, experiment_names{i}, '.txt'], 'r');
    display([data_path, experiment_names{i}, '.txt']);
    error = fscanf(fileID, '%f');
    fclose(fileID);
    N = length(error)/2;
    error = reshape(error, 2, N)';
    error = error(start_offset:N, :);
    errors{i} = error(:, 1);
    sequence_length = min(sequence_length, length(errors{i}));
end

%% Time
if (display_time)
    figure('units', 'normalized', 'outerposition', figure_size); hold on;
    for i = 1:length(experiment_names)
        for j = 1:length(estimation_types)
            k = strfind(experiment_names{i}, estimation_types{j});
            if (~isempty(k)), break; end;
        end
        plot(1:sequence_length - start_offset, errors{i}(1:sequence_length - start_offset, 1), 'lineWidth', 2, 'color', colors{1 + rem(j, length(colors))});
    end
    legend(legend_names);
    if display_title, title('average data-model distance'); end
    xlim([1, sequence_length - start_offset]);
    ylim([ylim_time_min, ylim_time_max]);
    xlabel('frame number');
    ylabel('metric');
    set(gca,'fontsize', 13);
end

%% Statistics
num_bins = 150;
thresholds = linspace(xlim_min, xlim_max, num_bins);
statistics = cell(length(experiment_names), 1);
for i = 1:length(experiment_names)
    statistics{i} = zeros(length(num_bins), 1);
    for j = 1:length(thresholds)
        statistics{i}(j) = numel(find(errors{i} < thresholds(j))) / numel(errors{i});
    end
end

%% Legend
dark_red = [178, 68, 117]/255;
light_red = [217, 154, 143]/255;
dark_green = [61, 131, 119]/255;
light_green = [144, 194, 171]/230;
grey = [0.7, 0.7, 0.7];
black = [0.3, 0.3, 0.3];
colors = {[0.3, 0.3, 0.3],  light_green, dark_green, [0.4940, 0.2840, 0.5560], light_red , dark_red, grey };

if (display_statistics)
    figure('units', 'normalized', 'outerposition', figure_size); hold on;
    for i = 1:length(experiment_names)
        plot(thresholds, statistics{i}, 'lineWidth', 1.5, 'color', colors{i});
    end
    %legend(legend_names);
    legend({'Sharp', 'Htrack', 'Taylor', 'MVS template', 'Kalman Standard',  'Kalman Extended', 'Online'});
    xlim([xlim_min, xlim_max]); ylim([0, 1]);
    if display_title
        xlabel('error threshold');
        ylabel('% frames with error < threshold');
        title('average data-model distance');
    end
    set(gca,'fontsize', 15, 'fontname', 'Cambria');
    box on; set(gca,'linewidth', 1.5);
    set(gca, 'xcolor', [0.2, 0.2, 0.2]);  set(gca, 'ycolor', [0.2, 0.2, 0.2]);
end


