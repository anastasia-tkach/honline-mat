colors = {[0.4940, 0.1840, 0.5560], [0.4660, 0.6740, 0.1880], [0.3010, 0.7450, 0.9330], [0.6350, 0.0780, 0.1840], [0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290 0.6940 0.1250]};

start_offset = 10;
end_offset = 2000;
line_width = 1.5;
figure_size = [0.25, 0.25, 0.5, 0.6];

display_title = true;

num_frames = 1318;

legend_names = cell(length(experiment_names), 1);
for i = 1:length(experiment_names)
    legend_names{i} = strrep(experiment_names{i}, '_', ' ');
end

legend_type_names = cell(length(estimation_types), 1);
for i = 1:length(estimation_types)
    legend_type_names{i} = strrep(estimation_types{i}, '_', ' ');
end

%% Data Hmodel
errors = cell(length(experiment_names), 1);
sequence_length = inf;
for i = 1:length(experiment_names)
    fileID = fopen([data_path, experiment_names{i}, '.txt'], 'r');
    display([data_path, experiment_names{i}, '.txt']);
    error = fscanf(fileID, '%f');
    fclose(fileID);
    if (weighted)
        N = length(error);
        errors{i} = error;
    else
        N = length(error)/2;
        error = reshape(error, 2, N)';
        error = error(start_offset:N - end_offset, :);
        errors{i} = error(:, 1);
    end
    errors{i} = sliding_window_averaging(errors{i}, half_window_size);
    errors{i} = errors{i}(half_window_size + 1:end - half_window_size - 1, :);
    errors{i} = errors{i}(start_offset:end);
    sequence_length = min(sequence_length, length(errors{i}));
    
end

%% Find limits
min_error = inf; max_error = 0;
for i = 1:length(errors)
    current_max_error = mean(errors{i}) + 3 * std(errors{i});
    if (current_max_error > max_error)
        max_error = current_max_error;
    end
    
    current_min_error = mean(errors{i}) - 3 * std(errors{i});
    if (current_min_error < min_error)
        min_error = current_min_error;
    end
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
    ylim([min_error, max_error]);
    ylim([ylim_time_min, ylim_time_max]);
    xlabel('frame number');
    ylabel('metric');
    set(gca,'fontsize', 13);
end

%% Statistics
num_bins = 100;
thresholds = linspace(max(min_error, 0), min(max_error, 5), num_bins);
statistics = cell(length(experiment_names), 1);
for i = 1:length(experiment_names)
    statistics{i} = zeros(length(num_bins), 1);
    for j = 1:length(thresholds)
        statistics{i}(j) = numel(find(errors{i} < thresholds(j))) / numel(errors{i});
    end
end

%% Statistics limits
min_index = inf; max_index = 0;
for i = 1:length(statistics)
    current_max_index = find(statistics{i} < 0.97, 1, 'last');
    if (current_max_index > max_index)
        max_index = current_max_index;
    end
    
    current_min_index = find(statistics{i} > 0, 1, 'first');
    if (current_min_index < min_index)
        min_index = current_min_index;
    end
end


%% Legend
if (display_statistics)
    figure('units', 'normalized', 'outerposition', figure_size); hold on;
    for i = 1:length(experiment_names)
        for j = 1:length(estimation_types)
            k = strfind(experiment_names{i}, estimation_types{j});
            if (~isempty(k)), break; end;
        end
        plot(thresholds, statistics{i}, 'lineWidth', line_width, 'color', colors{1 + rem(j, length(colors))});
    end
    legend(legend_names, 'Location','southeast');
    xlim([thresholds(min_index), thresholds(max_index)]);
    %xlim([xlim_min, xlim_max]); ylim([0, 1]);
    if display_title
        xlabel('error threshold');
        ylabel('% frames with error < threshold');
        title('average data-model distance');
    end
    set(gca,'fontsize', 13);
end

%% Compute mean and variance by name
means = cell(length(estimation_types), 1);
stds = cell(length(estimation_types), 1);

for j = 1:length(estimation_types)
    M = [];
    for i = 1:length(experiment_names)
        k = strfind(experiment_names{i}, estimation_types{j});
        if (~isempty(k))
            M = [M; statistics{i}];
        end
    end
    
    means{j} = mean(M);
    stds{j} = std(M);
end

if (display_variances)
    figure('units', 'normalized', 'outerposition', figure_size); hold on;
    for i = 1:length(estimation_types)
        plot(thresholds, means{i}, 'lineWidth', 3, 'color', colors{1 + rem(i, length(colors))});
    end
    legend(legend_type_names, 'Location','southeast');
    for i = 1:length(estimation_types)
        plot(thresholds, means{i} + stds{i}, 'lineWidth', 1, 'lineStyle', '-.', 'color', colors{1 + rem(i, length(colors))});
        plot(thresholds, means{i} - stds{i}, 'lineWidth', 1, 'lineStyle', '-.', 'color', colors{1 + rem(i, length(colors))});
    end
    xlim([xlim_min, xlim_max]);
    ylim([0, 1]);
    if display_title
        xlabel('error threshold');
        ylabel('% frames with error < threshold');
        title('average data-model distance');
    end
    set(gca,'fontsize', 13);
end