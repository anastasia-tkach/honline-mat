beta_true = [37.1409 28.436 14.0033 37.0552 20.5967 12.8329 40.4944 23.2687 15.9074 37.9263 23.9032 13.6009 31.9997 19.1319 12.8205 9.72575 3.87206 ...
    -7.16046 25.3963 50.8191 5.28707 8.13206 52.8925 10.5463 -5.93253 49.1046 10.7991 -18.1006 44.8872 6.8978 23.9371 45.629 -17.9329 38.572 ...
    7.31059 3.0928 -10.0083 1.01374 28.5967 41.231 4.05823 48.6988 12 6.73474 1.52305 15.6827 10.3118 7.30846 7.0311 8.44143 7.55251 ...
 5.85299 5.17427 7.68834 7.64302 5.67621 5.58432 7.26768 6.97092 5.01217 4.84959 7.84562 6.22559 4.77166 4.18002 9.50548 10.726 10.2172 8.98482 13.2994 13.6244 13.4193];

num_betas = 72;
num_thetas = 34;

%% Data Hmodel
errors = cell(length(experiment_names), 1);

for i = 1:length(experiment_names)
    display([data_path, experiment_names{i}, '.txt']);
    fileID = fopen([data_path, experiment_names{i} , '.txt'], 'r');
    
    thetas_betas = fscanf(fileID, '%f');
    N = length(thetas_betas)/(num_betas + num_thetas);
    thetas_betas = reshape(thetas_betas, num_betas + num_thetas, N)';
    betas = thetas_betas(start_offset:N, num_thetas + 1:end);
    errors{i} = abs(betas - repmat(beta_true, length(betas), 1));
    errors{i} = mean(errors{i}, 2);
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
figure('units', 'normalized', 'outerposition', figure_size); hold on;

for i = 1:length(errors)
    statistics = zeros(length(num_bins), 1);
    for j = 1:length(thresholds)
        statistics(j) = numel(find(errors{i} < thresholds(j))) / numel(errors{i});
    end
    plot(thresholds, statistics, 'lineWidth', 3);
end

legend(legend_names, 'Location','southeast');
xlabel('error threshold');
ylabel('% frames with error < threshold');
set(gca,'fontsize', 13);
title('average markers distance');
xlim([xlim_min, xlim_max]); ylim([0, 1]);


