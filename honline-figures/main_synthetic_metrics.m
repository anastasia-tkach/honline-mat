beta_true = [37.1409 28.436 14.0033 37.0552 20.5967 12.8329 40.4944 23.2687 15.9074 37.9263 23.9032 13.6009 31.9997 19.1319 12.8205 9.72575 3.87206 -7.16046 ...
    25.3963 50.8191 2.28707 8.13206 52.8925 7.5463 -5.93253 49.1046 7.7991 -18.1006 44.8872 3.8978 23.9371 45.629 -17.9329 38.572 7.31059 3.0928 ...
    -10.0083 1.01374 28.5967 41.231 4.05823 48.6988 12 6.73474 1.52305 15.6827 10.3118 7.30846 7.0311 8.44143 7.55251 5.85299 5.17427 7.68834 7.64302 ...
    5.67621 5.58432 7.26768 6.97092 5.01217 4.84959 7.84562 6.22559 4.77166 4.18002 9.50548 10.726 10.2172 8.98482 13.2994 13.6244 13.4193];

num_betas = 72;
num_thetas = 0;
num_iters = 10;


%% Data Hmodel
sequence_length = inf;
errors = cell(length(experiment_names), 1);
errors_norm = cell(length(experiment_names), 1);

for i = 1:length(experiment_names)
    display([data_path, experiment_names{i}, '.txt']);
    fileID = fopen([data_path, experiment_names{i} , '.txt'], 'r');
    
    thetas_betas = fscanf(fileID, '%f');
    N = length(thetas_betas)/(num_betas + num_thetas);
    thetas_betas = reshape(thetas_betas, num_betas + num_thetas, N)';
    betas = thetas_betas(start_offset:N, num_thetas + 1:end);
    
    errors{i} = betas - repmat(beta_true, length(betas), 1);
    errors_norm{i} = zeros(size(errors{i}, 1), 1);
    for j = 1:length(errors{i})
        errors_norm{i}(j) = norm(errors{i}(j, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 45:71] + 1));
    end
    fclose(fileID);
    
    errors{i} = errors{i}(1:num_iters:end, :);
    errors_norm{i} = errors_norm{i}(1:num_iters:end, :);
    sequence_length = min(sequence_length, length(errors_norm{i}));
end;

for i = 1:length(errors_norm), errors_norm{i} = errors_norm{i}(1:sequence_length); end

%% Plot data metric
figure_size = [0.25, 0.25, 0.5, 0.6];
figure_borders = [0.05 0.08 0.93 0.84];

if (display_time)
    figure('units', 'normalized', 'outerposition', figure_size); hold on;
    %plot(1:length(errors{i}), zeros(length(errors{i}), 1), 'lineWidth', 1, 'color', 'k');
    for i = 1:length(experiment_names)
        %for j = [0, 1, 2, 3, 4, 5, 19, 6, 7, 8, 22, 9, 10, 11, 25, 12, 13, 14, 28] + 1
        %for j = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 45:71] + 1
        %   plot((1:length(errors{i})), errors{i}(:, j), 'lineWidth', 2);
        %end
        plot((1:length(errors_norm{i})), errors_norm{i}, 'lineWidth', 3);
    end
    legend(legend_names);
    xlabel('frame number');
    ylabel('metric');
    set(gca,'position', figure_borders, 'units','normalized');
    xlim([xlim_min, xlim_max]); ylim([0, 30]);
end

%% Compute mean and variance by name
means = cell(length(estimation_types), 1);
stds = cell(length(estimation_types), 1);

for j = 1:length(estimation_types)
    M = [];
    for i = 1:length(experiment_names)
        k = strfind(experiment_names{i}, estimation_types{j});
        if (~isempty(k)), M = [M, errors_norm{i}]; end
    end
    means{j} = mean(M, 2);
    stds{j} = std(M, [], 2);
end

figure('units', 'normalized', 'outerposition', figure_size); hold on;
for i = 1:length(estimation_types), 
    plot((1:length(means{i})), means{i}, 'lineWidth', 3, 'color', colors{rem(i, 8) + 1}); 
end
legend(legend_type_names);
for i = 1:length(estimation_types)
    plot((1:length(means{i})), means{i}, 'lineWidth', 3, 'color', colors{rem(i, 8) + 1});
    plot((1:length(means{i})), means{i} - stds{i}, 'lineWidth', 2, 'color', colors{rem(i, 8) + 1}, 'lineStyle', '-.');
    plot((1:length(means{i})), means{i} + stds{i}, 'lineWidth', 2, 'color', colors{rem(i, 8) + 1}, 'lineStyle', '-.');
end

xlabel('frame number');
ylabel('metric');
set(gca,'position', figure_borders, 'units','normalized');
xlim([xlim_min, xlim_max]); ylim([0, 30]);

%% Initail pertrubation vs resuls
dark_red = [188, 58, 117]/255;
light_red = 0.8 * [217, 154, 143]/255 + 0.2 * [238, 198, 199]/255;
dark_green = [33, 159, 126]/255;
light_green = [144, 194, 171]/230;
grey = [0.75, 0.75, 0.75];
colors = {dark_red, light_red, grey, light_green, dark_green};

num_steps = 5;

interval_means = zeros(length(estimation_types), num_steps);
interval_stds = zeros(length(estimation_types), num_steps);
pertrubations = [0.05, 0.075, 0.10, 0.125, 0.15, 0.175, 0.20, 0.225, 0.25, 0.275, 0.30, 0.325, 0.35, 0.375, 0.40];
%intervals = round(linspace(1, sequence_length, num_steps));
intervals = [1, 10, 100, 300, 850];

for i = 1:length(estimation_types)
    for j = 1:num_steps
        interval_means(i, j) = means{i}(intervals(j)); 
        interval_stds(i, j) = stds{i}(intervals(j));
    end
end
figure('units', 'normalized', 'outerposition', figure_size); hold on;
for i = 1:num_steps
    plot(pertrubations, interval_means(:, i), 'lineWidth', 3, 'color', colors{i}); 
end
legend({'frame 0', 'frame 10', 'frame 100', 'frame 300', 'frame 850'}); legend('boxoff');

for i = 1:num_steps   
    plot(pertrubations, interval_means(:, i) - interval_stds(:, i), 'lineWidth', 1, 'lineStyle', '-.', 'color', colors{i});
    plot(pertrubations, interval_means(:, i) + interval_stds(:, i), 'lineWidth', 1, 'lineStyle', '-.', 'color', colors{i});
    plot(pertrubations, interval_means(:, i), 'lineWidth', 3, 'color', colors{i}); 
end
xlim([pertrubations(1), pertrubations(end)]); ylim([0, 40]);
box on; set(gca,'linewidth', 1.5); set(gca, 'fontsize', 14, 'fontname', 'Cambria');
set(gca, 'xcolor', [0.2, 0.2, 0.2]);
set(gca, 'ycolor', [0.2, 0.2, 0.2]); 
xlabel('template pertrubation \sigma');
ylabel('||\beta - \beta{true}||, mm');
return

%% Plot time
errors = cell(length(experiment_names), 1);
sequence_length = inf;
for i = 1:length(experiment_names)
    filename = [data_path, experiment_names{i}(1:end-10), '_weighted', '.txt']; weighted = true;
    %filename = [data_path, experiment_names{i}(1:end-10),'.txt'];
    fileID = fopen(filename, 'r');
    display(filename);
    error = fscanf(fileID, '%f');
    fclose(fileID);
    if (weighted)
        N = length(error);
        errors{i} = error;
    else
        N = length(error)/2;
        error = reshape(error, 2, N)';
        error = error(start_offset:N, :);
        errors{i} = error(:, 1);
    end
    
    errors{i} = errors{i}(num_iters:num_iters:end, :);
    
    errors{i} = sliding_window_averaging(errors{i}, half_window_size);
    errors{i} = errors{i}(half_window_size + 1:end - half_window_size - 1, :);
    errors{i} = errors{i}(start_offset:end);
    
    sequence_length = min(sequence_length, length(errors{i}));
end


figure('units', 'normalized', 'outerposition', figure_size); hold on;
for i = 1:length(experiment_names)
    for j = 1:length(estimation_types)
        k = strfind(experiment_names{i}, estimation_types{j});
        if (~isempty(k)), break; end;
    end
    plot(1:sequence_length - start_offset, errors{i}(1:sequence_length - start_offset, 1), 'lineWidth', 2, 'color', colors{1 + rem(j, length(colors))});
end
legend(legend_names);
xlim([1, sequence_length - start_offset]);
ylim([0, 5]);
xlabel('frame number');
ylabel('metric');
set(gca,'fontsize', 13);



