%{
clear; clc; close all;
N = 70;
B = 3;
beta_true = [3; 3; 3];
num_samples = 15;
measurement_noise_std = 0.07;
theta_noise_std = 0.15;
experiment_name = 'batch2';
dataset_test = false;

num_items = 100;
beta_results = zeros(num_items, B - 1);
for item = 1:num_items
    disp(item);
    dataset_entry_number = item;
    p = rand; 
    if p < 0.1
        dataset_test = true;
    else
        dataset_test = false;
    end
    %dataset_test = true;
    main_batch_with_latent;
    beta_results(item, :) = history{item}.mean(1:B - 1);
end

save(['saved_variables\', experiment_name], 'beta_results');
return
%}
%% Load and display
num_items = 70;
experiment_names = {'modified_measurement_noise_01', 'batch2'};%, 'no_filtering_damping_200', 'no_filtering_damping_20'};
%experiment_names = {'IEKF', 'modified_update_with_max', 'modified_update_no_pose_certainty', 'modified_measurement_noise_01', 'modified_measurement_noise_same'};
results = cell(length(experiment_names), 1);
mean_values = zeros(length(experiment_names), B - 1);
for i = 1:length(experiment_names)
    load(['saved_variables\', experiment_names{i}]);
    results{i} = beta_results;
    results{i} = results{i}(1:num_items, :);
    mean_values(i, 1) = mean(abs(beta_results(:, 1) - beta_true(1)));
    mean_values(i, 2) = mean(abs(beta_results(:, 2) - beta_true(2)));
end

figure('units', 'normalized', 'outerposition', [0.1, 0.2, 0.33, 0.7]); hold on;
for b = 1:B - 1
    h = subplot(B - 1, 1, b); hold on;
    p = get(h, 'pos'); set(h, 'pos', [0.07, p(2) - 0.06, 0.9, 0.43]);
    set(gca,'XTick',[]);
    
    bar_handle = bar(1:length(experiment_names), mean_values(:, b));
    set(bar_handle,'FaceColor', [196, 229, 230]/255, 'EdgeColor', 'none', 'lineWidth', 2);
    
    for i = 1:length(experiment_names)
        x_values = i + 0.05 * randn(num_items, 1);
        y_values = results{i}(:, b);
        scatter(x_values, abs(y_values - beta_true(b)), 20, [3, 155, 199]./255, 'filled');
    end
    
    %ylim([0, 0.7]);
    set(gca, 'fontSize', 12); title(['beta ', num2str(b)]);
end
