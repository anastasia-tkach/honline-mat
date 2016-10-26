%%{
clear; clc;
rng default;
num_runs = 100;
Histories = cell(num_runs, 1);
for run_index = 1:num_runs
    disp(run_index);
    main_batch;
    Histories{run_index} = history;
end
%%}5
beta_indices = 1:2;

%% Compute statistics

for beta_index = beta_indices
    means = zeros(N, N);
    standard_deviations = zeros(N, N);
    current_run_results = zeros(num_runs, 1);
    importance_means = zeros(N, N);
    importance_standard_deviations = zeros(N, N);
    current_run_importance = zeros(num_runs, 1);
    
    for i = 1:N
        for j = 1:i
            for run_index = 1:num_runs
                if j > i - settings.batch_size
                    indices = (B + T) * (settings.batch_size - i + j - 1) + beta_index;
                    current_run_results(run_index) = Histories{run_index}.x_batch(i, indices);
                    current_run_importance(run_index) = Histories{run_index}.h_batch(i, indices).^0.5;
                else
                    indices = j + settings.batch_size - 1;
                    current_run_results(run_index) = Histories{run_index}.x_batch(indices, beta_index);
                    current_run_importance(run_index) = Histories{run_index}.h_batch(indices, beta_index).^0.5;
                end
            end
            means(i, j) = mean(current_run_results);
            standard_deviations(i, j) = std(current_run_results);
            importance_means(i, j) = mean(current_run_importance);
            importance_standard_deviations(i, j) = std(current_run_importance);
        end
    end
    
    %% Set display parameters
    ylimit = [1.5, 4.5];
    if length(beta_indices) == 1
        frame_centrainty = [(thetas_true(:, 2) >= pi/4) .* (thetas_true(:, 3) >= pi/4), (thetas_true(:, 2) >= pi/4) .* (thetas_true(:, 3) >= pi/4)];
    else
        frame_centrainty = [thetas_true(:, 2) >= pi/4 , thetas_true(:, 3) >= pi/4];
    end
    
    %% Display empirical variance
    display_empirical_variance(means, standard_deviations, importance_means, importance_standard_deviations, beta_true(beta_index), ylimit, settings, N, frame_centrainty, 'sticks_finger', beta_indices, beta_index);
    
    %% Display history
    %display_history_with_variance(means, standard_deviations, importance_means, importance_standard_deviations, beta_true(beta_index), ylimit, settings, N, w2, frame_centrainty, 'sticks_finger', beta_indices, beta_index);
    
end