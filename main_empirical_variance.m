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
%%}
beta_indices = 1:3;

%% Compute statistics

for beta_index = beta_indices
    means = zeros(settings.num_frames, settings.num_frames);
    standard_deviations = zeros(settings.num_frames, settings.num_frames);
    current_ij_results = zeros(num_runs, 1);
    importance_means = zeros(settings.num_frames, settings.num_frames);
    importance_standard_deviations = zeros(settings.num_frames, settings.num_frames);
    current_ij_importance = zeros(num_runs, 1);
    
    results_history = zeros(num_runs, settings.num_frames, B);
    if settings.store_covariance
        covariance_history = zeros(num_runs, settings.num_frames, B, B);
    end
    
    for i = 1:N
        for j = 1:i
            for run_index = 1:num_runs
                if j > i - settings.batch_size
                    indices = (B + T) * (settings.batch_size - i + j - 1) + beta_index;
                    current_ij_results(run_index) = Histories{run_index}.x_batch(i, indices);
                    current_ij_importance(run_index) = Histories{run_index}.h_batch(i, indices).^0.5;
                else
                    indices = j + settings.batch_size - 1;
                    current_ij_results(run_index) = Histories{run_index}.x_batch(indices, beta_index);
                    current_ij_importance(run_index) = Histories{run_index}.h_batch(indices, beta_index).^0.5;
                end
                
                %% prepare for empirical covariance
                if (i == j)
                    results_history(run_index, i, :) = Histories{run_index}.x_batch(i, end - B - T + 1:end - T);
                    if settings.store_covariance
                        covariance_history(run_index, i, :, :) = Histories{run_index}.covariance(i, :, :);
                    end
                end
                
            end
            means(i, j) = mean(current_ij_results);
            standard_deviations(i, j) = std(current_ij_results);
            importance_means(i, j) = mean(current_ij_importance);
            importance_standard_deviations(i, j) = std(current_ij_importance);
        end
    end
    
    %% Set display parameters
    ylimit = [1.5, 4.5];
    if length(beta_indices) == 1
        v = [(thetas_true(:, 2) >= pi/4) .* (thetas_true(:, 3) >= pi/4), (thetas_true(:, 2) >= pi/4) .* (thetas_true(:, 3) >= pi/4)];
    end
    if length(beta_indices) == 2
        frame_certainty = [thetas_true(:, 2) >= pi/4 , thetas_true(:, 3) >= pi/4];
    end
    if length(beta_indices) == 3
        frame_certainty = [thetas_true(:, 2) >= pi/4 , (thetas_true(:, 2) >= pi/4) .* (thetas_true(:, 3) >= pi/4), thetas_true(:, 3) >= pi/4];
    end
    
    %% Display empirical variance
    display_empirical_variance(means, standard_deviations, importance_means, importance_standard_deviations, beta_true(beta_index), ylimit, settings, settings.num_frames, frame_certainty, 'sticks_finger', beta_indices, beta_index);
    
    %% Display history
    %display_history_with_variance(means, standard_deviations, importance_means, importance_standard_deviations, beta_true(beta_index), ylimit, settings, settings.num_frames, frame_certainty, 'sticks_finger', beta_indices, beta_index);
    
end

%% Display covariance 1-2
settings.num_runs = num_runs;
%display_covariance(settings, results_history, covariance_history, frame_certainty);
