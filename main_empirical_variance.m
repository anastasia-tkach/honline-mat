clear; clc;
rng default;
num_runs = 100;
Histories = cell(num_runs, 1);
for run_index = 1:num_runs
    disp(run_index);
    main_batch;
    Histories{run_index} = history;
end

beta_index = 2;

means = zeros(N, N);
standard_deviations = zeros(N, N);
current_run_results = zeros(num_runs, 1);
importance_means = zeros(N, N);
importance_standard_deviations = zeros(N, N);
current_run_importance = zeros(num_runs, 1);
for i = 1:N
    for j = 1:i
        for run_index = 1:num_runs
            current_run_results(run_index) = Histories{run_index}{i}.betas{j}(beta_index);
            current_run_importance(run_index) = Histories{run_index}{i}.JtJ((B + T) * (j - 1) + beta_index, (B + T) * (j - 1) + beta_index).^0.5;
        end
        means(i, j) = mean(current_run_results);
        standard_deviations(i, j) = std(current_run_results);
        importance_means(i, j) = mean(current_run_importance);
        importance_standard_deviations(i, j) = std(current_run_importance);
    end
end


%% Display empirical variance

ylimit = [1.5, 4.5];
frame_centrainty = (thetas_true(:, 2) >= pi/4) .* (thetas_true(:, 3) >= pi/4);

display_empirical_variance(means, standard_deviations, importance_means, importance_standard_deviations, beta_true(beta_index), ylimit, settings, N, w2, frame_centrainty, 'sticks_finger');
