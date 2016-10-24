clear; clc;
rng(36);
num_runs = 100;
Histories = cell(num_runs, 1);
for run_index = 1:num_runs
    disp(run_index);
    test_simple_problem;
    Histories{run_index} = history;
end

means = zeros(N, N);
standard_deviations = zeros(N, N);
current_run_results = zeros(num_runs, 1);
importance_means = zeros(N, N);
importance_standard_deviations = zeros(N, N);
current_run_importance = zeros(num_runs, 1);
for i = 1:N
    for j = 1:i
        for run_index = 1:num_runs
            current_run_results(run_index) = Histories{run_index}{i}.X(j);
            current_run_importance(run_index) = Histories{run_index}{i}.JtJ(j, j).^0.5;
        end
        means(i, j) = mean(current_run_results);
        standard_deviations(i, j) = std(current_run_results);
        importance_means(i, j) = mean(current_run_importance);
        importance_standard_deviations(i, j) = std(current_run_importance);
    end
end


%% Display empirical variance

frame_centrainty = T < 1.5;

display_empirical_variance(means, standard_deviations, importance_means, importance_standard_deviations, x_true, ylimit, settings, N, w2, frame_centrainty, 'simplest_problem', 1, 1);

