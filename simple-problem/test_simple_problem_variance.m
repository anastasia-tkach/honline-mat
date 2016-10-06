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

display_empirical_variance(means, standard_deviations, importance_means, importance_standard_deviations, x_true, ylimit, settings, N, w2, frame_centrainty);

%% Display last iterations
%{
figure('units', 'normalized', 'outerposition', [0.1, 0.3, 0.45, 0.45]); hold on;
set(gca,'position', [0.05 0.08 0.93 0.90], 'units','normalized');
for j = 1:N
    if (T(j) > 1.5)
        rectangle('Position',[j, ylimit(1), 1, ylimit(2) - ylimit(1)],'FaceColor',[1; 0.96; 0.93],'EdgeColor','none')
    else
        rectangle('Position',[j, ylimit(1), 1, ylimit(2) - ylimit(1)],'FaceColor',[0.96; 1; 0.90],'EdgeColor','none')
    end
end
plot(0:length(history), x_true * ones(length(history) + 1, 1), 'lineWidth', 2, 'color', [1, 0.7, 0.3], 'lineStyle', '-.');
plot(1:N, means(N, :), '.-', 'lineWidth', 2, 'markersize', 13, 'color', [1, 0.5, 0.4]);
plot(1:N, means(N, :) + standard_deviations(N, :), 'lineWidth', 2, 'color', [0.65, 0.8, 0.6]);
plot(1:N, means(N, :) - standard_deviations(N, :), 'lineWidth', 2, 'color', [0.65, 0.8, 0.6]);

set(gca, 'fontSize', 13); xlim([1, N]); ylim(ylimit);
return;
%}
