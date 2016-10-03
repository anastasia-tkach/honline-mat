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


%% Display last iteration
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

%% Display
w = N * 0.040909;
offset = 1/N;
line_color = [1, 0.85, 0.5];
point_color = [1.0 0.5 0.3]; %[0.3, 0.6, 0.8];
f = figure('units', 'normalized', 'outerposition', [0.1, 0.3, w, 0.55]); hold on;
set(gca,'position', [0.06 0.06 0.87 0.85], 'units','normalized');
for j = 1:N
    if (T(j) > 1.5)
        rectangle('Position',[j, ylimit(1), 1, ylimit(2) - ylimit(1)],'FaceColor',[1; 0.98; 0.95],'EdgeColor','none')
    else
        rectangle('Position',[j, ylimit(1), 1, ylimit(2) - ylimit(1)],'FaceColor',[0.96; 1; 0.93],'EdgeColor','none')
    end
end
plot(1:length(history) + 1, x_true * ones(length(history) + 1, 1), 'lineWidth', 2, 'color', [1, 0.5, 0.4], 'lineStyle', '-.');


max_value = max(importance_means(:) + importance_standard_deviations(:));
%importance_means = importance_means/max_value * 0.4;
%importance_standard_deviations = importance_standard_deviations/max_value * 0.4;

for j = 1:N
    
    for k = 1:j
        current_line_color = (1 - k / j ) * [0.95, 0.88, 0.88] + k / j * line_color;
        current_point_color = (1 - k / j ) * [0.85, 0.75, 0.75] + k / j * point_color;
        point_size = 20;
        line_width = 3.2;
        if k == j
            current_line_color = [1, 0.75, 0.3];
            current_point_color = [1.0 0.45 0.3];%[0.3, 0.7, 0.8];
            point_size = 40;
            line_width = 4.5;
        end
        
        yyaxis left; 
        myline([j + offset * k, means(j, k) + standard_deviations(j, k)], ...
            [j + offset * k, means(j, k) - standard_deviations(j, k)], current_line_color, line_width);
        mypoint([j + offset * k, means(j, k)], current_point_color, point_size);         
        
        if k == 1, myline([j, ylimit(1)], [j, ylimit(2)], [0.88, 0.88, 0.88], 2); end
        
        % importance  
        yyaxis right; 
        y_position = 0;
        myline([j + offset * k, y_position], ...
            [j + offset * k, y_position + importance_means(j, k) + importance_standard_deviations(j, k)], [0.75, 0.9, 0.7], 3.2);
        myline([j + offset * k, y_position], ...
            [j + offset * k, y_position + importance_means(j, k)], [0.65, 0.8, 0.6], 3.2);  
        
    end
end

set(gca, 'fontSize', 14); set(gca,'fontname','Cambria');
algorithm_name = '';
if quadratic_two == true, algorithm_name = 'quadratic-two'; end
if last_n == true, algorithm_name = 'last-n'; end
if kalman_like == true, algorithm_name = 'kalman-like'; end
if kalman == true, algorithm_name = 'kalman'; end
if quadratic_all == true, algorithm_name = 'quadratic-all'; end
if batch == true, algorithm_name = 'batch'; end
title_string = ['\color[rgb]{0.9 0.4 0.3}', algorithm_name, '\color[rgb]{0.25 0.25 0.25}'];
if batch || last_n, title_string = [title_string, '     batch = ', '\color[rgb]{0 0.6 0.3}', num2str(batch_size), '\color[rgb]{0.25 0.25 0.25}']; end
title_string = [title_string, '     \omega_2 = ', '\color[rgb]{0 0.6 0.3}', num2str(w2), '\color[rgb]{0.25 0.25 0.25}'];
if last_n, title_string = [title_string, '     \omega_3 = ', '\color[rgb]{0 0.6 0.3}', num2str(w3), '\color[rgb]{0.25 0.25 0.25}']; end
title_string = [title_string, '     noise = ', '\color[rgb]{0 0.6 0.3}', num2str(measurement_noise_std), '\color[rgb]{0.25 0.25 0.25}'];
title(title_string, 'interpreter','tex', 'FontWeight','Normal', 'fontsize', 20);

xlim([1, length(history) + 1]); 
yyaxis left; ax = gca; ax.YColor = [0.8 0.3 0]; ax.YLim = ylimit; ylabel('\mu(x) \pm \sigma(x)');
yyaxis right; ax = gca; ax.YColor = [0 0.5 0]; ax.YLim = [-0.015, max_value * 4]; ylabel('(J^TJ)^{1/2}');





