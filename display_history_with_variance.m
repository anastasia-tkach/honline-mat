function [] = display_history_with_variance(means, standard_deviations, importance_means, importance_standard_deviations, x_true, ylimit, settings, N, w2, frame_centrainty, problem_type)

w = 0.6;
line_color = [1, 0.7, 0.6];
point_color = [1.0 0.5 0.3]; 

certain_line_color = [0.75, 0.9, 0.7];
certain_point_color =  [0.25, 0.75, 0.35];

point_size = 30;
line_width = 3;

f = figure('units', 'normalized', 'outerposition', [0.1, 0.3, w, 0.55]); hold on;
set(gca,'position', [0.06 0.06 0.87 0.85], 'units','normalized');


max_value = max(importance_means(:) + importance_standard_deviations(:));

%% Put results to vector
online_means = zeros(N, 1);
for k = 1:N
    online_means(k) = means(k, k);
end

%% Print foreground
yyaxis left;
plot(0:N + 1, x_true * ones(N + 2, 1), 'lineWidth', 2.3, 'color', [1, 0.5, 0.4], 'lineStyle', '-.');
plot(online_means, 'lineWidth', 2.7, 'color', [245, 207, 184]/255, 'lineStyle', '-');


%% Plot
for k = 1:N
    if (frame_centrainty(k) == 0)
        current_line_color = (1 - k / N ) * [0.95, 0.88, 0.88] + k / N * line_color;
        current_point_color = (1 - k / N ) * [0.85, 0.75, 0.75] + k / N * point_color;
    else
        current_line_color = 0.5 * [0.95, 0.88, 0.88] + 0.5* certain_line_color;
        current_point_color = 0.5 * [0.85, 0.75, 0.75] + 0.5 * certain_point_color;
    end
    if k == N
        if (frame_centrainty(k) == 0)
            current_line_color = line_color;
            current_point_color = [1.0 0.45 0.3];
        else
            current_line_color =  [157, 216, 105]/255;
            current_point_color =  [0.1, 0.7, 0.3];
        end
    end
    
    % results
    yyaxis left;
    myline([k, means(k, k) + standard_deviations(k, k)], [k, means(k, k) - standard_deviations(k, k)], current_line_color, line_width);
    mypoint([k, means(k, k)], current_point_color, point_size);
    
    % importance
    yyaxis right;
    y_position = 0.08;
    myline([k, y_position], [k, y_position + importance_means(k, k) + importance_standard_deviations(k, k)], [0.75, 0.9, 0.7], line_width);
    myline([k, y_position], [k, y_position + importance_means(k, k)], [0.65, 0.8, 0.6], line_width);
    
end



%% Print algorithm parameters

display_algorithm_title(N, w2, ylimit, max_value, settings, problem_type);


xlim([0, N + 1]); 