function [] = display_empirical_variance(means, standard_deviations, importance_means, importance_standard_deviations, x_true, ylimit, settings, N, w2, frame_centrainty, problem_type)

w = N * 0.040909;
offset = 1/N;
line_color = [1, 0.85, 0.5];
point_color = [1.0 0.5 0.3]; %[0.3, 0.6, 0.8];

certain_line_color = [0.75, 0.9, 0.7];
certain_point_color =  [0.25, 0.75, 0.35];

f = figure('units', 'normalized', 'outerposition', [0.1, 0.3, w, 0.55]); hold on;
set(gca,'position', [0.06 0.06 0.87 0.85], 'units','normalized');
for j = 1:N
    if (frame_centrainty(j) == 0)
        rectangle('Position',[j, ylimit(1), 1, ylimit(2) - ylimit(1)],'FaceColor',[1; 0.98; 0.95],'EdgeColor','none')
    else
        rectangle('Position',[j, ylimit(1), 1, ylimit(2) - ylimit(1)],'FaceColor',[0.96; 1; 0.93],'EdgeColor','none')
    end
end
plot(1:N + 1, x_true * ones(N + 1, 1), 'lineWidth', 2, 'color', [1, 0.5, 0.4], 'lineStyle', '-.');

max_value = max(importance_means(:) + importance_standard_deviations(:));

for j = 1:N
    
    for k = 1:j
        if (frame_centrainty(k) == 0)
            current_line_color = (1 - k / j ) * [0.95, 0.88, 0.88] + k / j * line_color;
            current_point_color = (1 - k / j ) * [0.85, 0.75, 0.75] + k / j * point_color;
        else
            current_line_color = (1 - k / j ) * [0.95, 0.88, 0.88] + k / j * certain_line_color;
            current_point_color = (1 - k / j ) * [0.85, 0.75, 0.75] + k / j * certain_point_color;
        end
        point_size = 20;
        line_width = 3.2;
        if k == j
            if (frame_centrainty(k) == 0)
                current_line_color = [1, 0.75, 0.3];
                current_point_color = [1.0 0.45 0.3];
            else
                current_line_color =  [157, 216, 105]/255;
                current_point_color =  [0.1, 0.7, 0.3];
            end
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
        y_position = 0.08;
        myline([j + offset * k, y_position], ...
            [j + offset * k, y_position + importance_means(j, k) + importance_standard_deviations(j, k)], [0.75, 0.9, 0.7], 3.2);
        myline([j + offset * k, y_position], ...
            [j + offset * k, y_position + importance_means(j, k)], [0.65, 0.8, 0.6], 3.2);  
        
    end
end

%% Print algorithm parameters

display_algorithm_title(N, w2, ylimit, max_value, settings, problem_type);

