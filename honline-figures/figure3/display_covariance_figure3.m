function [] = display_covariance_figure3(settings, results_history, covariance_history, frame_certainty)
min_lim = -1; max_lim = 7;
chisquare_val = 1.5;%2.4477;
uncertain_background_color = [1; 0.98; 0.95];
certain_background_color = [0.96; 1; 0.93];

dark_red = [188, 58, 117]/255;
light_red = [230, 168, 169]/255;
dark_green = [61, 131, 119]/255;
light_green = [144, 194, 171]/240;
orange = [255, 173, 153]/255;

figure('units', 'normalized', 'outerposition', [0.1, 0.3, 0.8, 0.43]); hold on;
set(gca,'position', [0.05 0.05 0.95 0.95], 'units','normalized');
shifts = 0.02 + [0. 0.2, 0.4, 0.6, 0.8];

frame_indices = 1:5;

%% Display all data points
for i = 1:length(frame_indices)
    frame_index = frame_indices(i);
    
    %% set up subplot
    h = subplot('Position', [shifts(i), 0.05, 0.16, 0.85]); hold on; axis equal;
    
    %% display
    if sum(frame_certainty(frame_index, :)) >= 1
        color = certain_background_color;
    else
        color = uncertain_background_color;
    end
    set(gca,'color', color);
    set(gca,'color', 'w');
    
    mean_h = zeros(2, 2);
    mean_mu = zeros(2, 1);
    run_indices = 1:settings.num_runs;
    
    %% draw true value
    plot([min_lim, max_lim], [settings.betas_true(frame_index, 2), settings.betas_true(frame_index, 2)], 'lineWidth', 1.2, 'color',  [0.75, 0.75, 0.75]);
    plot([settings.betas_true(frame_index, 1), settings.betas_true(frame_index, 1)], [min_lim, max_lim], 'lineWidth', 1.2, 'color',  [0.75, 0.75, 0.75]);
    
    %% draw covariance elipces
    num_valid_points = length(run_indices);
    for run_index = run_indices
        mu = squeeze(results_history(run_index, frame_index, 1:2));
        sigma = 0.5 * squeeze(covariance_history(run_index, frame_index, 1:2, 1:2));
        h = inv(sigma);
        
        if (any(isnan(h))), num_valid_points = num_valid_points - 1; continue; end   
        
        mean_h = mean_h + h;
        mean_mu = mean_mu + mu;
    end
    
    
    %% plot data points
    scatter(results_history(run_indices, frame_index, 1), results_history(run_indices, frame_index, 2), 50, dark_green, 'o', 'filled');
    mean_h = mean_h / num_valid_points;
    mean_mu = mean_mu / num_valid_points;
    mean_sigma = inv(mean_h);
    [ellipse_points] = get_covarince_elipse(mean_sigma, chisquare_val);
    plot(ellipse_points(:,1) + mean_mu(1), ellipse_points(:,2) + mean_mu(2), '-', 'lineWidth', 4, 'color', light_green);    
   
    %% plot parameters
    xlim([min_lim, max_lim]); ylim([min_lim, max_lim]); box on; set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'linewidth', 1.5);
    set(gca, 'xcolor', [0.5, 0.5, 0.5]); set(gca, 'ycolor', [0.5, 0.5, 0.5]); 
end

%% Display title
% h = display_algorithm_title(settings, 'sticks_finger');
% set(h,'Position', [-18, 8.1,  0]);

