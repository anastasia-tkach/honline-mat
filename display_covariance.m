function [] = display_covariance(settings, results_history, covariance_history, frame_certainty)
min_lim = -1; max_lim = 7;
chisquare_val = 1.5;%2.4477;
uncertain_background_color = [1; 0.98; 0.95];
certain_background_color = [0.96; 1; 0.93];

figure('units', 'normalized', 'outerposition', [0.1, 0.3, 0.8, 0.43]); hold on;
set(gca,'position', [0.05 0.05 0.95 0.95], 'units','normalized');
shifts = 0.02 + [0. 0.2, 0.4, 0.6, 0.8];

if settings.num_frames <= 5
    frame_indices = 1:settings.num_frames;
else
    %frame_indices = [8, 9, 10, 12, 14];
    %frame_indices = [1, 2, 3, 7, 12];
    frame_indices = 1:5;
    %frame_indices = [2, 5, 7, 11, 13];
    %frame_indices = [70, 71, 73, 75, 146];
end

%% Find an outlier (look at the last frame)
data = squeeze(results_history(:, frame_indices(1), 1:2));
mean_data = mean(data);
data = data - repmat(mean_data, size(data, 1), 1);
sigma_data = cov(data);
distances = zeros(size(data, 1), 1);
for i = 1:size(data, 1)
    d = data(i, :)';
    distances(i) = sqrt(d' * inv(sigma_data) * d);
end
outlier_indices = find(distances > chisquare_val);
%outlier_indices = find(distances < 1.5);
%outlier_indices = find((distances > 2.3) .* (distances < chisquare_val) );
%outlier_indices = outlier_indices(randi([1, length(outlier_indices)], 1, 1));

%% Display all data points
for i = 1:length(frame_indices)
    frame_index = frame_indices(i);
    
    %% set up subplot
    %h = subplot(1, 5, i); hold on; p = get(h, 'pos'); disp(p);
    %0.1300    0.1100    0.1237    0.8150
    h = subplot('Position', [shifts(i), 0.05, 0.16, 0.85]); hold on; axis equal;
    
    %% display
    if sum(frame_certainty(frame_index, :)) >= 1
        color = certain_background_color;
    else
        color = uncertain_background_color;
    end
    set(gca,'color', color);

    mean_h = zeros(2, 2);
    mean_mu = zeros(2, 1);
    run_indices = 1:settings.num_runs;
    
    %% Draw true value
    plot([min_lim, max_lim], [settings.betas_true(frame_index, 2), settings.betas_true(frame_index, 2)], 'lineWidth', 1, 'color',  [0.75, 0.75, 0.75]);
    plot([settings.betas_true(frame_index, 1), settings.betas_true(frame_index, 1)], [min_lim, max_lim], 'lineWidth', 1, 'color',  [0.75, 0.75, 0.75]);
    
    %% draw inliers covariance elipces
    num_valid_points = length(run_indices);
    for run_index = run_indices
        mu = squeeze(results_history(run_index, frame_index, 1:2));
        sigma = 0.5 * squeeze(covariance_history(run_index, frame_index, 1:2, 1:2));
        h = inv(sigma);
        
        if (any(isnan(h)))
            num_valid_points = num_valid_points - 1;
            continue
        end
        
        [ellipse_points, ok] = get_covarince_elipse(sigma, chisquare_val);
        if ok
            plot(ellipse_points(:,1) + mu(1), ellipse_points(:,2) + mu(2), '-', 'lineWidth', 2, 'color',  [228, 244, 223]/255);
        end
        
        mean_h = mean_h + h;
        mean_mu = mean_mu + mu;
    end
    
    %% draw outliner covariance elipces
    %{
    for outlier_index = outlier_indices'
        mu = squeeze(results_history(outlier_index, frame_index, 1:2));
        sigma = 0.5 * squeeze(covariance_history(outlier_index, frame_index, 1:2, 1:2));
        [ellipse_points, ok] = get_covarince_elipse(sigma, chisquare_val);
        if ok
            plot(ellipse_points(:,1) + mu(1), ellipse_points(:,2) + mu(2), '-', 'lineWidth', 2, 'color',  [0.8, 0.8, 0.8]);
        end
    end
    %}
    
    %% plot data points
    scatter(results_history(run_indices, frame_index, 1), results_history(run_indices, frame_index, 2), 20, [1.0 0.45 0.3], 'o', 'filled');
    %scatter(results_history(outlier_indices, frame_index, 1), results_history(outlier_indices, frame_index, 2), 20, [0.5, 0.5, 0.5], 'o', 'filled');
    mean_h = mean_h / num_valid_points;
    mean_mu = mean_mu / num_valid_points;
    mean_sigma = inv(mean_h);
    [ellipse_points] = get_covarince_elipse(mean_sigma, chisquare_val);
    plot(ellipse_points(:,1) + mean_mu(1), ellipse_points(:,2) + mean_mu(2), '-', 'lineWidth', 2, 'color', [136, 187, 119]/255);
    
    %% compute data covariance
    if (settings.num_runs > 1)
        data = squeeze(results_history(:, frame_index, 1:2));
        if (size(data, 2) == 1) data = data'; end
        mean_data = mean(data, 1);
        data = data - repmat(mean_data, size(data, 1), 1);
        sigma_data = cov(data);
        [ellipse_points] = get_covarince_elipse(sigma_data, chisquare_val);
        plot(ellipse_points(:,1) + mean_data(1), ellipse_points(:,2) + mean_data(2), '-', 'lineWidth', 2, 'color', [255, 173, 135]/255);
    end
    
    %% plot parameters
    xlim([min_lim, max_lim]); ylim([min_lim, max_lim]);
    title(['frame ', num2str(frame_index)], 'FontWeight','Normal'); set(gca, 'fontSize', 10); set(gca,'fontname','Cambria');
    xlabel('\beta_1'); if i > 1, ylabel('\beta_2'); end
end

%% Display title
h = display_algorithm_title(settings, 'sticks_finger');
set(h,'Position', [-18, 8.1,  0]);

