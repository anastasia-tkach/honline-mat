function [] = display_covariance(settings, results_history, covariance_history, frame_certainty)

chisquare_val = 2.4477;
mean_h = zeros(2, 2);
mean_mu = zeros(2, 1);
uncertain_background_color = [1; 0.98; 0.95];
certain_background_color = [0.96; 1; 0.93];

figure('units', 'normalized', 'outerposition', [0.1, 0.3, 0.8, 0.4]); hold on;
set(gca,'position', [0.05 0.05 0.95 0.95], 'units','normalized');
shifts = 0.02 + [0. 0.2, 0.4, 0.6, 0.8];

frame_indices = [2, 5, 8, 11, 14];
for i = 1:length(frame_indices)
    frame_index = frame_indices(i);
    
    %% set up subplot       
    %h = subplot(1, 5, i); hold on; p = get(h, 'pos'); disp(p);
    %0.1300    0.1100    0.1237    0.8150
    h = subplot('Position', [shifts(i), 0.05, 0.16, 0.95]); hold on; axis equal;
    
    %% display
    if sum(frame_certainty(frame_index, :)) >= 1
        color = certain_background_color;
    else
        color = uncertain_background_color;
    end
    set(gca,'color', color);
    for run_index = 1:settings.num_runs
        mu = squeeze(results_history(run_index, frame_index, 1:2));
        h = squeeze(covariance_history(run_index, frame_index, 1:2, 1:2));
        sigma = inv(h);
        [r_ellipse] = get_covarince_elipse(sigma, chisquare_val);
        plot(r_ellipse(:,1) + mu(1), r_ellipse(:,2) + mu(2), '-', 'lineWidth', 2, 'color',  [228, 244, 223]/255);
        
        mean_h = mean_h + h;
        mean_mu = mean_mu + mu;
    end
    scatter(results_history(:, frame_index, 1), results_history(:, frame_index, 2), 20, [1.0 0.45 0.3], 'o', 'filled');
    mean_h = mean_h / settings.num_runs;
    mean_mu = mean_mu / settings.num_runs;
    mean_sigma = inv(h);
    [r_ellipse] = get_covarince_elipse(mean_sigma, chisquare_val);
    plot(r_ellipse(:,1) + mean_mu(1), r_ellipse(:,2) + mean_mu(2), '-', 'lineWidth', 2, 'color', [136, 187, 119]/255);
    
    %% compute data covariance
    data = squeeze(results_history(:, frame_index, 1:2));
    mean_data = mean(data);
    data = data - repmat(mean_data, size(data, 1), 1);
    sigma_data = cov(data);
    [r_ellipse] = get_covarince_elipse(sigma_data, chisquare_val);
    plot(r_ellipse(:,1) + mean_data(1), r_ellipse(:,2) + mean_data(2), '-', 'lineWidth', 2, 'color', [255, 173, 135]/255);
    
    %% plot parameters
    xlim([-1, 7]); ylim([-1, 7]); title(['frame ', num2str(frame_index)]); set(gca, 'fontSize', 12); set(gca,'fontname','Cambria');
    xlabel('\beta_1'); if i > 1, ylabel('\beta_2'); end
end


