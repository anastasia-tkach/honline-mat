function [] = display_covariance(settings, results_history, covariance_history, frame_certainty)

chisquare_val = 2.4477;
mean_h = zeros(2, 2);
mean_mu = zeros(2, 1);
uncertain_background_color = [1; 0.98; 0.95];
certain_background_color = [0.96; 1; 0.93];
for frame_index = 1:settings.num_frames
%for frame_index = [2, 5, 8, 11, 14]
    figure; hold on;  axis equal;
    if sum(frame_certainty(frame_index)) >= 1
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
        plot(r_ellipse(:,1) + mu(1), r_ellipse(:,2) + mu(2), '-', 'lineWidth', 2, 'color', 0.93 * color);
        
        mean_h = mean_h + h;
        mean_mu = mean_mu + mu;
    end
    scatter(results_history(:, frame_index, 1), results_history(:, frame_index, 2), 20, [0.1, 0.7, 0.5], 'o', 'filled');
    mean_h = mean_h / settings.num_runs;
    mean_mu = mean_mu / settings.num_runs;
    mean_sigma = inv(h);
    [r_ellipse] = get_covarince_elipse(mean_sigma, chisquare_val);
    plot(r_ellipse(:,1) + mean_mu(1), r_ellipse(:,2) + mean_mu(2), '-', 'lineWidth', 2, 'color', [1, 0.4, 0]);
    xlim([-2, 8]); ylim([-2, 8]); title(['frame ', num2str(frame_index)]); set(gca, 'fontSize', 12); set(gca,'fontname','Cambria');
    xlabel('\beta_1');
    ylabel('\beta_2');
end
