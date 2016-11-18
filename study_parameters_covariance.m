function [] = study_parameters_covariance(settings, Histories)

B = 3;
T = 3;
chisquare_val = 2.4477;
figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.7, 0.8]); hold on;
%set(gca,'position', [0.05 0.05 0.95 0.95], 'units','normalized');

frame_indices = 1:settings.num_frames;

index1 = 1; index2 = 2;
%index1 = 5; index2 = 1;

for i = 1:length(frame_indices)
    frame_index = frame_indices(i);   
   
    mean_h = zeros(2, 2);
    mean_mu = zeros(2, 1);
    run_indices = 1:settings.num_runs;
       
    %% draw inliers covariance elipces
    data = zeros(settings.num_runs, B + T);
    for run_index = run_indices
        data(run_index, :) = Histories{run_index}.x_batch(frame_index, :);
        mu = data(run_index, [index1, index2])';               
        H = squeeze(Histories{run_index}.full_covariance(frame_index, [index1, index2], [index1, index2]));
        sigma = 0.5 * inv(H); 
        h = inv(sigma);
        
        %[ellipse_points, ~] = get_covarince_elipse(sigma, chisquare_val);
        %plot(ellipse_points(:,1) + mu(1), ellipse_points(:,2) + mu(2), '-', 'lineWidth', 2, 'color',  [228, 244, 223]/255);
        
        mean_h = mean_h + h;
        mean_mu = mean_mu + mu;
    end
    
    %% plot data points
    scatter(data(:, index1), data(:, index2), 20, [1.0 0.45 0.3], 'o', 'filled');
    mean_h = mean_h / length(run_indices);
    mean_mu = mean_mu / length(run_indices);
    mean_sigma = inv(mean_h);
    [ellipse_points] = get_covarince_elipse(mean_sigma, chisquare_val);
    plot(ellipse_points(:,1) + mean_mu(1), ellipse_points(:,2) + mean_mu(2), '-', 'lineWidth', 2, 'color', [136, 187, 119]/255);
    
end
xlim([0, 6]); ylim([0, 6]); axis equal;
xlabel('\beta_1'); if i > 1, ylabel('\beta_2'); end

% xlim([-0.5, 2.2]); ylim([0, 6]); 
% xlabel('\beta_1'); if i > 1, ylabel('\beta_2'); end

