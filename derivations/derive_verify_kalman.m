clear; clc; close all;
num_betas = 30;

data_path = 'C:/Users/tkach/Desktop/Test/';

estimated_hessians = cell(0, 1);
measured_hessians = cell(0, 1);
estimated_values = cell(0, 1);
measured_values = cell(0, 1);

for frame_id = 0:30:540
    fileID = fopen([data_path, 'estimated_hessian-', num2str(frame_id), '.txt'], 'r');
    estimated_hessian = fscanf(fileID, '%f');
    estimated_hessian = reshape(estimated_hessian, num_betas, num_betas)';
    estimated_hessians{end + 1} = estimated_hessian;
    fclose(fileID);
    
    fileID = fopen([data_path, 'measured_hessian-', num2str(frame_id), '.txt'], 'r');
    measured_hessian = fscanf(fileID, '%f');
    measured_hessian = reshape(measured_hessian, num_betas, num_betas)';
    measured_hessians{end + 1} = measured_hessian;
    fclose(fileID);
    
    fileID = fopen([data_path, 'estimated_values-', num2str(frame_id), '.txt'], 'r');
    estimated_value = fscanf(fileID, '%f');
    estimated_values{end + 1} = estimated_value;
    fclose(fileID);
    
    fileID = fopen([data_path, 'measured_values-', num2str(frame_id), '.txt'], 'r');
    measured_value = fscanf(fileID, '%f');
    measured_values{end + 1} = measured_value;
    fclose(fileID);
end

N = length(measured_hessians);

indices = [3, 4, 5, 18, 19, 20] + 1;
estimated_hessian_test = eye(num_betas, num_betas);
estimated_value_test = [38.6827 30.8717 17.9383 42.5040 23.3138 13.7786 44.5005 26.3163 17.7812 38.6084 27.8816 16.7835 40.5208 27.0978 18.1370 14.7393 13.1366 -1.4419 31.4316 59.9664 8.6118 16.2903 55.7356 12.9243 -1.3718 55.0835 12.6442 -11.2736 45.1948 7.9830]';

%% sub-hessian
estimated_hessian_sub = estimated_hessian_test(indices, indices);
estimated_value_sub = estimated_value_test(indices);

EV = cell(N, 1);
MV = cell(N, 1);
EH = cell(N, 1);
MH = cell(N, 1);

for i = 1:length(measured_hessians)
    measured_hessian_sub = measured_hessians{i}(indices, indices);
    measured_values_sub = measured_values{i}(indices);
    
    estimated_value_sub = inv(estimated_hessian_sub + measured_hessian_sub) * (estimated_hessian_sub * estimated_value_sub + measured_hessian_sub * measured_values_sub);
    
    %estimated_value_sub = inv(diag(diag(estimated_hessian_sub)) + diag(diag(measured_hessian_sub))) * ...
    %    (diag(diag(estimated_hessian_sub)) * estimated_value_sub +  diag(diag(measured_hessian_sub)) * measured_values_sub);
    
    estimated_hessian_sub = estimated_hessian_sub + measured_hessian_sub;
    disp([estimated_values{i}(indices), estimated_value_sub]);
    
    EV{i} = estimated_value_sub;
    MV{i} = measured_values_sub;
    EH{i} = estimated_hessian_sub;
    MH{i} = measured_hessian_sub;
end

beta_id = 1;
for beta_id = 1:6
    estimated_values_history = zeros(length(measured_hessians), 1);
    measured_values_history = zeros(length(measured_hessians), 1);
    estimated_certainty_history = zeros(length(measured_hessians), 1);
    measured_certainty_history = zeros(length(measured_hessians), 1);
    for i = 1:length(measured_hessians)
        estimated_values_history(i) = EV{i}(beta_id);
        measured_values_history(i) = MV{i}(beta_id);
        estimated_certainty_history(i) = EH{i}(beta_id, beta_id);
        measured_certainty_history(i) = MH{i}(beta_id, beta_id);
    end
    
    %{
    figure; hold on;
    plot(estimated_values_history ./ max(estimated_values_history), 'lineWidth', 2);
    plot(measured_values_history ./ max(estimated_values_history), 'lineWidth', 2);
    plot(estimated_certainty_history ./ max(estimated_certainty_history), 'lineWidth', 2);
    plot(measured_certainty_history ./ max(measured_certainty_history), 'lineWidth', 2);
    %}
end

%% Plot covariance at first frame
frame_id = 1;
for i = [1, 4, 5, 6]
    for j = [1, 4, 5, 6]
        if (i >= j), continue; end
        figure; axis equal; hold on;
        
        h1 = EH{frame_id}([i, j], [i, j]);
        x1 = EV{frame_id}([i, j]);
        h2 = MH{frame_id + 1}([i, j], [i, j]);
        x2 = MV{frame_id + 1}([i, j]);        
        
        x3 = EV{frame_id + 1}([i, j]);
        h3 = EH{frame_id + 1}([i, j], [i, j]);
       
        x3_indep = inv(h1 + h2) * (h1 * x1 + h2 * x2);
        h3_indep = h1 + h2;
        
        indices = [1, 4, 5, 6];
        H1 = EH{frame_id}(indices, indices);
        X1 = EV{frame_id}(indices);
        H2 = MH{frame_id + 1}(indices, indices) + 1000 * eye(length(indices), length(indices));
        X2 = MV{frame_id + 1}(indices);
        
        
        X3 = inv(H1 + H2) * (H1 * X1 + H2 * X2);
        H3 = H1 + H2;
        XX3 = zeros(6, 1);
        XX3(indices) = X3;
        x3_regul = XX3([i, j]);
        HH3 = zeros(6, 6);
        HH3(indices, indices) = H3;
        h3_regul = HH3([i, j], [i, j]);
        
        %% display
        estimated_sigma = 200000 * inv(h1);
        estimated_mu = x1;
        ellipse_points = get_covarince_elipse(estimated_sigma, 1);
        plot(ellipse_points(:,1) + estimated_mu(1), ellipse_points(:,2) + estimated_mu(2), '-', 'lineWidth', 2, 'color',  'r');
        mypoint(estimated_mu, 'r', 30);
        
        measured_sigma = 200000 * inv(h2);
        measured_mu = x2;
        ellipse_points = get_covarince_elipse(measured_sigma, 1);
        plot(ellipse_points(:,1) + measured_mu(1), ellipse_points(:,2) + measured_mu(2), '-', 'lineWidth', 2, 'color',  'b');
        mypoint(measured_mu, 'b', 30);        
      
        estimated_sigma_new = 200000 * inv(h3);
        estimated_mu_new = x3;
        ellipse_points = get_covarince_elipse(estimated_sigma_new, 1);
        plot(ellipse_points(:,1) + estimated_mu_new(1), ellipse_points(:,2) + estimated_mu_new(2), '-', 'lineWidth', 2, 'color',  'g');
        mypoint(estimated_mu_new, 'g', 30);  
        
        %% regularized
        estimated_sigma_regul = 200000 * inv(h3_regul);
        estimated_mu_regul = x3_regul;
        ellipse_points = get_covarince_elipse(estimated_sigma_regul, 1);
        plot(ellipse_points(:,1) + estimated_mu_regul(1), ellipse_points(:,2) + estimated_mu_regul(2), '-', 'lineWidth', 2, 'color',  'm');
        mypoint(estimated_mu_regul, 'm', 30);  
        
        %% independent
        estimated_sigma_indep = 200000 * inv(h3_indep);
        estimated_mu_indep = x3_indep;
        ellipse_points = get_covarince_elipse(estimated_sigma_indep, 1);
        plot(ellipse_points(:,1) + estimated_mu_indep(1), ellipse_points(:,2) + estimated_mu_indep(2), '-', 'lineWidth', 2, 'color',  'c');
        mypoint(estimated_mu_indep, 'c', 30);  
        
        title([num2str(i), ', ', num2str(j)]);
    end
end







