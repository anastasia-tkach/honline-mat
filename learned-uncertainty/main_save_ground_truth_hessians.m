clear; clc;
rng default;
B = 3;
T = 3;
beta = [3; 3; 3];
settings.beta_noise_std = 0;
settings.theta_noise_std = 0;
settings.measurement_noise_std = 1e-4;
[segments0, joints] = segments_and_joints_2D();
settings.num_iters = 20;
settings.num_samples = 10;
settings.num_runs = 100;

settings.data_model_energy = true;
settings.model_data_energy = false;
settings.silhouette_energy = false;
settings.display_iterations = false;

num_steps = 101;

thetas = zeros(num_steps, num_steps, T);
true_hessians = zeros(num_steps, num_steps, B, B);

%% Experiments
steps = linspace(-pi/2, pi/2, num_steps);
for t2 = 1:num_steps
    disp(t2);
    for t3 = 1:num_steps        
        %disp(t3);
        theta = [0; steps(t2); steps(t3)];
        thetas(t2, t3, :) = theta';
        mean_hessian = zeros(B, B);
        for run_index = 1:settings.num_runs
            
            [segments] = shape_2D(segments0, beta);
            [segments] = pose_2D(segments, joints, theta);
            [data_points] = sample_2D(segments, settings.num_samples);
            for j = 1:length(data_points)
                data_points{j} = data_points{j} + settings.measurement_noise_std * randn(2, 1);
            end
            
            X = [beta; theta];
            [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_data(X, segments0, joints, data_points, settings, 'cpp'), X, settings.num_iters);
            H = J' * J;
            sigma = inv(H(1:B, 1:B));
            mean_hessian = mean_hessian + inv(sigma);
        end
        mean_hessian = mean_hessian / settings.num_runs;        
        true_hessians(t2, t3, :, :) = mean_hessian;
    end
end


%% display
rng default;
chisquare_val = 2.4477;
figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.7, 0.8]); hold on; %axis equal;
for i = 1:10:size(true_hessians, 1)
    for j = 1:10:size(true_hessians, 2)
        h = squeeze(true_hessians(i, j, :, :));
        
        sigma = 0.002 * inv(h(1:2, 1:2));
        mypoint(true_thetas(i, j, 2:3), [1.0 0.45 0.3], 20);
        [ellipse_points] = get_covarince_elipse(sigma, chisquare_val);
        plot(ellipse_points(:,1) + true_thetas(i, j, 2), ellipse_points(:,2) + true_thetas(i, j, 3), '-', 'lineWidth', 2, 'color', [136, 187, 119]/255);
    end
end

theta_lookup = [0.0161, -0.2371, 0.2548];%[0; pi * rand() - pi/2; pi * rand() - pi/2; ];

[h, theta_closest] = lookup_ground_truth_hessian(theta_lookup, true_thetas, true_hessians);

sigma = 0.002 * inv(h(1:2, 1:2));
mypoint(theta_lookup(2:3), 'b', 20);
[ellipse_points] = get_covarince_elipse(sigma, chisquare_val);
plot(ellipse_points(:,1) + theta_closest(2), ...
    ellipse_points(:,2) + theta_closest(3), '-', 'lineWidth', 2, 'color', 'b');


xlim([-4, 4]); ylim([-2.5, 2.5]);