clear; clc;
close all;
rng('default')
%% Parameters
settings.num_samples = 10;
B = 3;
T = 3;
settings.measurement_noise_std = 0.03;
settings.beta_bias = [0; 0; 0];
settings.beta_noise_std = 0.0;
settings.theta_noise_std = 0.0;
blocks = {[1, 2], [2, 3], [3, 4]};
[segments0, joints] = segments_and_joints_2D();

settings.data_model_energy = true;
settings.model_data_energy = false;
settings.silhouette_energy = false;
settings.display_converged = false;
settings.write_video = false;
settings.display_iterations = false;
settings.display_jacobian = false;

beta_true = [3; 3; 3];
theta_true = [0;-pi/3; 0]; %-pi/15, -pi/3
beta = beta_true + settings.beta_noise_std * randn;
theta = theta_true + settings.theta_noise_std * randn;
X = [beta; theta];

settings.num_frames = 1;
settings.num_iters = 100;

for N = 1:settings.num_frames
    
    %% Initialize
    theta_true = theta_true + settings.theta_noise_std * randn(size(theta));
    [data_segments] = shape_2D(segments0, beta_true);
    [data_segments] = pose_2D(data_segments, joints, theta_true);
    [data_points] = sample_2D(data_segments, settings.num_samples);
    for i = 1:length(data_points)
        data_points{i} = data_points{i} + settings.measurement_noise_std * randn(2, 1);
    end
    
    %% Optimize
    
    [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_data(X, segments0, joints, data_points, settings, 'cpp'), X, settings.num_iters);
    
    
    %% Display
    %display_jacobian(X, J, settings, N);
    
    if (settings.display_converged)
        beta = X(end - (B + T) + 1:end - T);
        theta = X(end - T + 1:end);
        data_points = data_points;
        [segments0, joints] = segments_and_joints_2D();
        [segments0] = shape_2D(segments0, beta);
        [segments] = pose_2D(segments0, joints, theta);
        [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
        clf; hold on; axis off; axis equal; set(gcf,'color','w');
        display_sticks_finger(segments, data_points, model_points);
        drawnow; pause(0.05); %waitforbuttonpress
        
    end
end

%% Display energy as a function of beta1
dark_red = [188, 58, 117]/255;
light_red = [217, 154, 143]/255;
dark_green = [33, 159, 126]/255;
light_green = [144, 194, 171]/230;
orange = [255, 119, 51] / 255;

X_opt = X;
beta1_opt = X_opt(1);
[f_opt, j_opt, h_opt] = sticks_finger_fg_data(X_opt, segments0, joints, data_points, settings, 'cpp');

H_opt = 2 * (j_opt' * j_opt);
for i = 1:size(f_opt, 1)
    H_opt = H_opt + 2 * f_opt(i) * squeeze(h_opt(i, :, :));
end

beta1 =linspace(X_opt(1) - 1.5, X_opt(1) + 1.5, 2000);
F = zeros(length(beta1), 1);
F_quad_approx = zeros(length(beta1), 1);
F_quad_approx2 = zeros(length(beta1), 1);
for i = 1:length(beta1)
    X(1) = beta1(i);
    [f, j, h] = sticks_finger_fg_data(X, segments0, joints, data_points, settings, 'cpp');
    F(i) = f' * f;
    
    F_quad_approx(i) = f_opt' * f_opt + (beta1(i) - beta1_opt)' * H_opt(1, 1) * (beta1(i) - beta1_opt);
end

figure; hold on;
plot(beta1, F_quad_approx, 'lineWidth', 2.2, 'color', light_red, 'linestyle', '-.');
plot(beta1, F, 'lineWidth', 3, 'color', dark_green);
scatter(beta1_opt, (f_opt' * f_opt), 50, 'filled', 'markerFaceColor', dark_red);
set(gca,'fontsize', 13);

f_minus = sticks_finger_fg_data([beta1_opt - 1; X_opt(2:end)], segments0, joints, data_points, settings, 'cpp');
f_plus = sticks_finger_fg_data([beta1_opt + 1; X_opt(2:end)], segments0, joints, data_points, settings, 'cpp');
scatter(beta1_opt - 1, (f_minus' * f_minus), 50, 'filled', 'markerFaceColor', light_red);
scatter(beta1_opt + 1, (f_plus' * f_plus), 50, 'filled', 'markerFaceColor', light_red);

xlim([X_opt(1) - 1.5, X_opt(1) + 1.5]);
ylim([0, 17]); box on; set(gca,'xtick',[]); set(gca,'ytick',[]); set(gca,'linewidth', 1.5);
set(gca, 'xcolor', [0.5, 0.5, 0.5]);
set(gca, 'ycolor', [0.5, 0.5, 0.5]); 


%% Display
theta = X_opt(end - T + 1:end); beta = X_opt(end - (B + T) + 1:end - T); 

beta(1) = X_opt(1);
display_sticks_finger_with_corresp(beta, theta, data_points, blocks)

beta(1) = X_opt(1) - 1; 
display_sticks_finger_with_corresp(beta, theta, data_points, blocks)

beta(1) = X_opt(1) + 1; 
display_sticks_finger_with_corresp(beta, theta, data_points, blocks)

