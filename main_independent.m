clear; clc;
close all;

%% Parameters
settings.num_samples = 10;
B = 3;
T = 3;
settings.measurement_noise_std = 0.03;
settings.beta_bias = [0; 0; 0];
settings.beta_noise_std = 0.5;
settings.theta_noise_std = 0.1;
blocks = {[1, 2], [2, 3], [3, 4]};
[segments0, joints] = segments_and_joints_2D();

settings.data_model_energy = true;
settings.model_data_energy = false;
settings.silhouette_energy = false;
settings.display_converged = false;
settings.display_iterations = false;
settings.display_jacobian = false;

beta_true = [3; 3; 3];
theta_true = [0; 0; 0];
beta = beta_true + settings.beta_noise_std * randn;
theta = theta_true + settings.theta_noise_std * randn;
X = [beta; theta];

settings.num_frames = 1;
settings.num_iters = 1;

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
    display_jacobian(X, J, settings, N);
    
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


