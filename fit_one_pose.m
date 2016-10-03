function [beta, theta] = fit_one_pose(beta_init, beta_true, theta_init, theta_true, measurement_noise_std, display)

%% Parameters
num_samples = 15;
B = length(beta_true);
T = length(theta_true);
blocks = {[1, 2], [2, 3], [3, 4]};
num_iters = 7;

%% Pose model
[segments0, joints] = segments_and_joints_2D();
[segments0] = shape_2D(segments0, beta_init);
[segments] = pose_2D(segments0, joints, theta_init);

%% Pose data
[data_segments] = shape_2D(segments0, beta_true);
[data_segments] = pose_2D(data_segments, joints, theta_true);
[data_points] = sample_2D(data_segments, num_samples);
for i = 1:length(data_points)
    data_points{i} = data_points{i} + measurement_noise_std * randn(2, 1);
end

%% Tracking
if (display), figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.8, 0.8]); axis off; axis equal; hold on; end
beta = beta_init;
theta = theta_init;
for iter = 1:num_iters
    
    %% Compute correspondences
    [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
    
    %% Display 
    if (display && iter == 1) 
        display_kalman_2D(segments, data_points, model_points, iter); 
    end
    
    %% Compute Jacobians
    [~, J_theta] = jacobian_ik_2D(segments, joints, model_points, data_points, segment_indices);
    [F, J_beta] = jacobian_shape_2D(segments, model_points, data_points, segment_indices);
    
    %% Solve
    J = [J_beta, J_theta];
    I = 0.1 * eye(T + B, T + B);
    I(end - T + 1:end, end -T + 1:end) = 500 * I(end - T + 1:end, end - T + 1:end);
    
    J1 = J; J2 = sqrt(I);
    F1 = F; F2 = zeros(B + T, 1);
    H = [J1; J2];  F = [F1; F2];
    
    delta = (H' * H + 0.01 * eye(T + B, T + B)) \ (H' * F);
    
    %% Update
    beta = beta + delta(1:B);
    theta = theta + delta(B + 1:end);
    [shaped_segments] = shape_2D(segments0, beta);
    [segments] = pose_2D(shaped_segments, joints, theta);
    
    if (display), disp(F' * F); end
end



