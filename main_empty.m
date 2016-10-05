%close all;
clear; clc;
rng(2572);

%% Parameters
num_samples = 15;
B = 3;
T = 3;
D = 3 * num_samples;
measurement_noise_std = 0.07;
beta_noise_std = 0.5;
theta_noise_std = 0.15;
blocks = {[1, 2], [2, 3], [3, 4]};

%% Get data
beta_true = [3; 3; 3];
beta_init = beta_true + beta_noise_std * randn;
theta_init = [0; 0; 0];
theta_true = theta_init;
N = 50;
num_iters = 7;

%% Initialization
beta = beta_init;
%beta = beta_true;
theta = theta_init;
[segments0, joints] = segments_and_joints_2D();
[segments0] = shape_2D(segments0, beta);
[segments] = pose_2D(segments0, joints, theta);


%% Tracking
%figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.8, 0.8]); axis off; axis equal; hold on;

for n = 1:N
    %disp(beta');
    %% Get data
    
    theta_true = theta_true + theta_noise_std * randn(size(theta));
    [data_segments] = shape_2D(segments0, beta_true);
    [data_segments] = pose_2D(data_segments, joints, theta_true);
    [data_points] = sample_2D(data_segments, num_samples);
    for i = 1:length(data_points)
        data_points{i} = data_points{i} + measurement_noise_std * randn(2, 1);
    end
    
    for iter = 1:num_iters
        
        %% Compute correspondences
        [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
        
        if (iter == num_iters)
            clf; hold on; axis off; axis equal; set(gcf,'color','w');
            display_sticks_finger(segments, data_points, model_points, iter);            
        end
        
        %% Compute Jacobians
        %%{
        [F, J_theta] = jacobian_ik_2D(segments, joints, model_points, data_points, segment_indices);
        %[F_, J_] = jacobian_pose_cpp_wrapper(segments, joints, model_points, data_points, segment_indices);       
        [F, J_beta] = jacobian_shape_2D(segments, model_points, data_points, segment_indices);
        
        J = [J_beta, J_theta];

        I = 0.1 * eye(T + B, T + B); I(end - T + 1:end, end -T + 1:end) = 50 * eye(T, T);
        
        LHS = J' * J + I;
        delta = LHS \ (J' * F);        
        
        beta = beta + delta(1:B);
        theta = theta + delta(B + 1:end);
        %%}
        
        %% Compute Jacobians
        %{
        [F, J_theta] = jacobian_ik_2D(segments, joints, model_points, data_points, segment_indices);       

        J = J_theta;
        I = 50 * eye(T, T);
        
        LHS = J' * J + I;
        delta = LHS \ (J' * F);        

        theta = theta + delta;
        %}
        
        %% Pose
        [shaped_segments] = shape_2D(segments0, beta);
        [segments] = pose_2D(shaped_segments, joints, theta);        
        
    end
    
    
end
