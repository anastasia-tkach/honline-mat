clear; close all; clc;
%% Parameters
num_samples = 15;
B = 3;
T = 3;
D = 3 * num_samples;
N = 80;
blocks = {[1, 2], [2, 3], [3, 4]};

to_display = true;
to_debug = false;
num_iters = 7;

beta_true = [3; 3; 3];
beta_init = beta_true;
theta = [0; 0; 0];
measurement_noise_std = 0.07;
theta_noise_std = 0.15;
r = 0.2;
theta_true = theta + theta_noise_std * randn(T, 1);

%% Initialization
beta = beta_init;
[segments0, joints] = segments_and_joints_2D();
[segments0] = shape_2D(segments0, beta);
[segments] = pose_2D(segments0, joints, theta);

%% Tracking
if (to_display), figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.8, 0.8]);
    axis off; axis equal; hold on;
end
if (to_debug),
    figure('units', 'normalized', 'outerposition', [0.1, 0.2, 0.33, 0.7]); hold on;
    set(gca,'position', [0.05 0.08 0.93 0.90], 'units','normalized');
end
history = {};
for n = 1:N
    %disp(beta');
    %% Get data
    %%{
    theta_true = theta_true + theta_noise_std * randn(size(theta));
    for i = 1:T, if (theta_true(i)) < -0.1, theta_true(i) = -0.1; end; end
    [data_segments] = shape_2D(segments0, beta_true);
    [data_segments] = pose_2D(data_segments, joints, theta_true);
    [data_points] = sample_2D(data_segments, num_samples);
    for i = 1:length(data_points)
        data_points{i} = data_points{i} + measurement_noise_std * randn(2, 1);
    end
    %%}
    beta0 = beta;
    theta0 = theta;
    
    for iter = 1:num_iters
        dx = [beta; theta] - [beta0; theta0];
        %disp(beta');
        %% Compute correspondences
        [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
        
        %% Display
        if (to_debug && iter == num_iters)
            figure(2); clf; hold on; axis off; axis equal; set(gcf,'color','w');
            display_history(history, beta_init, beta_variance_threshold, N, B, T);
        end
        if (to_display && iter == num_iters)
            figure(1); clf; hold on; axis off; axis equal; set(gcf,'color','w');
            display_kalman_2D(segments, data_points, model_points, iter);            
        end
        
        %% Data Jacobians
        [F, J] = jacobian_ik_2D(segments, joints, model_points, data_points, segment_indices);

        %% Joint Limits
        F_limits = zeros(T, 1);
        J_limits = zeros(T, T);
        w_limits = 10000;
        for j = 1:T
            if theta(j) < 0
                F_limits(j) = sqrt(w_limits) * (0 - theta(j));
                J_limits(j, j) = sqrt(w_limits) * (1);
            end
        end
        
        %% Solve
        I = 50 * eye(T, T);
        LHS = J' * J + J_limits' * J_limits + I + D;
        RHS = J' * F +  J_limits' * F_limits;
        delta = LHS \ RHS;
        
        %% Update
        theta = theta + delta;
        [shaped_segments] = shape_2D(segments0, beta);
        [segments] = pose_2D(shaped_segments, joints, theta);
        
        %% Stopping condition
        %disp(F' * F);
        %if F' * F < 0.1, break; end
    end    
    
    %% History
    h = length(history);
    history{h + 1}.mean = [beta; theta];
    history{h + 1}.true = [beta_true; theta_true];
    history{h + 1}.energy(1) = F(segment_indices == 1)' * F(segment_indices == 1);
    history{h + 1}.energy(2) = F(segment_indices == 2)' * F(segment_indices == 2);
    
    if (to_debug), waitforbuttonpress; end
end


