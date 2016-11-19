%{
M = cell(2, 1);
for i = 1:2
    M{i} = @(x) [x * i, x];
end
f = @(x) arrayfun(@(i) M{i}(x), (1:2)', 'uniformOutput', true);
%}
%{
f = @(x) [];
for i = 1:100
    f = @ (x) [f(x); i];
end
%}


close all;
clear; clc;
settings.data_model_energy = true;
settings.model_data_energy = false;
settings.silhouette_energy = false;
settings.display_iterations = false;
rng(11);

%% Parameters
num_samples = 3;
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
theta_init = [pi/3; pi/3; pi/3];
theta_true = theta_init;
N = 1;
num_iters = 20;

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
   
    %[F, G, H] = sticks_finger_eg_single([beta; theta], segments0, joints, data_points, 'numerical');
    [F_, G_, H_] = sticks_finger_fg_data([beta; theta], segments0, joints, data_points, settings, 'analytical');
    [F_cpp, G_cpp, H_cpp] = sticks_finger_fg_data([beta; theta], segments0, joints, data_points, settings, 'cpp');
    
    imagesc(G_ - G_cpp); colorbar;
    disp([G_, G_cpp]);
    
    %disp([F; F_; F_cpp]);
    
    %for i = 1:size(F, 1)
    %    disp([G(i, :); G_(i, :); G_cpp(i, :)]);
    %end
    %for i = 1:6
    %   disp([H(i, :); H_(i, :)]);
    %end
    
    %% lsqnonlin
    %{
    x = [beta; theta];
    [x, j] = my_lsqnonlin(@(x) sticks_finger_fg_single(x, segments0, joints, data_points), x, num_iters);
    
    %[F, J] = sticks_finger_fg_single(x, segments0, joints, data_points);
    %df = my_gradient(@(x) sticks_finger_fg_single(x, segments0, joints, data_points), x);
    
    beta = x(1:B); theta = x(B + 1:B + T);
    %}
    
    %% fminunc
    %{
    x = [beta; theta];
    options = optimoptions(@fminunc,'Algorithm','trust-region', 'SpecifyObjectiveGradient', true, 'HessianFcn', 'objective');
    [x, fval, exitflag, output, grad, hessian] = fminunc(@(x) sticks_finger_eg_single(x, segments0, joints, data_points), x, options);
    beta = x(1:B); theta = x(B + 1:B + T);   
    %}   
    
    %% display
    %{
    [segments0, joints] = segments_and_joints_2D();
    [segments0] = shape_2D(segments0, beta);
    [segments] = pose_2D(segments0, joints, theta);
    [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
    figure(1); clf; hold on; axis off; axis equal; set(gcf,'color','w');
    display_sticks_finger(segments, data_points, model_points, num_iters);
    %}
    
    %{
    for iter = 1:num_iters
        
        %% Compute correspondences
        [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
        
        %{
        if (iter == num_iters)
            clf; hold on; axis off; axis equal; set(gcf,'color','w');
            display_sticks_finger(segments, data_points, model_points, iter);
        end
        %}
        
        %% Compute Jacobians
        %%{
        %[F, J_theta] = jacobian_ik_2D(beta, theta, segments, joints, model_points, data_points, segment_indices);
        %[F, J_beta] = jacobian_shape_2D(beta, theta, segments, model_points, data_points, segment_indices);
        %J = [J_beta, J_theta];
        %F_ = F;
        %J_ = J;
        
        [F, J] = jacobian_shape_pose_cpp_wrapper(beta, theta, segments, joints, model_points, data_points, segment_indices);
                        

        I = 0.1 * eye(T + B, T + B); I(end - T + 1:end, end -T + 1:end) = 50 * eye(T, T);
        
        LHS = J' * J + I;
        delta = LHS \ (J' * F);
        
        beta = beta + delta(1:B);
        theta = theta + delta(B + 1:end);
        
        %% Pose
        [shaped_segments] = shape_2D(segments0, beta);
        [segments] = pose_2D(shaped_segments, joints, theta);
        
    end
    %}
    
end
