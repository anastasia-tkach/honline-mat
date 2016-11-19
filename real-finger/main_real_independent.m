rng default;
%rng(2);
clear; clc;
%close all;

%% Parameters
settings.num_samples = 10;
B = 3;
T = 3;
settings.measurement_noise_std = 0.01;
settings.beta_bias = [0; 0; 0];
settings.beta_noise_std = 1;
settings.theta_noise_std = 0.05;
blocks = {[1, 2], [2, 3], [3, 4]};
[segments0, joints] = segments_and_joints_3D();

settings.data_model_energy = true;
settings.model_data_energy = false;
settings.silhouette_energy = false;
settings.display_converged = true;
settings.display_iterations = false;
settings.display_jacobian = false;

beta_true = [3; 3; 3];
theta_true = [pi/8; pi/8; pi/8];
radii_true = [1.0; 0.8; 0.6; 0.5];
beta = beta_true + settings.beta_noise_std * randn;
theta = theta_true + settings.theta_noise_std * randn;
radii = radii_true;
X = [beta; theta];

settings.num_frames = 1;
settings.num_iters = 5;


for N = 1:settings.num_frames
    
    %% Initialize
    theta_true = theta_true + settings.theta_noise_std * randn(size(theta));
    [data_segments] = shape_3D(segments0, beta_true);
    [data_segments] = pose_3D(data_segments, joints, theta_true);
    [data_points] = sample_3D(data_segments, radii, blocks);
    for i = 1:length(data_points)
        data_points{i} = data_points{i} + settings.measurement_noise_std * randn(3, 1);
    end
    
    %% Optimize
    
    [X, J] = my_lsqnonlin(@(X) real_finger_data(X, segments0, joints, radii, blocks, data_points), X, settings.num_iters);
    %[F, J] = real_finger_data(X, segments0, joints, radii, blocks, data_points);
    %return;
    
end

%%{
beta = X(end - (B + T) + 1:end - T);
theta = X(end - T + 1:end);
[segments0, joints] = segments_and_joints_3D();
[segments0] = shape_3D(segments0, beta);
[segments] = pose_3D(segments0, joints, theta);
[model_points, segment_indices] = compute_correspondences_3D(segments, radii, blocks, data_points);
figure; clf; hold on; axis off; axis equal; set(gcf,'color','w');
display_finger_3D(segments, radii, blocks, data_points, model_points);
drawnow; pause(0.05);
%%}
