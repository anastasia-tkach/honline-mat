function [frames, thetas_true] = get_random_frames(N, num_samples, beta_true, measurement_noise_std, theta_noise_std)

B = 3;

theta_init = [0; 0; 0];
theta_true = theta_init;
thetas_true = zeros(N, B);

theta_limit = pi/12;
theta_certain = pi/4;
beat = 8;

% Generate random frames
theta_limit = 0;
for i = 1:N
    theta_true = theta_true + theta_noise_std * randn(B, 1); 
    theta_true(2) = max(theta_limit, theta_true(2));
    theta_true(3) = max(theta_limit, theta_true(3));
    thetas_true(i, :) = theta_true';
end

%{
% Generate random frames
for i = 1:(N - 4 * beat)/3
    theta_true = theta_true + theta_noise_std * randn(B, 1); 
    theta_true(2) = min(theta_limit, theta_true(2)); theta_true(2) = max(-theta_limit, theta_true(2));
    theta_true(3) = min(theta_limit, theta_true(3)); theta_true(3) = max(-theta_limit, theta_true(3));
    thetas_true(i, :) = theta_true';
end
last_frame_index = i;

% Add certain pose
thetas_true_2 = [linspace(theta_true(2), theta_certain, beat), linspace(theta_certain, 0, beat)];
thetas_true_3 = [linspace(theta_true(3), 0, beat), linspace(theta_certain, 0, beat)];
thetas_true(last_frame_index + 1:last_frame_index + 2 * beat, :) = [zeros(length(thetas_true_2), 1), thetas_true_2', thetas_true_3'];
last_frame_index = last_frame_index + 2 * beat;

% More random frames
theta_true = [0; 0; 0];
for i = last_frame_index + 1: 2 * (N - 2 * beat)/3
    theta_true = theta_true + theta_noise_std * randn(B, 1);    
    theta_true(2) = min(theta_limit, theta_true(2)); theta_true(2) = max(-theta_limit, theta_true(2));
    theta_true(3) = min(theta_limit, theta_true(3)); theta_true(3) = max(-theta_limit, theta_true(3));
    thetas_true(i, :) = theta_true';
end
last_frame_index = i;

% Add certain pose
thetas_true_2 = [linspace(theta_true(2), theta_certain, beat), linspace(theta_certain, 0, beat)];
thetas_true_3 = [linspace(theta_true(3), theta_certain, beat), linspace(theta_certain, 0, beat)];
thetas_true(last_frame_index + 1:last_frame_index + 2 * beat, :) = [zeros(length(thetas_true_2), 1), thetas_true_2', thetas_true_3'];
last_frame_index = last_frame_index + 2 * beat;

% More random frames
theta_true = [0; 0; 0];
for i = last_frame_index + 1:N
    theta_true = theta_true + theta_noise_std * randn(B, 1);    
    theta_true(2) = min(theta_limit, theta_true(2)); theta_true(2) = max(-theta_limit, theta_true(2));
    theta_true(3) = min(theta_limit, theta_true(3)); theta_true(3) = max(-theta_limit, theta_true(3));
    thetas_true(i, :) = theta_true';
end
%}

%% Generate data points
thetas_true(:, 1) = 0;
frames = cell(length(thetas_true), 1);
for i = 1:N
    [segments0, joints] = segments_and_joints_2D();
    [data_segments] = shape_2D(segments0, beta_true);
    [data_segments] = pose_2D(data_segments, joints, thetas_true(i, :)');
    [data_points] = sample_2D(data_segments, num_samples);
    for j = 1:length(data_points)
        data_points{j} = data_points{j} + measurement_noise_std * randn(2, 1);
    end
    frames{i} = data_points;
end