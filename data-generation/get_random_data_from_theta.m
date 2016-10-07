function [frames, beta_init, thetas_init] = get_random_data_from_theta(beta_true, thetas_true, beta_noise_std, theta_noise_std, measurement_noise_std, num_samples)

B = 3;
T = 3;

num_frames = length(thetas_true);
thetas_init = cell(num_frames, 1);
beta_init = beta_true + beta_noise_std * randn(B, 1);

thetas = cell(num_frames, 1);
betas = cell(num_frames, 1);
frames = cell(length(thetas_true), 1);
for i = 1:size(thetas_true, 1)
    thetas_init{i} = thetas_true(i, :)' + theta_noise_std * randn(T, 1);
    thetas{i} = thetas_init{i};
    betas{i} = beta_init;
    [segments0, joints] = segments_and_joints_2D();
    [data_segments] = shape_2D(segments0, beta_true);
    [data_segments] = pose_2D(data_segments, joints, thetas_true(i, :)');
    [data_points] = sample_2D(data_segments, num_samples);
    for j = 1:length(data_points)
        data_points{j} = data_points{j} + measurement_noise_std * randn(2, 1);
    end
    frames{i} = data_points;
end