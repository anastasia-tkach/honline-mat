function [frames, betas_true, betas_init, thetas_init] = get_random_data_from_theta(betas_true, thetas_true, settings)

B = 3;
T = 3;
num_frames = length(thetas_true);
betas_init = cell(num_frames, 1);
thetas_init = cell(num_frames, 1);

if (numel(betas_true) == 3)
    beta_init = betas_true + settings.beta_bias + settings.beta_noise_std * randn(B, 1);
    beta_init = max(beta_init, 1);
    for i = 1:num_frames
        betas_init{i} = beta_init;
    end  
    betas_true = repmat(betas_true', num_frames, 1);
else
    for i = 1:num_frames
        betas_init{i} = betas_true(i, :)' + settings.beta_noise_std * randn(B, 1);
    end
end


thetas = cell(num_frames, 1);
betas = cell(num_frames, 1);
frames = cell(length(thetas_true), 1);
for i = 1:size(thetas_true, 1)
    thetas_init{i} = thetas_true(i, :)' + settings.theta_noise_std * randn(T, 1);
    thetas{i} = thetas_init{i};
    betas{i} = betas_init{i};
    [segments0, joints] = segments_and_joints_2D();
    [data_segments] = shape_2D(segments0, betas_true(i, :)');
    [data_segments] = pose_2D(data_segments, joints, thetas_true(i, :)');
    [data_points] = sample_2D(data_segments, settings.num_samples);
    for j = 1:length(data_points)
        data_points{j} = data_points{j} + settings.measurement_noise_std * randn(2, 1);
    end
    frames{i} = data_points;
end