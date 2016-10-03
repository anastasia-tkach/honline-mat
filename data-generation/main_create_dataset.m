clear; clc; close all;
N = 70;
B = 3;
beta_true = [3; 3; 3];
num_samples = 15;
D = 3 * num_samples;
measurement_noise_std = 0.07; 
theta_noise_std = 0.15;
beta_noise_std = 0.7;

num_items = 100;

dataset_path = 'C:\Users\t-antka\OneDrive - Microsoft\Data\CalibrationDataset\';
for i = 1:num_items
    [frames, thetas_true] = get_random_frames(N, num_samples, beta_true, measurement_noise_std, theta_noise_std);
    beta_init = [0.3; 0.3; 0.3] + beta_noise_std * rand(B, 1);
    for k = 1:B
        v = randi([0, 1]);
        if v == 0, beta_init(k) = -1 * beta_init(k); end
    end
    beta_init = beta_true + beta_init;
    name_suffix = sprintf('%03d', i);
    save([dataset_path, 'frames_', name_suffix], 'frames');
    save([dataset_path, 'thetas_true_', name_suffix], 'thetas_true');
    save([dataset_path, 'beta_init_', name_suffix], 'beta_init');
end