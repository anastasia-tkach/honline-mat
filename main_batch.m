%close all;
clear; clc;
rng default;

%% Parameters
num_samples = 7;
B = 3;
T = 3;
D = 3 * num_samples;
measurement_noise_std = 0.07;  % 0.01 - converges to true beta, 0.1 - does not find true beta
beta_noise_std = 0.5; %3
theta_noise_std = 0.15;
blocks = {[1, 2], [2, 3], [3, 4]};
beta_true = [3; 3; 3];
theta_init = [0; 0; 0];

R = 0.2 * eye(D + B + T, D + B + T);
Q = diag([0 * [1; 1; 0]; theta_noise_std * ones(T, 1)]); %0.0001
P = diag([beta_noise_std * ones(B, 1); zeros(T, 1)]);

%% Scrip
%%{
theta_certain = [0, pi/3, pi/3];
theta_uncertain = [0, 0, 0];
thetas_true = [repmat(theta_uncertain, 3, 1); repmat(theta_certain, 4, 1); repmat(theta_uncertain, 8, 1)];
N = length(thetas_true);
[frames, beta_init, thetas_init] = get_random_data_from_theta(beta_true, thetas_true, ...
    beta_noise_std, theta_noise_std, measurement_noise_std, num_samples);

num_frames = length(thetas_true);
num_iters = 20;
to_display = false;

init_from_previous_frame = false;
with_proper_lm = true;
batch_size = num_frames;

run_kalman_filter = false;
r = 0.2;
R = r * eye(D + B + T, D + B + T);
Q = diag([zeros(B, 1); theta_noise_std * ones(T, 1)]);
P = diag([beta_noise_std * ones(B, 1); zeros(T, 1)]);
P_ = P + Q; C_ = inv(P_);

history = {};
J_previous = [];
F_previous = [];

if (to_display), figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.8, 0.8]);
    axis off; axis equal; hold on;
end

w2 = 10;

%% Tracking
for N = 1:num_frames
    disp(N);
    
    %% Initialization
    for i = max(1, N - batch_size + 1):N
        if init_from_previous_frame
            if i == 1
                betas{1} = beta_init;
                thetas{1} = theta_init;
            else
                betas{i} = history{N - 1}.betas{i - 1};
                thetas{i} = history{N - 1}.thetas{i - 1};
            end
        else
            thetas{i} = thetas_init{i};
            betas{i} = beta_init;
        end
    end
    beta0 = betas{N}; theta0 = thetas{N};
    
    %if (to_display), figure('units', 'normalized', 'outerposition', [0.1, 0.3, 0.8, 0.5]); axis off; axis equal; hold on; end
    
    %% Batch optimization
    X = zeros(N * (B + T), 1);
    for i = 1:N
        X((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = betas{i};
        X((B + T) * (i - 1) + B + 1:(B + T) * i) = thetas{i};
    end
    
    [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_batch(X, frames, N, D, batch_size, w2), X, num_iters);    
    
    betas = cell(N, 1);
    thetas = cell(N, 1);
    for i = 1:N
        betas{i} = X((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B);
        thetas{i} = X((B + T) * (i - 1) + B + 1:(B + T) * i);
    end
    
    %% Display
    if (to_display)
        data_points = frames{N};
        [segments0, joints] = segments_and_joints_2D();
        [segments0] = shape_2D(segments0, betas{i});
        [segments] = pose_2D(segments0, joints, thetas{i});
        [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
        figure(1); clf; hold on; axis off; axis equal; set(gcf,'color','w');
        display_sticks_finger(segments, data_points, model_points, iter);
    end
    
    %% Save history
    history{N}.P = P;
    history{N}.betas = betas;
    history{N}.thetas = thetas;
    history{N}.JtJ = J' * J;
    history{N}.mean = betas{N};
end

%figure; hold on; plot(energy(2:end), 'lineWidth', 2);
%figure; hold on; plot(beta(:, 1), 'lineWidth', 2); plot(beta(:, 2), 'lineWidth', 2); plot(beta(:, 3), 'lineWidth', 2);

%% Display online equivalent
%{
figure_borders = [0.05 0.08 0.93 0.90];
means = zeros(N, B);
vars = zeros(N, B);
trues = zeros(N, B);
%beta_std = zeros(N, B);
for i = 1:length(history)
    means(i, :) = history{i}.betas{i};
    trues(i, :) = beta_true;
    if (run_kalman_filter)
        vars(i, :) = diag(history{i}.P(1:B, 1:B));
    else
        vars(i, :) = diag(history{i}.JtJ((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B, (B + T) * (i - 1) + 1:(B + T) * (i - 1) + B));
    end
    %beta_std(i, :) = history{i}.beta_std;
end

figure('units', 'normalized', 'outerposition', [0.1, 0.2, 0.33, 0.7]); hold on;
set(gca,'position', figure_borders, 'units','normalized');
for i = 1:B - 1
    h = subplot(B - 1, 1, i); hold on;
    p = get(h, 'pos'); set(h, 'pos', [0.05, p(2) - 0.06, 0.9, 0.43]);
    set(gca,'XTick',[]);
    plot(1:length(history), means(:, i), '.-', 'lineWidth', 2, 'markersize', 13);
    plot(1:length(history), trues(:, i), 'lineWidth', 2);
    plot(1:length(history), beta_init(i) * ones(length(history), 1), 'lineWidth', 2, 'lineStyle', '-.');
    
    if (run_kalman_filter)
        plot(1:length(history), means(:, i) + vars(:, i).^0.5, 'lineWidth', 2, 'color', [0.55, 0.75, 0.8]);
        plot(1:length(history), means(:, i) - vars(:, i).^0.5, 'lineWidth', 2, 'color', [0.55, 0.75, 0.8]);
    end
    
    %plot(1:length(history), trues(:, i) + 0.5 * beta_variance_threshold * ones(length(history), 1), 'lineWidth', 1, 'color', [0.7, 0.9, 0.5]);
    %plot(1:length(history), trues(:, i) + 0.5 * beta_std(:, i).^2, 'lineWidth', 1, 'color', [0.4, 0.7, 0.5]);
    
    ylim([1.5, 4.5]); xlim([0, length(history)]);
    set(gca, 'fontSize', 13); title(['beta ', num2str(i)]);
end

if exist('video_writer', 'var'), video_writer.close(); end
%}

%% Display history
offset = 0.065;
figure('units', 'normalized', 'outerposition', [0.1, 0.2, 0.55, 0.7]); hold on;
set(gca,'position', [0.05 0.08 0.93 0.90], 'units','normalized');
for b = 1:B - 1
    h = subplot(B - 1, 1, b); hold on;
    p = get(h, 'pos'); set(h, 'pos', [0.05, p(2) - 0.06, 0.9, 0.43]);
    set(gca,'XTick',[]);
    plot(1:length(history) + 1, beta_true(b) * ones(length(history) + 1, 1), 'lineWidth', 2, 'color', [1, 0.5, 0.4], 'lineStyle', '-.');
    for j = 1:N
        myline([j, 0], [j, 8], [0.85, 0.85, 0.85], 2);
        for k = 1:j
            if (~run_kalman_filter)
                importance = history{j}.JtJ((B + T) * (k - 1) + b, (B + T) * (k - 1) + b);
                variance = min(5, 1/importance);
                myline([j + offset * k, history{j}.betas{k}(b) + 0.1 * sqrt(importance)], ...
                    [j + offset * k, history{j}.betas{k}(b) - 0.1 * sqrt(importance)], [1, 0.9, 0.3], 2.7);
            end
            mypoint([j + offset * k, history{j}.betas{k}(b)], [0.3, 0.6, 0.8], 15);
        end
    end
    set(gca, 'fontSize', 13); title(['beta ', num2str(b)]);
    xlim([1, length(history) + 1]);
    ylim([1.5, 5]);
    %legend({'beta-true', 'batch boundary', 'diag(JtJ(i, i))', 'beta'})
end
