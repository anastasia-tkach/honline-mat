%clear; clc;
%rng default;

%% Parameters
num_samples = 10;
B = 3;
T = 3;
D = 3 * num_samples;
settings.measurement_noise_std = 0.07;
settings.beta_noise_std = 0.5;
settings.theta_noise_std = 0.15;
blocks = {[1, 2], [2, 3], [3, 4]};
beta_true = [3; 3; 3];
theta_init = [0; 0; 0];

[segments0, joints] = segments_and_joints_2D();

%% Scrip
theta_certain = [0, pi/3, pi/3];
theta_uncertain = [0, 0, 0];
thetas_true = [repmat(theta_uncertain, 4, 1); repmat(theta_certain, 4, 1); repmat(theta_uncertain, 7, 1)];
settings.num_frames = size(thetas_true, 1);

[frames, beta_init, thetas_init] = get_random_data_from_theta(beta_true, thetas_true, ...
    settings.beta_noise_std, settings.theta_noise_std, settings.measurement_noise_std, num_samples);
X_init = zeros((B + T) * settings.num_frames, 1);
for i = 1:settings.num_frames
    X_init((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = beta_init;
    X_init((B + T) * (i - 1) + B + 1:(B + T) * i) = thetas_init{i};
end

to_display = false;

if (to_display), figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.8, 0.8]);
    axis off; axis equal; hold on;
end

settings.quadratic_one = false;
settings.quadratic_two = false;
settings.kalman_like = false;
settings.batch = true;
settings.independent = false;

settings.batch_size = 2;%settings.num_frames;
settings.num_iters = 20;

settings.batch_independent = false;
settings.batch_online = false;
settings.batch_online_robust = true;

[settings, history] = set_batch_size(settings);

w2 = 1;
h = [];

%% Tracking
for N = 1:settings.num_frames
    %disp(N);
    
    %% Quadratic-one
    if (settings.quadratic_one && N >= 3)
        
        X = X_init((B + T) * (N - 2) + 1:(B + T) * N);        
        x_1 = history.x_batch(N - 1, 1:B + T)';
        x0 = history.x_batch(N - 1, (B + T) + 1:end)';
        
        [X, J, h] = sticks_finger_quadratic_one(X, x0, x_1, h, segments0, joints, frames, N, w2, settings);
        
        H = h;       
        H(1:B + T, 1:B + T) = diag(history.h_batch(N - 1, (B + T) * (settings.batch_size - 1) + 1:(B + T) * settings.batch_size));
    end
    
    %% Laplace approximation
    if (settings.quadratic_two && N >= 3)
        X = X_init((B + T) * (N - 2) + 1:(B + T) * N);
        
        x_1 = history.x_batch(N - 1, 1:B + T)';
        x0 = history.x_batch(N - 1, (B + T) + 1:end)';
        [X, J, h] = sticks_finger_laplace_approx(X, x0, x_1, h, segments0, joints, frames, N, w2, settings.num_iters);
        
        H = zeros(2 * (B + T), 2 * (B + T));
        a = h(1:B, 1:B); b = h(1:B, B + T + 1:B + T + B); c = h(B + T + 1:B + T + B, 1:B); d = h(B + T + 1:B + T + B, B + T + 1:B + T + B);
        H(B + T + 1:B + T + B, B + T + 1:B + T + B) = d - c * inv(a) * b;
        H(1:B + T, 1:B + T) = diag(history.h_batch(N - 1, (B + T) * (settings.batch_size - 1) + 1:(B + T) * settings.batch_size));
    end
    
    %% Kalman-like
    if (settings.kalman_like)
        X = X_init((B + T) * (N - 1) + 1:(B + T) * N);
        if (N == 1)
            x0 = [];
            JtJ = zeros(B, B);
        else
            x0 = history.x_batch(N - 1, :)';
        end
        
        [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_kalman_like(X, x0, segments0, joints, frames{N}, JtJ, N, w2), X, settings.num_iters);
        
        J1 = J(1:D, 1:B);
        JtJ = JtJ + (J1'* J1);
        H = [JtJ, zeros(B, B); zeros(B, B), zeros(B, B)];
    end
    
    %% Separate optimization
    if (settings.independent)
        for i = 1:N
            X = X_init((B + T) * (i - 1) + 1:(B + T) * i);
            [X, j] = my_lsqnonlin(@(X) sticks_finger_fg_data(x, segments0, joints, frames{i}), X, settings.num_iters);
        end
        H = j' * j;
    end
    
    %% Batch optimization
    if (settings.batch || (N < 3 && (settings.quadratic_two || settings.quadratic_one)))
        if N <= settings.batch_size
            X = X_init(1:(B + T) * N);
            x0 = [];
        else
            X = X_init((B + T) * (N - settings.batch_size) + 1:(B + T) * N);
            x0 = history.x_batch(N - 1, 1:B + T)';
        end
        
        [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_batch(X, x0, segments0, joints, frames, N, D, settings, w2), X, settings.num_iters);
        H = J' * J;        
    end
    
    %% Save new historyory
    if N <= settings.batch_size
        history.x_batch(N, :) = [zeros((B + T) * (settings.batch_size - N), 1); X];
        history.h_batch(N, :) = [zeros((B + T) * (settings.batch_size - N), 1); diag(H)];
    else
        indices = 1:(B + T) *  settings.batch_size;
        history.x_batch(N, :) = X(indices);
        history.h_batch(N, :) = diag(H);
    end
    
    %% Display
    if (to_display)
        betas_N = X(end - (B + T) + 1:end - T);
        thetas_N = X(end - T + 1:end);
        data_points = frames{N};
        [segments0, joints] = segments_and_joints_2D();
        [segments0] = shape_2D(segments0, betas_N);
        [segments] = pose_2D(segments0, joints, thetas_N);
        [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
        clf; hold on; axis off; axis equal; set(gcf,'color','w');
        display_sticks_finger(segments, data_points, model_points, settings.num_iters);
    end
end



