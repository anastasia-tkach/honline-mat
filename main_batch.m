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
num_frames = size(thetas_true, 1);

[frames, beta_init, thetas_init] = get_random_data_from_theta(beta_true, thetas_true, ...
    settings.beta_noise_std, settings.theta_noise_std, settings.measurement_noise_std, num_samples);
X_init = zeros((B + T) * num_frames, 1);
for i = 1:num_frames
    X_init((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = beta_init;
    X_init((B + T) * (i - 1) + B + 1:(B + T) * i) = thetas_init{i};
end

to_display = false;

if (to_display), figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.8, 0.8]);
    axis off; axis equal; hold on;
end

settings.quadratic_one = false;
settings.quadratic_two = true;
settings.kalman_like = false;
settings.kalman = false;
settings.batch = false;
settings.independent = false;

settings.batch_size = 2;
settings.num_iters = 20;

settings.batch_independent = false;
settings.batch_robust = false;

w2 = 1;
X = [];
% history.X_current = zeros(num_frames, settings.batch_size);
% history.X_final = zeros(num_frames, 1);
% history.H_current = zeros(num_frames, settings.batch_size);
% history.H_final = zeros(num_frames, 1);

%% Tracking
for N = 1:num_frames
    %disp(N);
    
    %% Quadratic-one
    if (settings.quadratic_one)
        if N < 3
            X = X_init(1:(B + T) * N);
            [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_batch(X, segments0, joints, frames, N, D, settings, w2), X, settings.num_iters);
            h = zeros(3, 2, 2);
        else
            X_prev = [X; X_init((B + T) * (N - 1) + 1:(B + T) * N)];
            X0 = [X(1:(B + T) * (N - 2)); X_init((B + T) * (N - 2) + 1:(B + T) * N)];
            [xx_opt, J, h] = sticks_finger_quadratic_one(X0, X_prev, h, segments0, joints, frames, N, w2, settings);
            H = h;
            H(1:B + T, 1:B + T) = history{N - 1}.JtJ((B + T) * (N - 2) + 1:(B + T) * (N - 1), (B + T) * (N - 2) + 1:(B + T) * (N - 1));
            X = [X(1:(B + T) * (N - 2)); xx_opt];
        end
    end
    
    %% Laplace approximation
    if (settings.quadratic_two)
        if N < 3
            X = X_init(1:(B + T) * N);
            [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_batch(X, segments0, joints, frames, N, D, settings, w2), X, settings.num_iters);
            h = zeros(3, 2, 2);
        else
            X_prev = [X; X_init((B + T) * (N - 1) + 1:(B + T) * N)];
            X0 = [X(1:(B + T) * (N - 2)); X_init((B + T) * (N - 2) + 1:(B + T) * N)];
            
            [xx_opt, J, h] = sticks_finger_laplace_approx(X0, X_prev, h, segments0, joints, frames, N, w2, settings.num_iters);
            %H = h;
            H = zeros(2 * (B + T), 2 * (B + T));
            a = h(1:B, 1:B);
            b = h(1:B, B + T + 1:B + T + B);
            c = h(B + T + 1:B + T + B, 1:B);
            d = h(B + T + 1:B + T + B, B + T + 1:B + T + B);
            H(B + T + 1:B + T + B, B + T + 1:B + T + B) = d - c * inv(a) * b;
            H(1:B + T, 1:B + T) = history{N - 1}.JtJ((B + T) * (N - 2) + 1:(B + T) * (N - 1), (B + T) * (N - 2) + 1:(B + T) * (N - 1));
            
            X = [X(1:(B + T) * (N - 2)); xx_opt];
        end
    end
    
    %% Kalman-like
    if (settings.kalman_like)
        x = X_init((B + T) * (N - 1) + 1:(B + T) * N);
        if (N == 1)
            x_prev = [];
            JtJ = zeros(B, B);
        else
            x_prev = X((B + T) * (N - 2) + 1:(B + T) * (N - 1));
        end
        
        [x, J] = my_lsqnonlin(@(x) sticks_finger_fg_kalman_like(x, x_prev, segments0, joints, frames{N}, JtJ, N, w2), x, settings.num_iters);
        
        J1 = J(1:D, 1:B);
        JtJ = JtJ + (J1'* J1);
        X = [X; x];
    end
    
    %% Kalman
    if (settings.kalman)
        if (N == 1)
            x = X_init((B + T) * (N - 1) + 1:(B + T) * N);
            C = zeros(B + T, B + T);
        else
            x = X_init((B + T) * (N - 1) + 1:(B + T) * N);
            x(1:B) = X((B + T) * (N - 2) + 1:(B + T) * (N - 2) + B);
        end
        [x, C] = sticks_finger_kalman(x, C, segments0, joints, frames{N}, settings.num_iters);
        X = [X; x]; JtJ = C(1:B, 1:B);
    end
    
    %% Separate optimization
    if (settings.independent)
        X = zeros(N * (B + T), 1);
        J = zeros(N * (B + T), N * (B + T));
        for i = 1:N
            x = X_init((B + T) * (i - 1) + 1:(B + T) * i);
            [x, j] = my_lsqnonlin(@(x) sticks_finger_fg_data(x, segments0, joints, frames{i}), x, settings.num_iters);
            X((B + T) * (i - 1) + 1:(B + T) * i) = x;
            J(D * (i - 1) + 1:D * i, (B + T) * (i - 1) + 1:(B + T) * i) = j;
        end
    end
    
    %% Batch optimization
    if (settings.batch)        
        if N <= settings.batch_size
            X = X_init(1:(B + T) * N);
        else
            X = [X(1:(B + T) * (N - settings.batch_size)); X_init((B + T) * (N - settings.batch_size) + 1:(B + T) * N)];
        end        
        [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_batch(X, segments0, joints, frames, N, D, settings, w2), X, settings.num_iters);
        
    end
    %% Display
    if (to_display)
        betas_N = X((B + T) * (N - 1) + 1:(B + T) * (N - 1) + B);
        thetas_N = X((B + T) * (N - 1) + B + 1:(B + T) * N);
        data_points = frames{N};
        [segments0, joints] = segments_and_joints_2D();
        [segments0] = shape_2D(segments0, betas_N);
        [segments] = pose_2D(segments0, joints, thetas_N);
        [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
        clf; hold on; axis off; axis equal; set(gcf,'color','w');
        display_sticks_finger(segments, data_points, model_points, settings.num_iters);
    end
    
    %% Save history
    %history.X_current(N, :) = X;
    history{N}.X = X;
    if settings.batch || settings.independent || (N < 3 && (settings.quadratic_two || settings.quadratic_one))
        history{N}.JtJ = zeros((B + T) * N, (B + T) * N);
        if N > 1
            history{N}.JtJ(1:(B + T) * max(1, N - settings.batch_size), 1:(B + T) * max(1, N - settings.batch_size)) = ...
                history{N - 1}.JtJ(1:(B + T) * max(1, N - settings.batch_size), 1:(B + T) * max(1, N - settings.batch_size));
        end
        JtJ = J' * J;
        history{N}.JtJ((B + T) * max(0, N - settings.batch_size) + 1:(B + T) * N, (B + T) * max(0, N - settings.batch_size) + 1:(B + T) * N) = ...
            JtJ((B + T) * max(0, N - settings.batch_size) + 1:(B + T) * N, (B + T) * max(0, N - settings.batch_size) + 1:(B + T) * N);
    end
    if N >= 3 && (settings.quadratic_two || settings.quadratic_one)
        history{N}.JtJ = zeros((B + T) * N, (B + T) * N);
        history{N}.JtJ(1:(B + T) * (N - 2), 1:(B + T) * (N - 2)) = history{N - 1}.JtJ(1:(B + T) * (N - 2), 1:(B + T) * (N - 2));
        history{N}.JtJ((B + T) * (N - 2) + 1:(B + T) * N, (B + T) * (N - 2) + 1:(B + T) * N) = H;
    end
    if settings.kalman_like || settings.kalman
        history{N}.JtJ = zeros((B + T) * N, (B + T) * N);
        if (N > 1), history{N}.JtJ(1:(B + T) * (N - 1), 1:(B + T) * (N - 1)) = history{N - 1}.JtJ; end
        history{N}.JtJ((B + T) * (N - 1) + 1: (B + T) * (N - 1) + B, (B + T) * (N - 1) + 1: (B + T) * (N - 1) + B) = JtJ;
    end
end



