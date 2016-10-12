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

num_iters = 20;
to_display = false;

if (to_display), figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.8, 0.8]);
    axis off; axis equal; hold on;
end

settings.laplace_approx = false;
settings.last_n = false;
settings.kalman_like = false;
settings.kalman = false;
settings.quadratic_all = false;
settings.batch = true;
settings.independent = false;
settings.no_lm = false;

settings.batch_size = 5;

w2 = 1;

%% Tracking
for N = 1:num_frames
    %disp(N);
    
    %% Laplace approximation
    if (settings.laplace_approx)
        if N < 3
            X = X_init(1:(B + T) * N);
            [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_batch(X, segments0, joints, frames, N, D, settings.batch_size, w2), X, num_iters);
            h = zeros(3, 2, 2);
        else
            X_prev = [X; X_init((B + T) * (N - 1) + 1:(B + T) * N)];
            X0 = [X(1:(B + T) * (N - 2)); X_init((B + T) * (N - 2) + 1:(B + T) * N)];
            
            [xx_opt, J, h] = simple_problem_quadratic_two_lsqnonlin(X0, X_prev, h, frames, N, w2, num_iters);
            J = sqrtm(h); %disp([J' * J, h]);
            
            X = [X(1:(B + T) * (N - 2)); xx_opt];           
        end
    end
    
    %% Separate optimization
    if (settings.independent)
        X = zeros(N * (B + T), 1);
        J = zeros(N * (B + T), N * (B + T));
        for i = 1:N
            x = X_init((B + T) * (i - 1) + 1:(B + T) * i);
            [x, j] = my_lsqnonlin(@(x) sticks_finger_fg_single(x, segments0, joints, frames{i}), x, num_iters);
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
        
        if (settings.no_lm)
            lambdas = zeros((B + T) * N, 1);
            indices = [zeros(B, 1); ones(T, 1)]; indices = repmat(indices, N, 1);
            lambdas(indices == 0) = 0.1; lambdas(indices == 1) = 50;
            for iter = 1:num_iters
                [F, J] = sticks_finger_fg_batch(X, segments0, joints, frames, N, D, settings.batch_size, w2);
                delta = - (J' * J + diag(lambdas)) \ (J' * F);
                X = X + delta;
            end
        else
            [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_batch(X, segments0, joints, frames, N, D, settings.batch_size, w2), X, num_iters);
        end
        
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
        display_sticks_finger(segments, data_points, model_points, num_iters);
    end
    
    %% Save history
    history{N}.X = X;
    history{N}.JtJ = J' * J;
end


%% Display history
%{
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
            
            importance = history{j}.JtJ((B + T) * (k - 1) + b, (B + T) * (k - 1) + b);
            variance = min(5, 1/importance);
            myline([j + offset * k, history{j}.betas{k}(b) + 0.1 * sqrt(importance)], ...
                [j + offset * k, history{j}.betas{k}(b) - 0.1 * sqrt(importance)], [1, 0.9, 0.3], 2.7);
            
            mypoint([j + offset * k, history{j}.betas{k}(b)], [0.3, 0.6, 0.8], 15);
        end
    end
    set(gca, 'fontSize', 13); title(['beta ', num2str(b)]);
    xlim([1, length(history) + 1]);
    ylim([1.5, 5]);
    %legend({'beta-true', 'batch boundary', 'diag(JtJ(i, i))', 'beta'})
end
%}

