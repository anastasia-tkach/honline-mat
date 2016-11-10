% clear; clc;
% close all;
% rng default;
global video_writer;

%% Parameters
settings.num_samples = 10;
B = 3;
T = 3;
settings.measurement_noise_std = 0.07;
settings.beta_bias = [0; 0; 0];
settings.beta_noise_std = 0.5;
settings.theta_noise_std = 0.15;
blocks = {[1, 2], [2, 3], [3, 4]};
beta_true = [3; 3; 3];
theta_init = [0; 0; 0];
[segments0, joints] = segments_and_joints_2D();

%% Scrip
theta_certain_1 = [0, pi/3, 0];
theta_certain_2 = [0, 0, pi/3];
theta_certain_12 = [0, pi/3, pi/3];
theta_semicertain = [0, pi/60, pi/60];
theta_uncertain = [0, 0, 0];
tact = 3;
%thetas_true = [repmat(theta_certain_12, 1, 1); repmat(theta_uncertain, 3, 1)];
%thetas_true = [repmat(theta_certain_1, 7, 1); repmat(theta_certain_2, 7, 1)];
thetas_true = [repmat(theta_uncertain, tact, 1); repmat(theta_certain_1, tact, 1); repmat(theta_uncertain, tact, 1); repmat(theta_certain_2, tact, 1);  repmat(theta_uncertain, tact, 1)];

%thetas_true = [repmat(theta_uncertain, 5, 1); repmat(theta_certain_1, 2, 1); repmat(theta_certain_2, 2, 1);  repmat(theta_uncertain, 75, 1)];
%thetas_true = [repmat(theta_uncertain, 4, 1); repmat(theta_certain_12, 4, 1); repmat(theta_uncertain, 7, 1)];

settings.num_frames = size(thetas_true, 1);

[frames, beta_init, thetas_init] = get_random_data_from_theta(beta_true, thetas_true, settings);
X_init = zeros((B + T) * settings.num_frames, 1);
for i = 1:settings.num_frames
    X_init((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = beta_init;
    X_init((B + T) * (i - 1) + B + 1:(B + T) * i) = thetas_init{i};
end

%% Algorithm
settings.quadratic_one = false;
settings.quadratic_two = false;
settings.kalman_like = true;
settings.kalman_two = false;
settings.batch = false;
settings.batch_simulation = false;
settings.independent = false;

settings.batch_size = 2;

%% Parameters
settings.num_iters = 20;

settings.batch_independent = false;
settings.batch_online = true;
settings.batch_online_robust = false;
settings.batch_online_robust_tau = 1;
settings.quadratic_two_maximization = false;
settings.quadratic_two_marginalization = true;

settings.uniform_shape_prior = false;
settings.constant_sum_shape_prior = false;
settings.data_model_energy = true;
settings.model_data_energy = false;
settings.silhouette_energy = false;
settings.ground_truth_hessians = false;

settings.w1 = 1;
settings.w2 = 1;
settings.w4 = 1;

%% Display
settings.display_covariance = true;
settings.display_converged = false;
settings.display_iterations = false;
settings.display_jacobian = false;

settings.write_video = false;
if (settings.display_converged || settings.display_iterations),
    figure('units', 'normalized', 'outerposition', [0.25, 0.275, 0.45, 0.7]);
    axis off; axis equal; hold on;
end
if settings.write_video
    video_writer = VideoWriter('C:\Users\tkach\Desktop\newfile.avi'); % iter - 8, seq - 5
    video_writer.FrameRate = 5; video_writer.Quality = 100; open(video_writer);
end

%% Tracking
[settings, history] = set_batch_size(settings);
h = [];
for N = 1:settings.num_frames
    %disp(N);
    
    %% Quadratic-one
    if (settings.quadratic_one && N >= 3)
        
        X = X_init((B + T) * (N - 2) + 1:(B + T) * N);
        x_1 = history.x_batch(N - 1, 1:B + T)';
        x0 = history.x_batch(N - 1, (B + T) + 1:end)';
        
        [X, J, h] = sticks_finger_quadratic_one(X, x0, x_1, h, segments0, joints, frames, N, settings);
        H = h;
        H(1:B + T, 1:B + T) = diag(history.h_batch(N - 1, (B + T) * (settings.batch_size - 1) + 1:(B + T) * settings.batch_size));
    end
    
    %% Quadratic-two
    if (settings.quadratic_two && N >= 3)
        X = X_init((B + T) * (N - 2) + 1:(B + T) * N);
        
        x_1 = history.x_batch(N - 1, 1:B + T)';
        x0 = history.x_batch(N - 1, (B + T) + 1:end)';
        [X, J, h] = sticks_finger_laplace_approx(X, x0, x_1, h, segments0, joints, frames, N, settings);
        
        H = zeros(2 * (B + T), 2 * (B + T));
        a = h(1:B, 1:B); b = h(1:B, B + T + 1:B + T + B); c = h(B + T + 1:B + T + B, 1:B); d = h(B + T + 1:B + T + B, B + T + 1:B + T + B);
        
        if (settings.quadratic_two_marginalization)
            H(B + T + 1:B + T + B, B + T + 1:B + T + B) = d - c * inv(a) * b;
        end
        if (settings.quadratic_two_maximization)
            H(B + T + 1:B + T + B, B + T + 1:B + T + B) = d;
        end
        
        H(1:B + T, 1:B + T) = diag(history.h_batch(N - 1, (B + T) * (settings.batch_size - 1) + 1:(B + T) * settings.batch_size));
    end
    
    %% Kalman-like
    if (settings.kalman_like || (settings.kalman_two && N == 1))
        X = X_init((B + T) * (N - 1) + 1:(B + T) * N);
        if (N == 1)
            x0 = [];
            JtJ = zeros(B, B);
        else
            x0 = history.x_batch(N - 1, :)';
        end
        
        [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_kalman_like(X, x0, segments0, joints, frames{N}, JtJ, N, settings), X, settings.num_iters);
        
        if (settings.ground_truth_hessians)        
            J1 = sqrtm(theta_to_hessian_map(num2str(thetas_true(N, :))));
        else        
            J1 = J(1:end - B, 1:B);  
        end
        
        JtJ = JtJ + (J1'* J1);
        H = [JtJ, zeros(B, B); zeros(B, B), zeros(B, B)];
    end
    
    %% Kalman-two
    if (settings.kalman_two && N >= 2)
        
        X = X_init((B + T) * (N - 2) + 1:(B + T) * N);
        x0 = history.x_batch(N - 1, (B + T) + 1:end)';
        
        [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_kalman_two(X, x0, segments0, joints, frames{N}, JtJ, settings), X, settings.num_iters);
        
        if (settings.ground_truth_hessians)        
            J1 = sqrtm(theta_to_hessian_map(num2str(thetas_true(N, :))));
        else        
            J1 = J(1:settings.num_samples * 3, B + T + 1:B + T + B);  
        end
        
        H = diag([diag(JtJ); zeros(T, 1); diag(JtJ + (J1'* J1)); zeros(T, 1)]);
        JtJ = JtJ + (J1'* J1);
    end
    
    %% Separate optimization
    if (settings.independent)
        X = X_init((B + T) * (N - 1) + 1:(B + T) * N);
        [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_data(X, segments0, joints, frames{N}, settings, 'cpp'), X, settings.num_iters);
        
        % get hessian
        %[F, J, H] = sticks_finger_fg_data(X, segments0, joints, frames{N}, settings, 'numerical');
        %H = hessian_for_scalar_objective(F, J, H);
        
        H = J' * J;
   
    end
    
    %% Batch
    if (settings.batch || (N < 3 && (settings.quadratic_two || settings.quadratic_one)))
        if N <= settings.batch_size
            X = X_init(1:(B + T) * N);
            x0 = [];
        else
            X = X_init((B + T) * (N - settings.batch_size) + 1:(B + T) * N);
            x0 = history.x_batch(N - 1, 1:B + T)';
        end
        
        [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_batch(X, x0, segments0, joints, frames, N, settings), X, settings.num_iters);
        H = J' * J;
    end
    
    %% Batch simulation
    if (settings.batch_simulation)
        % independent data jacobian
        X = X_init((B + T) * (N - 1) + 1:(B + T) * N);
        [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_data(X, segments0, joints, frames{N}, settings, 'cpp'), X, settings.num_iters);
        H = J' * J;
        history.hessian_independent(N, :, :) = H(1:B, 1:B);
        history.mu_independent(N, :) = X(1:B);
        
        % batch simulation jacobian
        if N <= settings.batch_size
            X = X_init(1:(B + T) * N);
            x0 = [];
        else
            X = X_init((B + T) * (N - settings.batch_size) + 1:(B + T) * N);
            x0 = history.x_batch(N - 1, 1:B + T)';           
        end
        
        x_ = []; if (N > 1), x_ = history.x_batch(N - 1, end - B - T + 1:end)'; end
        
        [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_batch_simulation(X, x0, x_, segments0, joints, frames, N, settings, history), X, settings.num_iters);
        H = J' * J;
        
        %{
        if (N <= settings.batch_size)
           for M = 1:N
               history.hessian_independent(M, :, :) = H((B + T) * (M - 1) + 1:(B + T) * (M - 1) + B, (B + T) * (M - 1) + 1:(B + T) * (M - 1) + B);
               history.mu_independent(M, :) = X((B + T) * (M - 1) + 1:(B + T) * (M - 1) + B);
           end
        else
           for M = 1:settings.batch_size
               history.hessian_independent(N - settings.batch_size + M, :, :) = H((B + T) * (M - 1) + 1:(B + T) * (M - 1) + B, (B + T) * (M - 1) + 1:(B + T) * (M - 1) + B);
               history.mu_independent(N - settings.batch_size + M, :) = X((B + T) * (M - 1) + 1:(B + T) * (M - 1) + B);
           end
        end
        %}
        
    end
    
    
    %% Save new history
    if N <= settings.batch_size
        history.x_batch(N, :) = [zeros((B + T) * (settings.batch_size - N), 1); X];
        history.h_batch(N, :) = [zeros((B + T) * (settings.batch_size - N), 1); diag(H)];
    else
        history.x_batch(N, :) = X;
        history.h_batch(N, :) = diag(H);
    end
    
    history = get_last_frame_covariance(H, history, settings, N);
    
    %% Display
    display_jacobian(X, J, settings, N);
    
    if (settings.display_converged)
        beta = X(end - (B + T) + 1:end - T);
        theta = X(end - T + 1:end);
        data_points = frames{N};
        [segments0, joints] = segments_and_joints_2D();
        [segments0] = shape_2D(segments0, beta);
        [segments] = pose_2D(segments0, joints, theta);
        [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
        clf; hold on; axis off; axis equal; set(gcf,'color','w');
        display_sticks_finger(segments, data_points, model_points);
        drawnow; pause(0.05); %waitforbuttonpress
        
        if settings.write_video
            f = getframe();
            writeVideo(video_writer, f.cdata);
        end
        
    end
end

if exist('video_writer', 'var') && ~isempty(video_writer), video_writer.close(); end


