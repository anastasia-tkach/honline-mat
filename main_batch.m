clear; clc;
rng default;

%% Parameters
settings.num_samples = 10;
B = 3 ;
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
theta_semicertain = [0, pi/30, pi/30];
theta_uncertain = [0, 0, 0];
tact = 3;
thetas_true = [repmat(theta_uncertain, tact, 1); repmat(theta_certain_1, tact, 1); repmat(theta_uncertain, tact, 1); repmat(theta_certain_2, tact, 1);  repmat(theta_uncertain, tact, 1)];
%thetas_true = [repmat(theta_uncertain, 4, 1); repmat(theta_certain_12, 4, 1); repmat(theta_uncertain, 7, 1)];
settings.num_frames = size(thetas_true, 1);

[frames, beta_init, thetas_init] = get_random_data_from_theta(beta_true, thetas_true, settings);
X_init = zeros((B + T) * settings.num_frames, 1);
for i = 1:settings.num_frames
    X_init((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = beta_init;
    X_init((B + T) * (i - 1) + B + 1:(B + T) * i) = thetas_init{i};
end

settings.store_covariance = true;
settings.display = false;
if (settings.display), 
    figure('units', 'normalized', 'outerposition', [0.25, 0.275, 0.45, 0.7]);
    axis off; axis equal; hold on;
end

settings.quadratic_one = false;
settings.quadratic_two = false;
settings.kalman_like = false;
settings.batch = true;
settings.independent = false;

settings.batch_size = settings.num_frames;
settings.num_iters = 20;

settings.batch_independent = false;
settings.batch_online = true;
settings.batch_online_robust = false;
settings.batch_online_robust_tau = 1;

settings.shape_prior = false;
settings.data_model_energy = true;
settings.model_data_energy = true;
settings.silhouette_energy = false;

[settings, history] = set_batch_size(settings);

settings.w2 = 1;
settings.w4 = 1;
h = [];

%% Tracking
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
    
    %% Laplace approximation
    if (settings.quadratic_two && N >= 3)
        X = X_init((B + T) * (N - 2) + 1:(B + T) * N);
        
        x_1 = history.x_batch(N - 1, 1:B + T)';
        x0 = history.x_batch(N - 1, (B + T) + 1:end)';
        [X, J, h] = sticks_finger_laplace_approx(X, x0, x_1, h, segments0, joints, frames, N, settings);
        
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
        
        [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_kalman_like(X, x0, segments0, joints, frames{N}, JtJ, N, settings), X, settings.num_iters);
        
        J1 = J(1:end - B, 1:B);
        JtJ = JtJ + (J1'* J1);
        H = [JtJ, zeros(B, B); zeros(B, B), zeros(B, B)];
    end
    
    %% Separate optimization
    if (settings.independent)
        X = X_init((B + T) * (N - 1) + 1:(B + T) * N);
        [X, j] = my_lsqnonlin(@(X) sticks_finger_fg_data(X, segments0, joints, frames{N}, settings), X, settings.num_iters);
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
        
        [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_batch(X, x0, segments0, joints, frames, N, settings), X, settings.num_iters);
        H = J' * J;
    end
    
    %% Save new history
    if N <= settings.batch_size
        history.x_batch(N, :) = [zeros((B + T) * (settings.batch_size - N), 1); X];
        history.h_batch(N, :) = [zeros((B + T) * (settings.batch_size - N), 1); diag(H)];
    else
        indices = 1:(B + T) *  settings.batch_size;
        history.x_batch(N, :) = X(indices);
        history.h_batch(N, :) = diag(H);
    end
    if settings.store_covariance
        history.covariance(N, :, :) = H(end - B - T + 1:end - T,  end - B - T + 1:end - T);
    end
    
    %% Covariance   
    %draw_covariance_matrix(X(end - B - T  + 1:end - B - T  + 2), inv(H(end - B - T  + 1:end - B - T  + 2, end - B - T  + 1:end - B - T  + 2)), 0);
    %xlim([-4, 4]); ylim([-4, 4]); title(['frame ', num2str(N)]); set(gca, 'fontSize', 12); set(gca,'fontname','Cambria');
    
    %% Display
    if (settings.display)
        beta = X(end - (B + T) + 1:end - T);        
        theta = X(end - T + 1:end);
        data_points = frames{N};
        [segments0, joints] = segments_and_joints_2D();
        [segments0] = shape_2D(segments0, beta);
        [segments] = pose_2D(segments0, joints, theta);
        [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
        clf; hold on; axis off; axis equal; set(gcf,'color','w');
        display_sticks_finger(segments, data_points, model_points);
        pause(0.1);
        %waitforbuttonpress
    end
end



