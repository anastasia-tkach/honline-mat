%close all;
clear; clc;
rng(2572);
load betas_prior_std; load thetas_init; thetas_prior_init = thetas_init;
dataset_path = 'C:\Users\t-antka\OneDrive - Microsoft\Data\CalibrationDataset\';

%% Parameters
num_samples = 15;
B = 3;
T = 3;
D = 3 * num_samples;
measurement_noise_std = 0.07;
beta_noise_std = 0.5;
theta_noise_std = 0.15;
blocks = {[1, 2], [2, 3], [3, 4]};
beta_variance_threshold = 0.15;

r = 0.2;
R = r * eye(D + B + T, D + B + T);
Q = diag([0 * [1; 1; 0]; theta_noise_std * ones(T, 1)]);
P = diag([beta_noise_std * ones(B, 1); zeros(T, 1)]);

%% Algorithm type
to_display = false;
to_debug = false;
init_from_previous_frame = false;
get_random_sequence = false;
get_scripted_sequence = true;
get_dataset_sequence = false;

writing_video = false;
if (writing_video)
    video_writer = VideoWriter('C:\Users\t-antka\Desktop\newfile.avi');
    video_writer.FrameRate = 7; video_writer.Quality = 100; open(video_writer);
end

num_iters = 20;

with_pose_certainty = false;
mimicking_objective_function = false;
nothing = false;
kalman_filter = true;
extended_kalman_filter = false;

information_filter = false;
extended_information_filter = false;

pose_prior_for_system_noise = false;

if (kalman_filter), num_iters = 1; end
if (information_filter), num_iters = 1; end
P_ = P + Q;
C_ = inv(P_);


%% Get data
beta_true = [3; 3; 3];
theta_init = [0; 0; 0];

if (get_dataset_sequence)
    beta_init = beta_true + beta_noise_std * randn(B, 1);
    theta_init = [0; 0; 0];
    theta_true = theta_init + theta_noise_std * randn(T, 1);
    dataset_entry_number = 50;
    name_suffix = sprintf('%03d', dataset_entry_number);
    load([dataset_path, 'frames_', name_suffix]);
    load([dataset_path, 'thetas_true_', name_suffix]);
    load([dataset_path, 'beta_init_', name_suffix]);
end

if (get_scripted_sequence)
    thetas_true = [0, 0, 0; 0, 0, 0; 0, 0, 0;
        0, pi/4, pi/4; 0, pi/4, pi/4; 0, pi/4, pi/4;
        0, 0, 0; 0, 0, 0; 0, 0, 0;
        0, pi/4, pi/4; 0, pi/4, pi/4; 0, pi/4, pi/4;
        0, 0, 0; 0, 0, 0; 0, 0, 0;];
    thetas_true = [linspace(0, 0, 17); 0, 0, 0, linspace(0, pi/4, 5), pi/4, linspace(pi/4, 0, 5), 0, 0, 0; ...
        0, 0, 0, linspace(0, pi/4, 5), pi/4, linspace(pi/4, 0, 5), 0, 0, 0]'; 
    [frames, beta_init, thetas_init] = get_random_data_from_theta(beta_true, thetas_true, ...
        beta_noise_std, theta_noise_std, measurement_noise_std, num_samples);
end

N = length(frames);

%% Initialization
if (init_from_previous_frame)
    beta = beta_init;
    theta = theta_init;
    [segments0, joints] = segments_and_joints_2D();
    [segments0] = shape_2D(segments0, beta);
    [segments] = pose_2D(segments0, joints, theta);
end

%% Tracking
if (to_display), figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.8, 0.8]);
    axis off; axis equal; hold on;
end
if (to_debug),
    figure('units', 'normalized', 'outerposition', [0.1, 0.2, 0.33, 0.7]); hold on;
    set(gca,'position', [0.05 0.08 0.93 0.90], 'units','normalized');
end
history = {};
for n = 1:N
    %disp(beta');
    %% Get data
    if (get_random_sequence)
        theta_true = theta_true + theta_noise_std * randn(size(theta));
        [data_segments] = shape_2D(segments0, beta_true);
        [data_segments] = pose_2D(data_segments, joints, theta_true);
        [data_points] = sample_2D(data_segments, num_samples);
        for i = 1:length(data_points)
            data_points{i} = data_points{i} + measurement_noise_std * randn(2, 1);
        end
    end
    if (get_scripted_sequence && ~init_from_previous_frame)
        theta = thetas_init{n};
        beta = beta_init;
        [segments0, joints] = segments_and_joints_2D();
        [segments0] = shape_2D(segments0, beta);
        [segments] = pose_2D(segments0, joints, theta);
    end
    
    theta_true = thetas_true(n, :)';
    data_points = frames{n};
    beta0 = beta;
    theta0 = theta;
    
    for iter = 1:num_iters
        dx = [beta; theta] - [beta0; theta0];
        %disp(beta');
        %% Compute correspondences
        [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
        
        %% Display
        if (to_debug && iter == num_iters)
            figure(2); clf; hold on; axis off; axis equal; set(gcf,'color','w');
            display_history(history, beta_init, beta_variance_threshold, N, B, T);
        end
        if (to_display && iter == num_iters)
            figure(1); clf; hold on; axis off; axis equal; set(gcf,'color','w');
            if (~nothing)
                display_posed_model(beta, theta, 5 * (ones(B + 1, 1) + 5 * [diag(P(1:B, 1:B)); 0]), [1.0, 0.87, 0.8]);
            end
            display_kalman_2D(segments, data_points, model_points, iter);
            
            if (writing_video)
                f = getframe();
                writeVideo(video_writer, f.cdata);
            end
            %waitforbuttonpress;
            %break;
        end
        
        %% Compute Jacobians
        [~, J_theta] = jacobian_ik_2D(segments, joints, model_points, data_points, segment_indices);
        [F_, J_beta] = jacobian_shape_2D(segments, model_points, data_points, segment_indices);
        
        J = [J_beta, J_theta];
        I = 0.1 * eye(T + B, T + B);
        I(end - T + 1:end, end -T + 1:end) = 500 * I(end - T + 1:end, end - T + 1:end);
        
        J1 = J; J2 = sqrt(I);
        F1 = F_; F2 = zeros(B + T, 1);
        H = [J1; J2];  F = [F1; F2];
        %disp(F' * F);
        
        %% Filtering
        if (kalman_filter)
            P_ = P + Q;
            K = (P_ * H') / (H * P_ * H' + R );
            P = (eye(B + T, B + T) - K * H) * P_;
            delta = K * F;
        end
        if (information_filter)
            P_ = P + Q;
            P = inv(inv(P_) + H' * inv(R) * H);
            K = P * H' * inv(R);
            delta = K * F;
        end
        if (extended_kalman_filter)
            K = (P_ * H') / (H * P_ * H' + R );
            delta = K * F + (eye(B + T, B + T) - K * H) * ([beta0; theta0] - [beta; theta]);
        end
        if (extended_information_filter)
            P = inv(inv(P_) + H' * inv(R) * H);
            K = P * H' * inv(R);
            delta = K * F + (eye(B + T, B + T) - K * H) * ([beta0; theta0] - [beta; theta]);
        end
        if (mimicking_objective_function)
            if (with_pose_certainty)
                beta_variance = squeeze(prior_lookup(theta0, betas_prior_std, thetas_init)).^2;
                indicator = beta_variance(1:B - 1) < beta_variance_threshold;
                d =  - 0.7 * indicator .* diag(r * C_(1:B - 1, 1:B - 1));
                w = r * ones(B + T, 1);
                h = [0; 0]; if (n > 1), h = history{n - 1}.energy; end
                w(1) = 0.02 + r * beta_variance(1)^0.5 + 0 * r * h(1);
                w(2) = 0.02 + r * beta_variance(2)^0.5 + 0 * r * h(2);
                R = diag(w);
            else
                R = diag(r * ones(B + T, 1));
            end
            LHS = J' * J + I + R * C_;
            %delta = LHS \ (J' * F_);
            delta = -dx + LHS \ (H' * (F + H * dx));
            
            if (iter == num_iters)
                C = C_ + 1/r * (J' * J + I);
                %C = C_ + 1/r * diag(diag(J' * J + I));
                %C = diag(max(diag(C_), 1/r * diag(J' * J + I)));
                %C = C_ + 0.1 /r * diag(diag(J' * J + I));
                %C = C_ + 0.1 * inv(R) * diag(diag(J' * J + I));
                
                C_ = inv(inv(C) + Q);
                P = inv(C);
            end
            
        end
        if (nothing)
            D = diag([0 * ones(B, 1); zeros(T, 1)]);
            LHS = J' * J + I + D;
            delta = LHS \ (J' * F_);
            %delta = LHS \ (J' * F_) + (LHS \ (J' * J + I) - eye(B + T, B + T)) * ([beta; theta] - [beta0; theta0]);
            P = 1/r * diag(diag(D).^(-1));
        end
        
        %% Update
        beta = beta + delta(1:B);
        theta = theta + delta(B + 1:end);
        [shaped_segments] = shape_2D(segments0, beta);
        [segments] = pose_2D(shaped_segments, joints, theta);
        
        %% Stopping condition
        %disp(F' * F);
        %if F' * F < 0.1, break; end
    end
    
    %% extended_kalman_filter Filter - measurement and time update
    if (pose_prior_for_system_noise)
        Q(1:B, 1:B) = 0.1 * diag(squeeze(prior_lookup(theta, betas_prior_std, thetas_prior_init)));
    end
    if (extended_kalman_filter)
        P = (eye(B + T, B + T) - K * H) * P_;
        P_ = P + Q;
    end
    if (extended_information_filter)
        P_ = P + Q;
    end
    
    %% History
    h = length(history);
    history{h + 1}.P = P;
    history{h + 1}.mean = [beta; theta];
    history{h + 1}.true = [beta_true; theta_true];
    history{h + 1}.beta_std = squeeze(prior_lookup(theta0, betas_prior_std, thetas_prior_init));
    history{h + 1}.energy(1) = F(segment_indices == 1)' * F(segment_indices == 1);
    history{h + 1}.energy(2) = F(segment_indices == 2)' * F(segment_indices == 2);
    
    if (to_debug), waitforbuttonpress; end
    %close all;
    %pause(0.05);
    %waitforbuttonpress;
end

%return
%close;

%% Display history
if (~to_debug && (~exist('dataset_test', 'var') || dataset_test == true)),
    figure('units', 'normalized', 'outerposition', [0.1, 0.2, 0.33, 0.7]); hold on;
    set(gca,'position', [0.05 0.08 0.93 0.90], 'units','normalized');
    display_history(history, beta_init, beta_variance_threshold, N, B, T);
end

if exist('video_writer', 'var'), video_writer.close(); end
