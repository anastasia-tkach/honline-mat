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
with_proper_lm = false;
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

%% Tracking
for N = 1:num_frames
    disp(N);
    E_previous = inf; lambda = 1; energy = []; beta = zeros(0, B);
    
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
    for iter = 1:num_iters
        j1 = cell(N, 1);
        f1 = cell(N, 1);
        
        for i = max(1, N - batch_size + 1):N
            data_points = frames{i};
            
            %% Initialize
            [segments0, joints] = segments_and_joints_2D();
            [segments0] = shape_2D(segments0, betas{i});
            [segments] = pose_2D(segments0, joints, thetas{i});
            
            %% Compute correspondences
            [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
            
            %% Display
            if (to_display && iter == num_iters && i == N)
                %rows = 5;  colums = 3;
                %column = rem(i - 1, rows);
                %row = colums - floor((i - 1)/rows) - 1;
                %h = subplot('Position', [column * (1/rows), row * (1/colums), 1/rows - 0.005, 1/colums - 0.005]); hold on;
                %set(gca,'XTick',[], 'YTick', [], 'XColor', [1, 1, 1], 'YColor', [1, 1, 1]);
                figure(1); clf; hold on; axis off; axis equal; set(gcf,'color','w');
                display_kalman_2D(segments, data_points, model_points, iter);
            end
            
            %% Compute Jacobians
            [~, j_theta] = jacobian_ik_2D(segments, joints, model_points, data_points, segment_indices);
            [f1_i, j_beta] = jacobian_shape_2D(segments, model_points, data_points, segment_indices);
            
            j1{i} = [j_beta, j_theta];
            f1{i} = f1_i;
        end
        
        %% Build data jacobian
        F1 = zeros(D * N, 1);
        J1 = zeros(D * N, N * (B + T));
        for i = max(1, N - batch_size + 1):N
            F1(D * (i - 1) + 1:D * i) = f1{i};
            J1(D * (i - 1) + 1:D * i, (B + T) * (i - 1) + 1:(B + T) * i) = j1{i};
        end
        
        %% Compute closeness
        F2 = zeros(2, 1);
        J2 = zeros(2, N * (B + T));
        count = 1;
        for i = max(1, N - batch_size):N - 1
            coefficients = [1; 1; 1];
            j = i + 1;
            F2(B * (count - 1) + 1: B * count) = diag(coefficients) * (betas{i} - betas{i + 1});
            if i > N - batch_size
                J2(B * (count - 1) + 1: B * count, (B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = - diag(coefficients) * eye(B, B);
            end
            J2(B * (count - 1) + 1: B * count, (B + T) * (j - 1) + 1:(B + T) * (j - 1) + B) = diag(coefficients) * eye(B, B);
            count = count + 1;
        end
        
        %% Solve
        I = zeros((B + T) * N, 1);
        indices = [zeros(B, 1); ones(T, 1)];
        indices = repmat(indices, N, 1);
        I(indices == 0) = 0.1;
        I(indices == 1) = 50;
        I = diag(I);
        w2 = 10;
        LHS = J1' * J1 + w2 * (J2' * J2);
        RHS = J1' * F1 + w2 * (J2' * F2);
        delta = (LHS + lambda * I) \ RHS;
        
        %% Update lambda
        E = (F1' * F1) + w2 * (F2' * F2);
        energy(iter) = E; beta(iter, 1:B) = betas{N};
        if (with_proper_lm)
            if (E > E_previous)
                lambda = lambda * 5;
                LHS = LHS_previous;
                RHS = RHS_previous;
                delta = (LHS + lambda * I) \ RHS;
                thetas = thetas_previous;
                betas = betas_previous;
                energy(iter) = energy(iter - 1); beta(iter, 1:B) = beta(iter - 1, 1:B);
            else
                lambda = lambda / 1.5;
                E_previous = E;
                LHS_previous = LHS;
                RHS_previous = RHS;
                betas_previous = betas;
                thetas_previous = thetas;
            end
        end
        
        %% Kalman
        if (run_kalman_filter)
            F1 = f1{N};
            J1 = j1{N};
            I = diag([0.1 * ones(B, 1); 50 * ones(T, 1)]);
            R = diag(r * ones(B + T, 1));
            dx = [betas{N}; thetas{N}] - [beta0; theta0];
            J2 = sqrt(I);  H = [J1; J2];
            F2 = zeros(B + T, 1); F = [F1; F2];
            LHS = H' * H + R * C_;
            %delta = LHS \ (J1' * F1);
            delta = -dx + LHS \ (H' * (F + H * dx));
            if (iter == num_iters)
                C = C_ + 1/r * (J1' * J1 + I);
                C_ = inv(inv(C) + Q);
                P = inv(C);
            end
            delta = [zeros((N - 1) * (B + T), 1); delta];
        end
        
        %% Update
        for i = 1:N
            betas{i} = betas{i} + delta((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B);
            thetas{i} = thetas{i} + delta((B + T) * (i - 1) + B + 1:(B + T) * i);
        end
        
    end
    history{N}.P = P;
    history{N}.betas = betas;
    history{N}.thetas = thetas;
    history{N}.JtJ = J1' * J1;
    history{N}.mean = betas{N};
end

%figure; hold on; plot(energy(2:end), 'lineWidth', 2);
%figure; hold on; plot(beta(:, 1), 'lineWidth', 2); plot(beta(:, 2), 'lineWidth', 2); plot(beta(:, 3), 'lineWidth', 2);

%% Display online equivalent
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
