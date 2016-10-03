%close all;
%clear; clc;
rng(457645);

%% Parameters
num_samples = 15;
B = 3;
T = 3;
D = 3 * num_samples;
measurement_noise_std = 0.07;  % 0.01 - converges to true beta, 0.1 - does not find true beta
beta_noise_std = 0.5; %3
theta_noise_std = 0.15;
blocks = {[1, 2], [2, 3], [3, 4]};
beta_true = [3; 3; 3];

R = 0.2 * eye(D + B + T, D + B + T);
Q = diag([0 * [1; 1; 0]; theta_noise_std * ones(T, 1)]); %0.0001
P = diag([beta_noise_std * ones(B, 1); zeros(T, 1)]);

replay = true;

%% Scrip
%%{
beat = 5;
%thetas_true_2 = [linspace(0, 0, beat), linspace(0, pi/3, beat), linspace(pi/3, 0, beat), linspace(0, 0, beat),    linspace(0, 0, beat),     linspace(0, pi/3, beat), linspace(pi/3, 0, beat), linspace(0, 0, 2 * beat)];
%thetas_true_3 = [linspace(0, 0, beat), linspace(0, 0, beat),    linspace(0, 0, beat),    linspace(0, pi/3, beat), linspace(pi/3, 0, beat),  linspace(0, pi/3, beat), linspace(pi/3, 0, beat), linspace(0, 0, 2 * beat)];
%thetas_true = [zeros(length(thetas_true_2), 1), thetas_true_2', thetas_true_3'];

thetas_true = [
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    pi/3, pi/3, pi/3;
    pi/3, pi/3, pi/3;
    pi/3, pi/3, pi/3;
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    pi/3, pi/3, pi/3;
    pi/3, pi/3, pi/3;
    pi/3, pi/3, pi/3;
    0, 0, 0;
    0, 0, 0;
    0, 0, 0;
    ];

N = length(thetas_true);
num_frames = length(thetas_true);

thetas_init = cell(num_frames, 1);
beta_init = beta_true + beta_noise_std * randn(B, 1);
theta_init = [0; 0; 0];
thetas = cell(num_frames, 1);
betas = cell(num_frames, 1);
frames = cell(length(thetas_true), 1);
for i = 1:length(thetas_true)
    thetas_init{i} = thetas_true(i, :)' + theta_noise_std * randn(T, 1);
    thetas{i} = thetas_init{i};
    betas{i} = beta_init;
    [segments0, joints] = segments_and_joints_2D();
    [data_segments] = shape_2D(segments0, beta_true);
    [data_segments] = pose_2D(data_segments, joints, thetas_true(i, :)');
    [data_points] = sample_2D(data_segments, num_samples);
    for j = 1:length(data_points)
        data_points{j} = data_points{j} + measurement_noise_std * randn(2, 1);
    end
    frames{i} = data_points;
end
%%}
%{
%dataset_entry_number = 50;
dataset_path = 'C:\Users\t-antka\OneDrive - Microsoft\Data\CalibrationDataset\';
name_suffix = sprintf('%03d', dataset_entry_number);
load([dataset_path, 'frames_', name_suffix]);
load([dataset_path, 'thetas_true_', name_suffix]);
load([dataset_path, 'beta_init_', name_suffix]);
theta_init = [0; 0; 0];
%}
num_iters = 15;
num_frames = length(thetas_true);
start_frame = 1;
to_display = false;
with_covariance = false;
with_many_previous_frames = false;
init_from_previous_frame = false;
with_latent_variables = false; beta_latent = beta_init;
certainty_threshold = 5;
history = {};
%% Tracking
for N = 1:num_frames
    disp(N);
    for i = start_frame:N
        
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
    if (to_display), figure('units', 'normalized', 'outerposition', [0.1, 0.3, 0.8, 0.5]); axis off; axis equal; hold on; end
    for iter = 1:num_iters
        j1 = cell(N, 1);
        f1 = cell(N, 1);
        for i = 1:N
            data_points = frames{i};
            
            %% Initialize all frames
            [segments0, joints] = segments_and_joints_2D();
            [segments0] = shape_2D(segments0, betas{i});
            [segments] = pose_2D(segments0, joints, thetas{i});
            
            %% Compute correspondences
            [segment_indices, model_points] = compute_correspondences_2D(segments, blocks, data_points);
            
            %% Display
            if (to_display && iter == num_iters)
                rows = 5;  colums = 3;
                column = rem(i - 1, rows);
                row = colums - floor((i - 1)/rows) - 1;
                h = subplot('Position', [column * (1/rows), row * (1/colums), 1/rows - 0.005, 1/colums - 0.005]); hold on;
                set(gca,'XTick',[], 'YTick', [], 'XColor', [1, 1, 1], 'YColor', [1, 1, 1]);
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
        for i = 1:N
            F1(D * (i - 1) + 1:D * i) = f1{i};
            J1(D * (i - 1) + 1:D * i, (B + T) * (i - 1) + 1:(B + T) * i) = j1{i};
        end
        
        
        if (~with_latent_variables)
            %% Compute closeness
            F2 = zeros(2, 1);
            J2 = zeros(2, N * (B + T));
            count = 1;
            for i = 1:N-1
                if (with_covariance)
                    coefficients = diag(history{N - 1}.JtJ((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B, (B + T) * (i - 1) + 1:(B + T) * (i - 1) + B));
                else
                    coefficients = [1; 1; 1];
                end
                if with_many_previous_frames && (coefficients(1) > certainty_threshold && coefficients(2) > certainty_threshold)
                    for j = i + 1:N
                        F2(B * (count - 1) + 1: B * count) = diag(coefficients) * (betas{i} - betas{i + 1});
                        J2(B * (count - 1) + 1: B * count, (B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = - diag(coefficients) * eye(B, B);
                        J2(B * (count - 1) + 1: B * count, (B + T) * (j - 1) + 1:(B + T) * (j - 1) + B) = diag(coefficients) * eye(B, B);
                        count = count + 1;
                    end
                else
                    j = i + 1;
                    F2(B * (count - 1) + 1: B * count) = diag(coefficients) * (betas{i} - betas{i + 1});
                    J2(B * (count - 1) + 1: B * count, (B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = - diag(coefficients) * eye(B, B);
                    J2(B * (count - 1) + 1: B * count, (B + T) * (j - 1) + 1:(B + T) * (j - 1) + B) = diag(coefficients) * eye(B, B);
                    count = count + 1;
                end
            end
            
            %% Solve
            I = zeros((B + T) * N, 1);
            indices = [zeros(B, 1); ones(T, 1)];
            indices = repmat(indices, N, 1);
            I(indices == 0) = 0.1;
            I(indices == 1) = 50;
            I = diag(I);
            w2 = 0;
            LHS = J1' * J1 + w2 * (J2' * J2) + I;
            RHS = J1' * F1 + w2 * (J2' * F2);
            delta = LHS \ RHS;
            
            %% Update
            for i = 1:N
                betas{i} = betas{i} + delta((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B);
                thetas{i} = thetas{i} + delta((B + T) * (i - 1) + B + 1:(B + T) * i);
            end
            
        end
        
        %% With latent variables
        if (with_latent_variables)
            F2 = zeros(N * (B - 1), 1);
            J2 = zeros(N * (B - 1), N * (B + T) + 2);
            count = 1;
            for i = 1:N
                F2(2 * (i - 1) + 1: 2 * i) = betas{i}(1:2) - beta_latent(1:2);
                J2(2 * (i - 1) + 1: 2 * i, (B + T) * (i - 1) + 1:(B + T) * (i - 1) + 2) = - eye(2, 2);
                J2(2 * (i - 1) + 1: 2 * i,  N * (B + T) + 1: N * (B + T) + 2) = eye(2, 2);
            end
            
            I = zeros((B + T) * N + 2, 1);
            indices = [zeros(B, 1); ones(T, 1)];
            indices = repmat(indices, N, 1);
            I(indices == 0) = 0.1;
            I(indices == 1) = 50;
            I(end - 1:end) = 0.1;
            I = diag(I);
            w2 = 0;
            J1 = [J1, zeros(size(J1, 1), 2)];
            LHS = J1' * J1 + w2 * (J2' * J2) + I;
            RHS = J1' * F1 + w2 * (J2' * F2);
            delta = LHS \ RHS;
            
            for i = 1:N
                betas{i} = betas{i} + delta((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B);
                thetas{i} = thetas{i} + delta((B + T) * (i - 1) + B + 1:(B + T) * i);
            end
            beta_latent(1:2) = beta_latent(1:2) + delta( N * (B + T) + 1: N * (B + T) + 2);
        end
        
        
    end
    history{N}.betas = betas;
    history{N}.thetas = thetas;
    history{N}.JtJ = J1' * J1;
    history{N}.mean = betas{N};
end

%% Display online equivalent
figure_borders = [0.05 0.08 0.93 0.90];
means = zeros(N, B);
vars = zeros(N, B);
trues = zeros(N, B);
%beta_std = zeros(N, B);
for i = 1:length(history)
    means(i, :) = history{i}.betas{i};
    trues(i, :) = beta_true;
    vars(i, :) = diag(history{i}.JtJ((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B, (B + T) * (i - 1) + 1:(B + T) * (i - 1) + B));
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
    
    %plot(1:length(history), means(:, i) + vars(:, i).^0.5, 'lineWidth', 2, 'color', [0.55, 0.75, 0.8]);
    %plot(1:length(history), means(:, i) - vars(:, i).^0.5, 'lineWidth', 2, 'color', [0.55, 0.75, 0.8]);
    
    %plot(1:length(history), trues(:, i) + 0.5 * beta_variance_threshold * ones(length(history), 1), 'lineWidth', 1, 'color', [0.7, 0.9, 0.5]);
    %plot(1:length(history), trues(:, i) + 0.5 * beta_std(:, i).^2, 'lineWidth', 1, 'color', [0.4, 0.7, 0.5]);
    
    ylim([2, 4]); xlim([0, length(history)]);
    set(gca, 'fontSize', 13); title(['beta ', num2str(i)]);
end

if exist('video_writer', 'var'), video_writer.close(); end

%% Display history
%{
offset = 0.07;
figure('units', 'normalized', 'outerposition', [0.1, 0.2, 0.55, 0.7]); hold on;
set(gca,'position', [0.05 0.08 0.93 0.90], 'units','normalized');
for b = 1:B - 1
    h = subplot(B - 1, 1, b); hold on;
    p = get(h, 'pos'); set(h, 'pos', [0.05, p(2) - 0.06, 0.9, 0.43]);
    set(gca,'XTick',[]);
    plot(1:length(history) + 1, beta_true(b) * ones(length(history) + 1, 1), 'lineWidth', 2, 'color', [1, 0.5, 0.4], 'lineStyle', '-.');
    for j = start_frame:N
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