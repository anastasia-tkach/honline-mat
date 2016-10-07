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
N = length(thetas_true);
[frames, beta_init, thetas_init] = get_random_data_from_theta(beta_true, thetas_true, ...
    settings.beta_noise_std, settings.theta_noise_std, settings.measurement_noise_std, num_samples);

num_frames = length(thetas_true);
num_iters = 20;
to_display = false;

init_from_previous_frame = false;

if (to_display), figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.8, 0.8]);
    axis off; axis equal; hold on;
end

settings.quadratic_two = false;
settings.last_n = false;
settings.kalman_like = false;
settings.kalman = false;
settings.quadratic_all = false;
settings.batch = true;
settings.independent = false;
settings.no_lm = false;

settings.batch_size = 1;

w2 = 0;

%% Tracking
for N = 1:num_frames
    %disp(N);
    
    %if (to_display), figure('units', 'normalized', 'outerposition', [0.1, 0.3, 0.8, 0.5]); axis off; axis equal; hold on; end
    
    %% Separate optimization
    if (settings.independent)
        X = zeros(N * (B + T), 1);
        J = zeros(N * (B + T), N * (B + T));
        for i = 1:N
            x = [beta_init; thetas_init{i}];
            [x, j] = my_lsqnonlin(@(x) sticks_finger_fg_single(x, segments0, joints, frames{i}), x, num_iters);
            betas{i} = x(1:B);
            thetas{i} = x(B + 1:B + T);
            X((B + T) * (i - 1) + 1:(B + T) * i) = x;
            J(D * (i - 1) + 1:D * i, (B + T) * (i - 1) + 1:(B + T) * i) = j;
        end
    end
    
    %% Batch optimization
    if (settings.batch)
        X = zeros(N * (B + T), 1);
        for i = 1:max(1, N - settings.batch_size + 1) - 1
            X((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = betas{i};
            X((B + T) * (i - 1) + B + 1:(B + T) * i) = thetas{i};
        end
        for i = max(1, N - settings.batch_size + 1):N
            X((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = beta_init;
            X((B + T) * (i - 1) + B + 1:(B + T) * i) = thetas_init{i};
        end
        
        if (settings.no_lm)
            lambdas = zeros((B + T) * N, 1);
            indices = [zeros(B, 1); ones(T, 1)]; indices = repmat(indices, N, 1);
            lambdas(indices == 0) = 0.1; lambdas(indices == 1) = 50;
            for iter = 1:num_iters
                [F, J] = sticks_finger_fg_batch(X, segments0, joints, frames, N, D, settings.batch_size, w2);
                %delta = - (J' * J + lambda * eye(size(J, 2), size(J, 2))) \ (J' * F);
                delta = - (J' * J + diag(lambdas)) \ (J' * F);
                X = X + delta;
            end
        else
            [X, J] = my_lsqnonlin(@(X) sticks_finger_fg_batch(X, segments0, joints, frames, N, D, settings.batch_size, w2), X, num_iters);
        end
        
        %betas = cell(N, 1);
        %thetas = cell(N, 1);
        for i = max(1, N - settings.batch_size + 1):N
            betas{i} = X((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B);
            thetas{i} = X((B + T) * (i - 1) + B + 1:(B + T) * i);
        end
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
    history{N}.betas = betas;
    history{N}.thetas = thetas;
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