% close all;
% clc; clear;
%rng default;

ylimit = [-1.8, -0.2]; X = []; num_iters = 50;
t_start = 0.01; t_end = 4.0;
settings.measurement_noise_std = 0.08;
x_true = -1;

x_init = 0.4 * randn;
% x_init = 0.3 + 0.2 * rand; v = randi([0, 1]);
% if randi([0, 1]) == 0, x_init = -1 * x_init; end
x_init = x_init + x_true;

f_ = @(x, t) exp(x * t)^2;
j_ = @(x, t) 2 * t * exp(x * t)^2;

%% Generate data
certain = t_start + 0.5;
uncertain = t_end - 0.2;
%T = linspace(t_start + 0.3, t_end - 0.1, 11);
%T = [linspace(t_end - 0.1, t_end - 0.3, 3)'; linspace(t_start + 0.4, t_start + 0.6, 4)'; linspace(t_end - 0.35, t_end - 0.15, 7)';]; % MAIN ROBLEM
T = [uncertain * ones(4, 1); certain * ones(4, 1); uncertain * ones(7, 1);];

num_data = length(T);
Y = zeros(num_data, 1);
for i = 1:num_data
    Y(i) = f_(x_true, T(i)) + settings.measurement_noise_std * randn();
    while Y(i) < 0,
        Y(i) = f_(x_true, T(i)) + settings.measurement_noise_std * randn();
    end
end

%% Display
T_plot = linspace(t_start, t_end);
F_plot = exp(x_true * T_plot).^2;
line_colors = {[0.3, 0.8, 1.0], [1, 0.6, 0.1]};
point_colors = {[0.7, 0.1, 0.6], [1, 0.4, 0.1]};

%% Optimize
display = false;

settings.laplace_approx = false;
settings.last_n = false;
settings.kalman_like = false;
settings.kalman = false;
settings.quadratic_all = false;
settings.batch = true;
settings.independent = false;

w2 = 1; w3 = 1;

settings.batch_independent = false;
settings.no_lm = false;
settings.batch_size = 3;
recompute_egh = true;

X = [];
for N = 1:num_data
    %% Quadratic two
    if (settings.laplace_approx)
        if N < 3
            X0 = x_init * ones(N, 1);
            [X, J] = my_lsqnonlin(@(X) simple_problem_fg_batch(X, Y, T, N, settings, w2), X0, num_iters);            
        else
            X_prev = [X; x_init];
            X0 = [X(1:N - 2); x_init * ones(2, 1)];
            %[xx_opt, J] = simple_problem_quadratic_two(X0, Y, T, N, w2);
            if N == 3, h = zeros(3, 2, 2); end
            [xx_opt, J, h] = simple_problem_laplace_approx(X0, X_prev, h, Y, T, N, w2, num_iters);            
            J = sqrtm(h); %disp([J' * J, h]);
            
            X = [X(1:N - 2); xx_opt];
            
            %% Draw covariance
            if display
                frame_certainty = T(N) < 1.5;
                draw_covariance_matrix(xx_opt, inv(h), frame_certainty); 
                xlim([-5, 10]); ylim([-5, 5]);
                mypoint([-1; -1] - xx_opt, [1, 0.7, 0], 50);
            end
                        
        end
        JtJ = J'* J;
    end
    
    %% Last n
    if (settings.last_n)
        y = Y(N); t = T(N);
        if N == 1
            JtJ = 0;
            history = {};
        end
        X0 = [X; x_init];
        [X, ~] = my_lsqnonlin(@(X) simple_problem_fg_last_n(X, X0, history, y, t, N, w2, w3, settings.batch_size), X0, num_iters);
        J1 = 2 * t * exp(X(N) * t)^2;
        
        %J2 = 1; if N == 1, J2 = 0; end
        %JtJ = J1' * J1 + J2' * J2;
        JtJ = JtJ + J1' * J1;
    end
    
    %% Kalman-like
    if (settings.kalman_like)
        x = x_init; y = Y(N); t = T(N);
        if (N == 1)
            x_prev = 0;
            JtJ = 0;
        else
            x_prev = X(N - 1);
        end
        %if T(N) < 1.0, w2 = 0.05; else w2 = 0.5; end
        [x, ~] = my_lsqnonlin(@(x) simple_problem_fg_kalman_like(x, x_prev, y, t, JtJ, N, w2), x, num_iters);
        
        J1 = 2 * t * exp(x * t)^2;
        JtJ = JtJ + (J1'* J1);        
        X = [X; x];
    end
    
    %% Quadratic all
    if (settings.quadratic_all && N > 1)
        settings.batch_size = num_data;
        options = optimoptions(@fminunc,'Algorithm','quasi-newton');
        x = x_init;
        X_init = X;
        
        E1 = history{N - 1}.E1; E2 = history{N - 1}.E2;
        dE1 = history{N - 1}.dE1; dE2 = history{N - 1}.dE2;
        ddE1 = history{N - 1}.ddE1; ddE2 = history{N - 1}.ddE2;
        E = E1 + E2; dE = dE1 + dE2; ddE = ddE1 + ddE2;
        
        K = [zeros(1, N - 2), 1];
        
        %% Optimize x
        %options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', 'Jacobian','on', 'OptimalityTolerance', 1e-10);
        %x_opt = lsqnonlin(@(x) simple_problem_fg_quadratic_all(x, X_init, Y, T, E, dE, ddE, N, w2), x, [], [], options);
        [x_opt, ~] = my_lsqnonlin(@(x) simple_problem_fg_quadratic_all(x, X_init, Y, T, E, dE, ddE, N, w2), x, num_iters);
        
        X_opt_before = - (2 * w2 * K' * K  +  ddE + ddE') \ (- 2 * w2 * x' * K + dE - 2 * X_init' * ddE)';
        X_opt = - (2 * w2 * K' * K  +  ddE + ddE') \ (- 2 * w2 * x_opt' * K + dE - 2 * X_init' * ddE)';
        
        X = [X_opt; x_opt];
        
        if (recompute_egh)
            [E1, dE1, ddE1] = simple_problem_egh_1_batch(X, Y, T, N, settings.batch_size, w2);
            [E2, dE2, ddE2] = simple_problem_egh_2_batch(X, Y, T, N, settings.batch_size, w2);
        else
            [E1, dE1, ddE1] = simple_problem_egh_1_cumulative(E1, dE1, ddE1, X, Y, T, N, w2);
            [E2, dE2, ddE2] = simple_problem_egh_2_cumulative(E2, dE2, ddE2, X, Y, T, N, w2);
        end
        
        %simple_problem_plot_gradient(X, Y, T, dE1, N);
        %simple_problem_debug_recompute(x, X_init, x_opt, X_opt, Y, T, history{N - 1}.E, history{N - 1}.dE, history{N - 1}.ddE, N, w2);
    end
    
    %% settings.batch
    if (settings.batch || N == 1)
        if N <= settings.batch_size
            X0 = x_init * ones(N, 1);
        else
            X0 = [X(1:N - settings.batch_size); x_init * ones(settings.batch_size, 1)];
        end
        %options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', 'Jacobian','on', 'OptimalityTolerance', 1e-10);
        %X = lsqnonlin(@(X) simple_problem_fg_batch(X, Y, T, N, w2), X0, [], [], options);
        
        [X, J] = my_lsqnonlin(@(X) simple_problem_fg_batch(X, Y, T, N, settings, w2), X0, num_iters);
    end
    
    %% Kalman
    if (settings.kalman)
        if (N == 1), x = x_init;
        else x = history{N - 1}.X(N - 1); end
        if N == 1, C_ = 0.001; end
        [x, C_] = simple_problem_kalman(x, Y(N), T(N), C_, num_iters);
        if (N == 1), X = x;
        else X = [history{N - 1}.X; x]; end
    end
    
    %% Display
    if 0; %N == num_data
        figure; hold on;
        plot(T_plot, F_plot, 'lineWidth', 3, 'color', [0.55, 0.75, 0.8]);
        point_colors = {[1, 0.6, 0.1],  [1, 0.8, 0]};
        xlim([t_start, t_end]);
        for i = 1:N
            mypoint([T(i), Y(i)], point_colors{1}, 50);
            plot(T_plot, exp(X(i) * T_plot).^2, 'lineWidth', 1, 'color', point_colors{1});
        end
        title(['x = ', num2str(X(i))]);
        ylim([-0.2, 1]);
    end
    
    %% Record history
    history{N}.X = X;
    if (settings.laplace_approx)
        if N < 3
            history{N}.JtJ = JtJ;
        else
            history{N}.JtJ = zeros(N, N); 
            history{N}.JtJ(1:N - 2, 1:N - 2) = history{N - 1}.JtJ(1:N - 2, 1:N - 2);
            history{N}.JtJ(N-1:N, N-1:N) = JtJ;
        end
    end
    if (settings.kalman_like || settings.last_n)
        history{N}.JtJ = zeros(N, N);
        if (N > 1), history{N}.JtJ(1:N-1, 1:N-1) = history{N - 1}.JtJ; end
        history{N}.JtJ(N, N) = JtJ;
        %history{N}.JtJ = JtJ';
    end
    if (settings.batch || settings.quadratic_all)
        [F, J] = simple_problem_fg_batch(X, Y, T, N, settings, w2);
        %J1 = J(1:N, :); history{N}.JtJ = J1' * J1;
        J1 = J(1:N, :); history{N}.JtJ = J' * J;
    end
    if settings.quadratic_all
        if N == 1
            history{N}.E1 = F' * F; history{N}.E2 = 0;
            history{N}.dE1 = 2 * F' * J; history{N}.dE2 = 0;
            history{N}.ddE1 = compute_hessian(X, Y, T, N, w2); history{N}.ddE2 = 0;
        else
            history{N}.E1 = E1; history{N}.E2 = E2;
            history{N}.dE1 = dE1; history{N}.dE2 = dE2;
            history{N}.ddE1 = ddE1; history{N}.ddE2 = ddE2;
        end
    end
    if (settings.kalman)
        if (N > 1), history{N}.JtJ(1:N-1, 1:N-1) = history{N - 1}.JtJ; end
        history{N}.JtJ(N, N) = 0.2 * 0.2 * C_;
    end
end

%% Display history
if display
    w = N * 0.040909;
    offset = 1/N;
    line_color = [1, 0.9, 0.3];
    point_color = [0.3, 0.6, 0.8];
    figure('units', 'normalized', 'outerposition', [0.1, 0.3, w, 0.55]); hold on;
    set(gca,'position', [0.05 0.08 0.93 0.90], 'units','normalized');
    for j = 1:N
        if (T(j) > 1.5)
            rectangle('Position',[j, ylimit(1), 1, ylimit(2) - ylimit(1)],'FaceColor',[1; 0.98; 0.95],'EdgeColor','none')
        else
            rectangle('Position',[j, ylimit(1), 1, ylimit(2) - ylimit(1)],'FaceColor',[0.96; 1; 0.93],'EdgeColor','none')
        end
    end
    plot(1:length(history) + 1, x_true * ones(length(history) + 1, 1), 'lineWidth', 2, 'color', [1, 0.5, 0.4], 'lineStyle', '-.');
    
    for j = 1:N
        myline([j, ylimit(1)], [j, ylimit(2)], [0.85, 0.85, 0.85], 1.5);
        for k = 1:j
            current_line_color = (1 - k / j ) * [0.85, 0.85, 0.85] + k / j * line_color;
            current_point_color = (1 - k / j ) * [0.75, 0.75, 0.75] + k / j * point_color;
            point_size = 20;
            line_width = 3.5;
            if k == j
                current_line_color = [1, 0.75, 0.3];
                current_point_color = [0.3, 0.7, 0.8];
                point_size = 40;
                line_width = 4.5;
            end
            importance = history{j}.JtJ(k, k);
            myline([j + offset * k, history{j}.X(k) + 0.1 * sqrt(importance)], ...
                [j + offset * k, history{j}.X(k) - 0.1 * sqrt(importance)], current_line_color, line_width);
            mypoint([j + offset * k, history{j}.X(k)], current_point_color, point_size);
        end
    end
    set(gca, 'fontSize', 13); xlim([1, length(history) + 1]); ylim(ylimit);
end


