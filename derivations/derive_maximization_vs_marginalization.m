%close all;
clc; clear;
%rng(2346);

t_start = 0.1;
t_end = 5;
measurement_noise_std = 0.1;
x_true = -1;

x_init = 0.3 + 0.2 * rand;
v = randi([0, 1]);
if randi([0, 1]) == 0, x_init = -1 * x_init; end
x_init = x_init + x_true;

f_ = @(x, t) exp(x * t)^2;
j_ = @(x, t) 2 * t * exp(x * t)^2;

%% Generate data
%T = linspace(t_start, t_end, num_data);
%T = [linspace(5, 4, 3)'; linspace(1.5, 0.5, 4)'; linspace(3.75, 4.75, 3)'];
T = [1.5; 1.0; 1.2; 1.3];
num_data = length(T);
Y = zeros(num_data, 1);
for i = 1:num_data
    Y(i) = f_(x_true, T(i)) + measurement_noise_std * randn();
end

%% Display
T_plot = linspace(0, t_end);
F_plot = exp(x_true * T_plot).^2;
line_colors = {[0.3, 0.8, 1.0], [1, 0.6, 0.1]};
point_colors = {[0.7, 0.1, 0.6], [1, 0.4, 0.1]};

%%{
figure; hold on;
X = linspace(-4, -0.2, 15);
for i = 1:length(X)
    x = X(i);
    plot(T_plot, exp(x * T_plot).^2, 'lineWidth', 1.5);
end
ylim([-0.2, 1]); xlim([0, 4]);
% plot(T_plot, F_plot, 'lineWidth', 3);
% xlim([0, t_end]);
% for i = 1:num_data
%     mypoint([T(i), Y(i)], point_colors{2}, 50);
% end
%%}
return
%% Optimize
num_iter = 25;
for N = 1:num_data
    X = x_init * ones(N, 1) + 0.01 * randn(N, 1);
    for iter = 1:num_iter
        
        %% Data term
        F1 = zeros(N, 1);
        J1 = zeros(N, N);
        for i = 1:N
            F1(i) = f_(X(i), T(i)) - Y(i);
            J1(i, i) = - j_(X(i), T(i));
        end
        
        %% Closeness term
        F2 = zeros(N - 1, 1);
        J2 = zeros(N - 1, N);
        for i = 1:N - 1
            F2(i) = X(i) - X(i + 1);
            J2(i, i) = -1;
            J2(i, i + 1) = 1;
        end
        
        %% Solve
        w2 = 0;
        delta = (J1' * J1 + w2 * J2' * J2 + 0.2 * eye(N, N)) \ (J1' * F1 + w2 * J2' * F2);
        X = X + delta;
        %disp([x, f' * f]);
    end
    %{
    figure; hold on;
    plot(T_plot, F_plot, 'lineWidth', 3);
    xlim([0, t_end]);
    for i = 1:N
        mypoint([T(i), Y(i)], point_colors{2}, 50);
        plot(T_plot, exp(X(i) * T_plot).^2, 'lineWidth', 1, 'color', [0.55, 0.75, 0.8]);
    end
    %}
    history{N}.X = X;
    history{N}.JtJ = J1' * J1;
end

%%{
%% Optimize online
for N = 1:num_data
    X = x_init * ones(N, 1) + 0.02 * randn(N, 1);
    for iter = 1:num_iter + 1
        
        %% Data term
        F1 = zeros(N, 1);
        J1 = zeros(N, N);
        for i = 1:N
            F1(i) = f_(X(i), T(i)) - Y(i);
            J1(i, i) = j_(X(i), T(i));
        end
        
        %% Closeness term
        F2 = zeros(N - 1, 1);
        J2 = zeros(N - 1, N);
        for i = 1:N - 1
            F2(i) = X(i) - X(i + 1);
            J2(i, i) = -1;
            J2(i, i + 1) = 1;
        end
        
        %% Recursion
        if N == 4
            w2 = 0.1;
            
            J1_prev = history{N - 1}.J1;
            F1_before = F1;
            
            F1 = @(X, N) exp(diag(T(1:N)) * X(1:N)).^2 - Y(1:N);
            dF1 = @(X, N) diag(2 * diag(T(1:N)) * exp(diag(T(1:N)) * X(1:N)).^2);
            for i = 1:N, ddF1{i} = @(X, N) diag([zeros(i - 1, 1); 4 * T(i).^2 * exp(diag(T(i)) * X(i)).^2; zeros(N - i, 1)]); end
            %dF1_N = @(X) dF1(X, N);
            %ddF1_num = my_gradient(dF1_N, X);
            %for i = 1:size(ddF1_num, 1)
            %    disp(['f(', num2str(i), ')'])
            %    disp(squeeze(ddF1_num(i, :, :))); disp(ddF1{i}(X, N));
            %end
            
            F2 = @(X, N) X(2:N) - X(1:N-1);
            dF2 = @(X, N) [zeros(N - 1, 1), eye(N - 1, N - 1)] - eye(N-1, N);
            for i = 1:N, ddF2{i} = @(X, N) zeros(N - 1, N); end
            %dF2_N = @(X) dF2(X, N);
            %ddF2_num = my_gradient(dF2_N, X);
            %for i = 1:size(ddF2_num, 1)
            %   disp(['f(', num2str(i), ')'])
            %   disp(squeeze(ddF2_num(i, :, :))); disp(ddF2{i}(X, N));
            %end
            F = @(X, N) [F1(X, N); sqrt(w2) * F2(X, N)];
            dF = @(X, N) [dF1(X, N); sqrt(w2) * dF2(X, N)];
            E = @(X, N) F(X, N)' * F(X, N);
            dE = @(X, N) 2 *  F(X, N)' * dF(X, N);
            
            ddE_an = 2 *  dF(X, N)' * dF(X, N);
            for i = 1:N
                ddE_an(i, :) = ddE_an(i, :) + 2 * F(X, N)' * [ddF1{i}(X, N); ddF2{i}(X, N)];
            end
            
            %E_N = @(X) E(X, N);
            %dE_N = @(X) dE(X, N);
            %dE_dX_num = my_gradient(E_N, X);
            %disp([dE(X, N); dE_dX_num]);
            %ddE_num = my_gradient(dE_N, X);
            
            %ddE_fun = compute_hessian(X, Y, T, N, w2);
            %for i = 1:size(ddE_num, 1)
            %    disp(['ddf(', num2str(i), ')'])
            %    disp(squeeze(ddE_num(i, :, :))); disp(ddE_fun(i, :));
            %end
            
            x = X(N);
            X_ = X(1:N - 1);
            
            e2_sqrt = @(x) exp(T(N) * x).^2 - Y(N);
            de2_sqrt_dx = @(x) 2 * T(N) * exp(T(N) * x).^2;
            
            e12_sqrt = @(X_, x) sqrt(w2) * (x - [zeros(1, N - 2), 1] * X_);
            de12_sqrt_dX_ = @(X_, x) - sqrt(w2) * [zeros(1, N - 2), 1];
            
            E1_sqrt = @(X_) [F1(X_, N - 1); sqrt(w2) * F2(X_, N - 1)];
            dE1_sqrt_dX_ = @(X_) [dF1(X_, N - 1); sqrt(w2) * dF2(X_, N - 1)];
            %dE1_sqrt_dX_num = my_gradient(E1_sqrt, X_);
            %disp([dE1_sqrt_dX_(X_), dE1_sqrt_dX_num]);
            
            e2 = @(x) e2_sqrt(x)' * e2_sqrt(x);
            de2_dX_ = @(X_) zeros(1, N - 1);
            de2_dx = @(x)  2 * e2_sqrt(x)' * de2_sqrt_dx(x);
            
            e12 = @(X_, x) e12_sqrt(X_, x)' * e12_sqrt(X_, x);
            de12_dX_ = @(X_, x)  2 * e12_sqrt(X_, x)' * de12_sqrt_dX_(X_, x);
            %e12_X = @(X_) e12(X_, x);
            %de12_dX_num = my_gradient(e12_X, X_);
            %disp([de12_dX_(X_, x), de12_dX_num]);
            
            E1 = @(X_) E1_sqrt(X_)' * E1_sqrt(X_);
            dE1_dX_ = @(X_)  2 * E1_sqrt(X_)' * dE1_sqrt_dX_(X_);
            %dE1_dX_num = my_gradient(E1, X_);
            %disp([dE1_dX_(X_), dE1_dX_num]);
            
            E2 = @(X_, x) e2(x) + e12(X_, x) + E1(X_);
            dE2_dX_ = @(X_, x) de2_dX_(x) + de12_dX_(X_, x) + dE1_dX_(X_);
            %disp([E(X, N) E2(X_, x)]);
            %E2_X_ = @(X_) E2(X_, x);
            %dE2_dX_num = my_gradient(E2_X_, X_);
            %disp([dE2_dX_(X_, x), dE2_dX_num]);
            
            %% Optimize X_
            
            X_init = X_;
            A = E1(X_init);
            B = dE1_dX_(X_init);
            C = compute_hessian(X_init, Y, T, N - 1, w2);
            
            %disp(C); disp(B' * B);
            
            K = [zeros(1, N - 2), 1];
            X_opt = X_init;
            E2_lin = @(X_opt) e2(x) + e12(X_opt, x) + A + B * (X_opt - X_init) + (X_opt - X_init)' * C * (X_opt - X_init);
            o2 = @(X_opt) e2(x) + w2 * (x - K * X_opt)' * (x - K * X_opt) + A + B * (X_opt - X_init) + (X_opt' * C * X_opt - 2 * X_init' * C * X_opt + X_init' * C * X_init);
            o3 = @(X_opt) e2(x) + w2 * (x' * x  - 2 * x' * K * X_opt + X_opt' * K' * K * X_opt) + A + B * (X_opt - X_init) + (X_opt' * C * X_opt - 2 * X_init' * C * X_opt + X_init' * C * X_init);
            %disp([E2_lin(X_opt); o3(X_opt)]);
            
            do3_dX_opt = @(X_opt)  - 2 * w2 * x' * K + 2 * w2 * X_opt' * K' * K + B +  X_opt' * (C + C') - 2 * X_init' * C;
            do3_dX_opt2 = @(X_opt)  (- 2 * w2 * x' * K + B - 2 * X_init' * C) + X_opt' * (2 * w2 * K' * K  +  C + C');
            do3_num = my_gradient(o3, X_opt);
            %disp([do3_num; do3_dX_opt2(X_opt)]);
            
            X_opt = - (2 * w2 * K' * K  +  C + C') \ (- 2 * w2 * x' * K + B - 2 * X_init' * C)';
            %disp(do3_dX_opt(X_opt));
            
            %% Optimize x
            X_opt = @ (x) - (2 * w2 * K' * K  +  C + C') \ (- 2 * w2 * x' * K + B - 2 * X_init' * C)';
            dX_opt = @(x) - (2 * w2 * K' * K  +  C + C') \ (- 2 * w2 * K)';
            
            e12_sqrt_opt = @(x) sqrt(w2) * (x - [zeros(1, N - 2), 1] * X_opt(x));
            de12_sqrt_opt = @(x) sqrt(w2) * (1 - [zeros(1, N - 2), 1] * dX_opt(x));
            e12_opt = @(x) e12_sqrt_opt(x)' * e12_sqrt_opt(x);
            de12_opt_dx = @(x) 2 * e12_sqrt_opt(x)' * de12_sqrt_opt(x);
            
            %de12_opt_num = my_gradient(e12_opt, x);
            %disp([de12_opt_num, de12_opt(x)]);
            
            E2_opt0 = @(x) e2(x) + e12_opt(x) + A + B * (X_opt(x) - X_init) + (X_opt(x) - X_init)' * C * (X_opt(x) - X_init);
            E2_opt = @(x) e2(x) + e12_opt(x) + B * X_opt(x) + X_opt(x)' * C * X_opt(x) - 2 * X_init' * C * X_opt(x) + A - B * X_init + X_init' * C * X_init;
            dE2_opt = @(x) de2_dx(x) + de12_opt_dx(x) + B * dX_opt(x) + 2 * X_opt(x)' * C * dX_opt(x) - 2 * X_init' * C * dX_opt(x);
            
            M1 = - (2 * w2 * K' * K  +  C + C') \ (- 2 * w2 * K)';
            M2 = - (2 * w2 * K' * K  +  C + C') \ (B - 2 * X_init' * C)';
            
            E2_opt1 = @(x) norm(exp(T(N) * x).^2 - Y(N))^2 + ...
                w2 * norm(x - [zeros(1, N - 2), 1] * X_opt(x))^2 + ...
                B *(M1 * x + M2) + (M1 * x + M2)' * C * (M1 * x + M2) - 2 * X_init' * C * (M1 * x + M2) + A - B * X_init + X_init' * C * X_init;
            
            E2_opt2 = @(x) norm(exp(T(N) * x).^2 - Y(N))^2 + ...
                w2 * norm(x - [zeros(1, N - 2), 1] * X_opt(x))^2 + ...
                B * M1 * x + B * M2 + x' * M1' * C * M1 * x + M2' * C  * M1 * x + x' * M1' * C * M2 + M2' * C * M2 - 2 * X_init' * C * M1 * x  - 2 * X_init' * C *  M2 + A - B * X_init + X_init' * C * X_init;
            
            E2_opt3 = @(x) norm(exp(T(N) * x).^2 - Y(N))^2 + ...
                w2 * norm(x - [zeros(1, N - 2), 1] * X_opt(x))^2 + ...
                x' * M1' * C * M1 * x + ...
                (B * M1 + 2 * M1' * C * M2 - 2 * X_init' * C * M1) * x  + ...
                - 2 * X_init' * C *  M2 + A - B * X_init + X_init' * C * X_init + B * M2 + M2' * C * M2;
            
            m3 = B * M1 + 2 * M1' * C * M2 - 2 * X_init' * C * M1;
            m4 = - 2 * X_init' * C *  M2 + A - B * X_init + X_init' * C * X_init + B * M2 + M2' * C * M2;
            
            [U, S, V] = eig(C);
            S1 = S.^0.5;
            if any(imag(S1))
                disp('imaginary');
            end
            M5 = S1' * U' * M1;
            m6 = (-0.5 * m3 / M5)';
            
            E2_opt4 = @(x) norm(exp(T(N) * x).^2 - Y(N))^2 + ...
                w2 * norm(x - [zeros(1, N - 2), 1] * (M1 * x + M2))^2 + ...
                x' * (M1' * U * E * U' * M1) * x + m3 * x  + m4;
            
            E2_opt5 = @(x) norm(exp(T(N) * x).^2 - Y(N))^2 + ...
                w2 * norm(x - [zeros(1, N - 2), 1] * (M1 * x + M2))^2 + ...
                x' * (M1' * U * S1 * S1' * U' * M1) * x + m3 * x  + m4;
            
            E2_opt6 = @(x) norm(exp(T(N) * x).^2 - Y(N))^2 + ...
                w2 * norm(x - [zeros(1, N - 2), 1] * (M1 * x + M2))^2 + ...
                x' * (M5' * M5) * x - 2 * m6' * M5 * x  + m6' * m6 - m6' * m6 + m4;
            
            E2_opt7 = @(x) norm(exp(T(N) * x).^2 - Y(N))^2 + ...
                w2 * norm(x - [zeros(1, N - 2), 1] * (M1 * x + M2))^2 + ...
                norm(M5 * x - m6)^2;
            
            F2_lin = @(x) [exp(T(N) * x).^2 - Y(N); sqrt(w2) * (x - [zeros(1, N - 2), 1] * (M1 * x + M2)); M5 * x - m6];
            
            %[X_opt(x)'; (M1 * x + M2)']
            %disp([E2_lin(X_init), E2_opt(x)]);
            %disp([E2_opt(x); E2_opt4(x)]);
            
            %% Optimize with fminunc
            %options = optimoptions(@fminunc,'Algorithm','trust-region', 'SpecifyObjectiveGradient',true);
            options = optimoptions(@fminunc,'Algorithm','quasi-newton');
            [x_opt, fval, exitflag, output] = fminunc(E2_opt4, x, options);
            
            return
            
            %% Plot E2
            if N == 2
                X_values = -2:0.01:-0.5;
                E2_values = zeros(length(X_values), 1);
                for l = 1:length(X_values)
                    E2_values(l) = E2(X_values(l), x);
                end
                figure; hold on; plot(X_values, E2_values, 'lineWidth', 2);
                mypoint([X_, E2(X_, x)], 'g', 30);
                mypoint([X_opt, E2(X_opt, x)], 'r', 30);
                %mypoint([x, E2(x, x)], 'k', 30);
            end
            return;
            
        end
        
        %% Solve
        if iter == num_iter, break; end
        
        w2 = 10;
        delta = (J1' * J1 + w2 * J2' * J2 + 0.2 * eye(N, N)) \ (J1' * F1 + w2 * J2' * F2);
        X = X + delta;
        
    end
    history{N}.J1 = J1;
    history{N}.J2 = J2;
    history{N}.F1 = F1;
    history{N}.F2 = F2;
end
%%}
return
%% Display history
offset = 0.09;
figure('units', 'normalized', 'outerposition', [0.1, 0.3, 0.45, 0.55]); hold on;
set(gca,'position', [0.05 0.08 0.93 0.90], 'units','normalized');
plot(1:length(history) + 1, x_true * ones(length(history) + 1, 1), 'lineWidth', 2, 'color', [1, 0.5, 0.4], 'lineStyle', '-.');
for j = 1:N
    myline([j, -1.5], [j, -0.5], [0.85, 0.85, 0.85], 2);
    for k = 1:j
        importance = history{j}.JtJ(k, k);
        myline([j + offset * k, history{j}.X(k) + 0.2 * sqrt(importance)], ...
            [j + offset * k, history{j}.X(k) - 0.2 * sqrt(importance)], [1, 0.9, 0.3], 3.5);
        mypoint([j + offset * k, history{j}.X(k)], [0.3, 0.6, 0.8], 20);
    end
end
set(gca, 'fontSize', 13);



