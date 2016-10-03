function [] = simple_problem_debug_recompute(x, X_, x_opt, X_opt, Y, T, E1, dE1, ddE1, N, w2)

K = [zeros(1, N - 2), 1];
M1 = - (2 * w2 * K' * K  +  ddE1 + ddE1') \ (- 2 * w2 * K)';
M2 = - (2 * w2 * K' * K  +  ddE1 + ddE1') \ (dE1 - 2 * X_' * ddE1)';

m3 = dE1 * M1 + 2 * M1' * ddE1 * M2 - 2 * X_' * ddE1 * M1;
m4 = - 2 * X_' * ddE1 *  M2 + E1 - dE1 * X_ + X_' * ddE1 * X_ + dE1 * M2 + M2' * ddE1 * M2;

[U, S, ~] = eig(ddE1);
S1 = S.^0.5;
if any(imag(S1))
    disp('imaginary');
end
M5 = S1' * U' * M1;
if (M5 == 0)
    m6 = 0;
else
    m6 = (-0.5 * m3 / M5)';
end

E2_lin =  @(x) norm(exp(T(N) * x).^2 - Y(N))^2 + w2 * norm(x - [zeros(1, N - 2), 1] * (M1 * x + M2))^2 + ...
    x' * (M1' * U * S1 * S1' * U' * M1) * x + m3 * x  + m4;

E2_real = @(x) norm(exp(T(N) * x).^2 - Y(N))^2 + w2 * norm(x - [zeros(1, N - 2), 1] * (M1 * x + M2))^2;

X_init = X_;
X_opt = @(x) - (2 * w2 * K' * K  +  ddE1 + ddE1') \ (- 2 * w2 * x' * K + dE1 - 2 * X_init' * ddE1)';
E2_expansion = @(x)  E1 + dE1 * (X_opt(x) - X_init) + (X_opt(x) - X_init)' * ddE1 * (X_opt(x) - X_init);

x_values = -2:0.01:-0.5;
E2_lin_values = zeros(length(x_values), 1);
E2_real_values = zeros(length(x_values), 1);
E2_expansion_values = zeros(length(x_values), 1);
for l = 1:length(x_values)
    E2_lin_values(l) = E2_lin(x_values(l));
    E2_real_values(l) = E2_real(x_values(l));
    E2_expansion_values(l) = E2_expansion(x_values(l));
end
figure; hold on;
plot(x_values, E2_lin_values, 'lineWidth', 2);
plot(x_values, E2_real_values, 'lineWidth', 2);
plot(x_values, E2_expansion_values, 'lineWidth', 2);
mypoint([x, E2_lin(x)], 'g', 30);
mypoint([x_opt, E2_lin(x_opt)], 'r', 30);

disp(' ');

%% Plot E(X_)
return

% Verify Xopt
E2_X = @(X_opt) norm(exp(T(N) * x).^2 - Y(N))^2 + w2 * norm(x - K * X_opt)^2 + ...
    history{N - 1}.E1 + history{N - 1}.dE1 * (X_opt - X_init) + (X_opt - X_init)' * history{N - 1}.ddE1 * (X_opt - X_init);
[X_opt_fminunc, ~, ~, ~] = fminunc(E2_X, X_init, options);

if 0%(N == 2)
    X_values = -2:0.01:-0.5;
    E2_X_values = zeros(length(X_values), 1);
    for l = 1:length(X_values)
        E2_X_values(l) = E2_X(X_values(l));
    end
    figure; hold on;
    plot(X_values, E2_X_values, 'lineWidth', 2, 'color', [0.55, 0.75, 0.8]);
    mypoint([X_opt_fminunc, E2_X(X_opt_fminunc)], 'r', 30);
    mypoint([X_init, E2_X(X_init)], [1, 0.6, 0.1], 30);
    %ylim([0, 0.15]);
end

if 0;%(N == 3)
    values = -2:0.05:-0.5;
    [X1_values, X2_values] = meshgrid(values, values);
    E2_X_values = zeros(length(values), 1);
    for l = 1:length(values)
        for m = 1:length(values)
            %disp([values(l), values(m)]);
            E2_X_values(l, m) = E2_X([values(l); values(m)]);
        end
    end
    figure; hold on;
    mesh(X1_values, X2_values, E2_X_values);
    mypoint([X_opt_fminunc; E2_X(X_opt_fminunc)], 'r', 30);
    mypoint([X_init; E2_X(X_init)], [1, 0.6, 0.1], 30);
    %mypoint([[-1.3623; -1.1778]; E2_X([-1.3623; -1.1778])], 'g', 30);
    view([90, 0]);
    disp(' ');
    %ylim([0, 0.15]);
end