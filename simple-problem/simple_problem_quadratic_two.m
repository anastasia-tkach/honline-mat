function [xx_opt, dE_opt] = simple_problem_quadratic_two(X_prev, Y, T, N, w2)

y_ = Y(N - 1);
t_ = T(N - 1);
y = Y(N);
t = T(N);
x_0 = X_prev(N - 1);
if (N > 2)
    x__0 = X_prev(N - 2);
end

%% Data term
f1 =  exp(x_0 * t_)^2 - y_;
j1 = 2 * t_ * exp(x_0 * t_)^2;
h1 = 4 * t_.^2 * exp(x_0 * t_).^2;

%% Closeness term
f2 = sqrt(w2) * (x_0 - x__0);
j2 = sqrt(w2) * 1;
h2 = 0;

%% Energies
E1_ = f1' * f1; E2_ = f2' * f2;
dE1_ = 2 * f1' * j1; dE2_ = 2 * f2' * j2;
ddE1_ = 2 * j1' * j1 + 2 * f1 * h1; 
ddE2_ = 2 * j2' * j2 + 2 * f2 * h2;
E_ = E1_ + E2_; dE_ = dE1_ + dE2_; ddE_ = ddE1_ + ddE2_;

%% Quadratic approximation
xx0 = X_prev(N-1:N);
options = optimoptions(@fminunc,'Algorithm','quasi-newton', 'Display', 'off', 'SpecifyObjectiveGradient', true);
[xx_opt] = fminunc(@(xx) simple_problem_fg_quadratic_two(xx, x_0, y, t, E_, dE_, ddE_, w2), xx0, options);

%% Return gradients
[E_opt, dE_opt, ~, ~, ~, ~, ~, ~] = simple_problem_fg_quadratic_two(xx_opt, x_0, y, t, E_, dE_, ddE_, w2);

%% Plot gradient wrt x1
return
for x_index = 1:2

x_values = -1.8:0.001:-0.2;
E_values = zeros(length(x_values), 1); dE_values = zeros(length(x_values), 1);
Q_values = zeros(length(x_values), 1); dQ_values = zeros(length(x_values), 1);
e1_values =  zeros(length(x_values), 1); de1_values =  zeros(length(x_values), 1);
e2_values =  zeros(length(x_values), 1); de2_values =  zeros(length(x_values), 1);
for i = 1:length(x_values)
    if x_index == 1
        xx = [x_values(i), xx_opt(2)];
    else
        xx = [xx_opt(1), x_values(i)];
    end
    [E, dE, Q, dQ, e1, de1, e2, de2] = simple_problem_fg_quadratic_two(xx, x_0, y, t, E_, dE_, ddE_, w2);
    E_values(i) = E;
    dE_values(i) = dE(x_index);    
    Q_values(i) = Q;
    dQ_values(i) = dQ(x_index);
    e1_values(i) = e1;
    de1_values(i) = de1(x_index); 
    e2_values(i) = e2;
    de2_values(i) = de2(x_index); 
end
figure; hold on;
plot(x_values, E_values, 'lineWidth', 2); plot(x_values, dE_values, 'lineWidth', 2, 'lineStyle', '-.');
plot(x_values, Q_values, 'lineWidth', 1); plot(x_values, dQ_values, 'lineWidth', 1, 'lineStyle', '-.');
plot(x_values, e1_values, 'lineWidth', 1); plot(x_values, de1_values, 'lineWidth', 1, 'lineStyle', '-.');
plot(x_values, e2_values, 'lineWidth', 1); plot(x_values, de2_values, 'lineWidth', 1, 'lineStyle', '-.');
mypoint([xx_opt(x_index); dE_opt(x_index)], 'r', 20);
mypoint([xx_opt(x_index); E_opt], 'b', 20);
%myline([-1.2; 0], [-0.2; 0], [0.75, 0.75, 0.75]);
legend({'E', 'dE', 'Q', 'dQ', 'e1', 'de1', 'e2', 'de2'});
end










