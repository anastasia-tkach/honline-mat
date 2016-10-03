function [] = simple_problem_plot_gradient(X, Y, T, dE1, N)

figure; hold on;
for i = 1:N
    x_values = -1.8:0.01:-0.3;
    dE1_values = zeros(length(x_values), 1);
    for l = 1:length(x_values)
        dE1_values(l) =  2 * (exp(x_values(l) * T(i))^2 - Y(i))' * (2 * T(i) * exp(x_values(l) * T(i))^2);
    end
    plot(x_values, dE1_values, 'lineWidth', 2);
    mypoint([X(i), dE1(i)], 'r', 30);
    ylim([-0.1, 0.5]);
    xlim([-1.8, -0.3]);
end
