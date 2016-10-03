function [x, C_] = simple_problem_kalman(x, y, t, C_, num_iters)

r = 0.1; R = r; Q = 0;
x0 = x;
for iter = 1:num_iters
    dx = x - x0;
    f = exp(x * t)^2 - y;
    j = - 2 * t * exp(x * t)^2;
    LHS = j' * j + R * C_;
    delta = -dx + LHS \ (j' * (f + j * dx));
    x = x + delta;
end
C = C_ + 1/r * (j' * j);
C_ = inv(inv(C) + Q);
