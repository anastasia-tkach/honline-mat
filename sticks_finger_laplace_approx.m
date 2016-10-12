function [xx_opt, J, h] = sticks_finger_laplace_approx(X0, X_prev, h_, segments0, joints, frames, N, w2, num_iters)

B = 3; T = 3;

x_0 = X_prev((B + T) * (N - 2) + 1:(B + T) * (N - 1));

if N == 3
    [~, ~, h_] = sticks_finger_fg_laplace_approx(X_prev((B + T) * (N - 3) + 1:(B + T) * (N - 1)), [], segments0, joints, frames{N - 1}, h_, w2);
end

xx0 = X0((B + T) * (N - 2) + 1:(B + T) * N);

[xx_opt] = my_lsqnonlin(@(xx) sticks_finger_fg_laplace_approx(xx, x_0, segments0, joints, frames{N}, h_, w2), xx0, num_iters);
[~, J, h] = sticks_finger_fg_laplace_approx(xx_opt, x_0, segments0, joints, frames{N}, h_, w2);


