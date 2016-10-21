function [xx_opt, J, h] = sticks_finger_laplace_approx(xx, x0, x_1, h_, segments0, joints, frames, N, w2, num_iters)

B = 3; T = 3;

if N == 3
    [~, ~, h_] = sticks_finger_fg_laplace_approx([x_1; x0], [], segments0, joints, frames{N - 1}, h_, w2);
end

%xx = X((B + T) * (N - 2) + 1:(B + T) * N);

[xx_opt] = my_lsqnonlin(@(xx) sticks_finger_fg_laplace_approx(xx, x0, segments0, joints, frames{N}, h_, w2), xx, num_iters);
[~, J, h] = sticks_finger_fg_laplace_approx(xx_opt, x0, segments0, joints, frames{N}, h_, w2);


