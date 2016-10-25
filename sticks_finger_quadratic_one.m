function [xx_opt, J, h] = sticks_finger_quadratic_one(X, x0, x_1, h_, segments0, joints, frames, N, w2, settings)

B = 3; T = 3;

if N == 3  
    [~, ~, h_] = sticks_finger_fg_quadratic_one(x0, x_1, segments0, joints, frames{N - 1}, [], w2, settings);   
end

x2 = X((B + T) + 1:2 * (B + T));
[x2_opt] = my_lsqnonlin(@(x2) sticks_finger_fg_quadratic_one(x2, x0, segments0, joints, frames{N}, h_, w2, settings), x2, settings.num_iters);

W2 = w2 * eye(B, B);

W3 = h_(B + T + 1:B + T + B, B + T + 1:B + T + B);

x1_opt = x0;
x1_opt(1:B) = inv(W2 + W3) * W2 * x2_opt(1:B) + inv(W2 + W3) * W3 * x0(1:B);

xx_opt = [x1_opt; x2_opt];

[~, J, h] = sticks_finger_fg_quadratic_one(x2_opt, x0, segments0, joints, frames{N}, h_, w2, settings);

