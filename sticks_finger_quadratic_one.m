function [xx_opt, J, h] = sticks_finger_quadratic_one(X0, X_prev, h_, segments0, joints, frames, N, w2, settings)

B = 3; T = 3;

x_0 = X_prev((B + T) * (N - 2) + 1:(B + T) * (N - 1));

if N == 3  
    if settings.quadratic_one_marginalization
        [~, ~, h_] = sticks_finger_fg_laplace_approx(X_prev((B + T) * (N - 3) + 1:(B + T) * (N - 1)), [], segments0, joints, frames{N - 1}, h_, w2);
    else
        [~, ~, h_] = sticks_finger_fg_quadratic_one(X_prev((B + T) * (N - 2) + 1:(B + T) * (N - 1)), X_prev((B + T) * (N - 3) + 1:(B + T) * (N - 2)), segments0, joints, frames{N - 1}, [], w2, settings);
    end
    
end

x2 = X0((B + T) * (N - 1) + 1:(B + T) * N);
[x2_opt] = my_lsqnonlin(@(x2) sticks_finger_fg_quadratic_one(x2, x_0, segments0, joints, frames{N}, h_, w2, settings), x2, settings.num_iters);

W2 = w2 * eye(B, B);

a = h_(1:B, 1:B);
b = h_(1:B, B + T + 1:B + T + B);
c = h_(B + T + 1:B + T + B, 1:B);
d = h_(B + T + 1:B + T + B, B + T + 1:B + T + B);
if settings.quadratic_one_marginalization
    W3 = d - c * inv(a) * b;
else
    W3 = d;
end

x1_opt = x_0;
x1_opt(1:B) = inv(W2 + W3) * W2 * x2_opt(1:B) + inv(W2 + W3) * W3 * x_0(1:B);

xx_opt = [x1_opt; x2_opt];

if settings.quadratic_one_marginalization    
    [~, J, h] = sticks_finger_fg_laplace_approx(xx_opt, x_0, segments0, joints, frames{N}, h_, w2);
else
    [~, J, h] = sticks_finger_fg_quadratic_one(x2_opt, x_0, segments0, joints, frames{N}, h_, w2, settings);
end
