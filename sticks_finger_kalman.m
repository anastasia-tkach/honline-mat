function [] = sticks_finger_kalman()

F1 = f1{N};
J1 = j1{N};
I = diag([0.1 * ones(B, 1); 50 * ones(T, 1)]);
R = diag(r * ones(B + T, 1));
dx = [betas{N}; thetas{N}] - [beta0; theta0];
J2 = sqrt(I);  H = [J1; J2];
F2 = zeros(B + T, 1); F = [F1; F2];
LHS = H' * H + R * C_;
%delta = LHS \ (J1' * F1);
delta = -dx + LHS \ (H' * (F + H * dx));
if (iter == num_iters)
    C = C_ + 1/r * (J1' * J1 + I);
    C_ = inv(inv(C) + Q);
    P = inv(C);
end
delta = [zeros((N - 1) * (B + T), 1); delta];
