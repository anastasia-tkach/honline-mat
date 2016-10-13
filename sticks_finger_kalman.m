function [x, C] = sticks_finger_kalman(x, C, segments0, joints, data_points, num_iters)

DOES NOT WORK

B = 3; T = 3;
r = 0.2;
R = diag(r * ones(B + T, 1));
x0 = x;

I = diag([0.5 * ones(B, 1); 50 * ones(T, 1)]);
for iter = 1:num_iters
    dx = x - x0;
    [F1, J1, ~] = sticks_finger_fg_data(x, segments0, joints, data_points);

    J2 = sqrt(I);
    J = [J1; J2];
    F2 = zeros(B + T, 1);
    F = [F1; F2];
    
    %LHS = J' * J;
    %delta = - LHS \ (J' * F);
    
    LHS = J' * J + R * C;
    delta = -dx + LHS \ (J' * (F + J * dx));
    
    x = x + delta;
end

C(1:B, 1:B) = C(1:B, 1:B) + 1/r * (J1(:, 1:B)' * J1(:, 1:B) + I(1:B, 1:B));



