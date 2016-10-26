function [F, J] = sticks_finger_fg_kalman_like(x, x_prev, segments0, joints, data_points, JtJ, N, settings)

B = 3; T = 3;

%% Data term
[F1, J1, ~] = sticks_finger_fg_data(x, segments0, joints, data_points, settings);

%% Closeness term
F2 = zeros(B, 1);
J2 = zeros(B, B + T);
if (N > 1)
    F2 = sqrt(settings.w2) * sqrtm(JtJ) * (x(1:B) - x_prev(1:B));
    J2 = sqrt(settings.w2) * sqrtm(JtJ) * [eye(B, B), zeros(B, T)];
end

%% All terms
F = [F1; F2];
J = [J1; J2];


