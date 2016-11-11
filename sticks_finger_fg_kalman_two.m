function [F, J] = sticks_finger_fg_kalman_two(xx, x0, segments0, joints, data_points, H, settings)

B = 3;  T = 3;
M = B + T;
num_data = length(data_points);

%% Data term
[F1, dF1_dx2, ~] = sticks_finger_fg_data(xx(M + 1:2 * M), segments0, joints, data_points, settings, 'cpp');
dF1 = zeros(num_data, 2 * M);
dF1(:, M + 1:2 * M) = dF1_dx2;

%% Closeness term
F2 = sqrt(settings.w2) * (xx(M + 1:M + B) - xx(1:B));
dF2_dx1 = [- sqrt(settings.w2) * eye(B, B), zeros(B, T)];
dF2_dx2 = [ sqrt(settings.w2) * eye(B, B), zeros(B, T)];
dF2 = [dF2_dx1, dF2_dx2];


F = [F1; F2];
J = [dF1; dF2];

%% Quadratic approx   
if ~isempty(x0)
    Q = sqrtm(H) * (xx(1:B) - x0(1:B));
    dQ_dx1 = [sqrtm(H), zeros(B, T)];
    dQ_dx2 = zeros(B, M);
    dQ = [dQ_dx1, dQ_dx2];

    F = [F1; F2; Q];
    J = [dF1; dF2; dQ];
end


