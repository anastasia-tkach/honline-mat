function [F, J, h] = sticks_finger_fg_quadratic_one(x2, x0, segments0, joints, data_points, h_, settings)

B = 3;  T = 3;

%% Recursive term
if ~isempty(h_)
    W3 = h_(B + T + 1:B + T + B, B + T + 1:B + T + B);
else
    W3 = zeros(B, B);
end

W2 = settings.w2 * eye(B, B);

%% Optimal value of x1
M = inv(W2 + W3) * W2;
m = inv(W2 + W3) * W3 * x0(1:B);
x1_opt = M * x2(1:B) + m;

%% Energy terms
I = eye(B, B);

% data energy
[F1, dF1_dx2, ~] = sticks_finger_fg_data(x2, segments0, joints, data_points);

% closeness energy
F2 = sqrtm(W2) * (x2(1:B) - M * x2(1:B) - m);
dF2_dx2 = [sqrtm(W2) * (I - M),  zeros(B, T)];

% recursive energy
if ~isempty(h_)
    F3 = real(sqrtm(W3)) * (M * x2(1:B) + m - x0(1:B));
    dF3_dx2 = [real(sqrtm(W3)) * M, zeros(B, T)];
else
    F3 = zeros(B, 1);
    dF3_dx2 = zeros(B, B + T);
end

% shape prior
w4 = 1;
Q = ones(1, B);
W4 = sqrt(w4) * (eye(B, B) -  1/B * (Q' * Q));
F4 = W4 * x2(1:B);
dF4_dx2 = [W4,  zeros(B, T)];

%% Combining together
F = [F3; F1; F2];
J = [dF3_dx2; dF1_dx2; dF2_dx2];

if settings.shape_prior
    F = [F3; F1; F2; F4];
    J = [dF3_dx2; dF1_dx2; dF2_dx2; dF4_dx2];
end

h = 2 * J' * J;

h = [eye(B + T, B + T), eye(B + T, B + T); eye(B + T, B + T), h];






