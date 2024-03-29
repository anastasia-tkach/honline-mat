function [F, J, h] = sticks_finger_fg_laplace_approx(xx, x_0, segments0, joints, data_points, h_, settings)

B = 3;  T = 3;
M = B + T;
num_data = length(data_points);

%% Data term
[F1, dF1_dx2, ddF1_dx2_dx2] = sticks_finger_fg_data(xx(M + 1:2 * M), segments0, joints, data_points, settings, 'cpp');

% jacobian
dF1 = zeros(num_data, 2 * M);
dF1(:, M + 1:2 * M) = dF1_dx2;

% hesian
ddF1 = zeros(num_data, 2 * M, 2 * M);
ddF1(:, M + 1:2 * M, M + 1:2 * M) = ddF1_dx2_dx2;


%% Closeness term
F2 = sqrt(settings.w2) * (xx(M + 1:M + B) - xx(1:B));
dF2_dx1 = [- sqrt(settings.w2) * eye(B, B), zeros(B, T)];
dF2_dx2 = [ sqrt(settings.w2) * eye(B, B), zeros(B, T)];
dF2 = [dF2_dx1, dF2_dx2];

ddF2 = zeros(M, 2 * M, 2 * M);

%% Quadratic approx
if isempty(x_0)
    Q = zeros(B, 1);
    dQ = zeros(B, 2 * M);
    ddQ = zeros(B, 2 * M, 2 * M);
else
    a = h_(1:B, 1:B);
    b = h_(1:B, M + 1:M + B);
    c = h_(M + 1:M + B, 1:B);
    d = h_(M + 1:M + B, M + 1:M + B);
    
    if (settings.quadratic_two_marginalization)
        ddQ_dx2_dx2 = d - c * inv(a) * b;
    end
    if (settings.quadratic_two_maximization)
        ddQ_dx2_dx2 = d;
    end

    ddQ_dx2_dx2_sqrt = real(sqrtm(ddQ_dx2_dx2));
    
    Q = ddQ_dx2_dx2_sqrt * (xx(1:B) - x_0(1:B));
    dQ_dx1 = [ddQ_dx2_dx2_sqrt, zeros(B, T)];
    dQ_dx2 = zeros(B, M);
    dQ = [dQ_dx1, dQ_dx2];
    ddQ = zeros(B, 2 * M, 2 * M);
end

F = [Q; F1; F2];
J = [dQ; dF1; dF2];
H = [ddQ; ddF1; ddF2];

%% Hessian for scalar objective
df = F' * F;
j = 2 * F' * J;
%h = hessian_for_scalar_objective(F, J, H);
h = 2 * J' * J;

%% Expression for h
%{
h = 2 * (dQ(:, [1:B, M + 1:M + B])' * dQ(:, [1:B, M + 1:M + B]) + dF1(:, [1:B, M + 1:M + B])' * dF1(:, [1:B, M + 1:M + B]) + dF2(:, [1:B, M + 1:M + B])' * dF2(:, [1:B, M + 1:M + B]));
a = 2 * dQ(:, 1:B)' * dQ(:, 1:B)' + 2 * settings.w2 * eye(B, B);
d = 2 * dF1(:, M + 1:M + B)' *  dF1(:, M + 1:M + B) + 2 * settings.w2 * eye(B, B);
b = -2 * settings.w2 * eye(B, B);
c = -2 * settings.w2 * eye(B, B);
h = [a, b; c, d];
%}

%% Compare with quadratic-one
%{
if ~isempty(x_0)
    x2 = xx(B + T + 1:2 * (B + T));
    
    Fo = @(x1)  [ddQ_dx2_dx2_sqrt * (x1(1:B) - x_0(1:B)); F1; sqrt(w2) * (x2(1:B) - x1(1:B));];
    x1 = xx(1:B + T);
    [x1_opt] = lsqnonlin(Fo, x1);
    disp(Fo(x1_opt)' * Fo(x1_opt));
end
%}

