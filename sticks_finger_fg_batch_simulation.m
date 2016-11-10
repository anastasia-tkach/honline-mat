function [F, J] = sticks_finger_fg_batch_simulation(X, x0, x_, segments0, joints, frames, N, settings, history)

B = 3; T = 3;
L = min(N, settings.batch_size);

%% Data term
F1 = zeros(B * L, 1);
J1 = zeros(B * L, (B + T) * L);
for i = 1:L    
    
    if (N <= settings.batch_size), frame_i = i; 
    else frame_i = N - settings.batch_size + i; end
    
    beta_i = X((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B);
    H = squeeze(history.hessian_independent(frame_i, :, :));

    mu = squeeze(history.mu_independent(frame_i, :))';
    %if (i < L)
    %  mu = x_(1:B);
    %else
    %  mu = squeeze(history.mu_independent(frame_i, :))';
    %end
    F1(B * (i - 1) + 1: B * i, 1) = sqrtm(H) * (beta_i - mu);
    J1(B * (i - 1) + 1: B * i, (B + T) * (i - 1) + 1:(B + T) * i) = [sqrtm(H), zeros(B, T)];
end

%% Closeness term
F2 = zeros(B * (L - 1), 1);
J2 = zeros(B * (L - 1), L * (B + T));
for i = 1:L - 1
    j = i + 1;
    
    beta_i = X((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B);
    beta_j = X((B + T) * (j - 1) + 1:(B + T) * (j - 1) + B);
    
    F2(B * (i - 1) + 1: B * i) = beta_i - beta_j;
    J2(B * (i - 1) + 1: B * i, (B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = eye(B, B);
    J2(B * (i - 1) + 1: B * i, (B + T) * (j - 1) + 1:(B + T) * (j - 1) + B) = - eye(B, B);
    
end

%% Online batch

if settings.batch_online && N > settings.batch_size && ~settings.batch_simulation_kalman
    beta_1 = X(1:B);
    beta_0 = x0(1:B);
    F2(B * (L - 1) + 1: B * L) = (beta_0 - beta_1);
    J2(B * (L - 1) + 1: B * L, 1:B) = - eye(B, B);
end

%% Online batch kalman
if settings.batch_simulation_kalman && N > settings.batch_size
    
    K = zeros(B, B);
    for i = 1:N - settings.batch_size
        K = K + squeeze(history.hessian_independent(i, :, :));
    end
    
    beta_1 = X(1:B);
    beta_0 = x0(1:B);
    F2(B * (L - 1) + 1: B * L) = sqrtm(K) * (beta_0 - beta_1);
    J2(B * (L - 1) + 1: B * L, 1:B) = - sqrtm(K);
end


%% Assemble
F = [F1; sqrt(settings.w2) * F2];
J = [J1; sqrt(settings.w2) * J2];

if settings.uniform_shape_prior || settings.constant_sum_shape_prior
    F = [F1; sqrt(settings.w2) * F2; sqrt(settings.w4) * F4];
    J = [J1; sqrt(settings.w2) * J2; sqrt(settings.w4) * J4];
end

