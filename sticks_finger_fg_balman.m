function [F, J] = sticks_finger_fg_balman(X, x0, x_, segments0, joints, frames, N, settings, history)

B = 3; T = 3;
L = min(N, settings.batch_size);

%% Data term
if settings.balman_solve_all
    F1 = zeros(0, 1);
    J1 = zeros(0, (B + T) * L);
    count = 0;
    for i = 1:L
        [f1, j1] = sticks_finger_fg_data(X((B + T) * (i - 1) + 1:(B + T) * i), segments0, joints, frames{N - L + i}, settings, 'cpp');
        num_points = size(f1, 1);
        F1(count + 1:count + num_points, 1) = f1;
        J1(count + 1:count + num_points, (B + T) * (i - 1) + 1:(B + T) * i) = j1;
        count = count + num_points;
    end
end

if settings.balman_simulate || settings.balman_solve_last
    F1 = zeros(B * L, 1);
    J1 = zeros(B * L, (B + T) * L);
    for i = 1:L
        if (N <= settings.batch_size), frame_i = i;
        else frame_i = N - settings.batch_size + i; end
        
        beta_i = X((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B);
        
        H = squeeze(history.hessian_independent(frame_i, :, :));
        mu = squeeze(history.mu_independent(frame_i, :))';
        
        F1(B * (i - 1) + 1: B * i, 1) = sqrtm(H) * (beta_i - mu);
        J1(B * (i - 1) + 1: B * i, (B + T) * (i - 1) + 1:(B + T) * i) = [sqrtm(H), zeros(B, T)];
    end
end

if settings.balman_solve_last
    [f1, j1] = sticks_finger_fg_data(X((B + T) * (L - 1) + 1:(B + T) * L), segments0, joints, frames{N}, settings, 'cpp');
    num_points = size(f1, 1);
    F1(B * (L - 1) + 1:B * (L - 1) + num_points, 1) = f1;
    J1(B * (L - 1) + 1:B * (L - 1) + num_points, (B + T) * (L - 1) + 1:(B + T) * L) = j1;
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

%% Prior term
if N > settings.batch_size
    
    if settings.balman_uniform_prior
        P = eye(B, B);
    end
    
    if settings.balman_kalman_prior
        P = zeros(B, B);
        for i = 1:N - settings.batch_size
            P = P + squeeze(history.hessian_independent(i, :, :));
        end
    end
    
    beta_1 = X(1:B);
    beta_0 = x0(1:B);
    F2(B * (L - 1) + 1: B * L) = sqrtm(P) * (beta_0 - beta_1);
    J2(B * (L - 1) + 1: B * L, 1:B) = - sqrtm(P);
end

%% Assemble
F = [F1; sqrt(settings.w2) * F2];
J = [J1; sqrt(settings.w2) * J2];


