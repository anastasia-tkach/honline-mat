function [F, J] = sticks_finger_fg_batch(X, segments0, joints, frames, N, D, settings, w2)

B = 3; T = 3;

%% Data term

F1 = zeros(D * N, 1);
J1 = zeros(D * N, N * (B + T));
for i = max(1, N - settings.batch_size + 1):N    
    [f1, j1] = sticks_finger_fg_data(X((B + T) * (i - 1) + 1:(B + T) * i), segments0, joints, frames{i});    
    F1(D * (i - 1) + 1:D * i) = f1;
    J1(D * (i - 1) + 1:D * i, (B + T) * (i - 1) + 1:(B + T) * i) = j1;
end

%% Closeness term

F2 = zeros(2, 1);
J2 = zeros(2, N * (B + T));
for i = max(1, N - settings.batch_size):N - 1
    j = i + 1;
    
    beta_i = X((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B);
    beta_j = X((B + T) * (j - 1) + 1:(B + T) * (j - 1) + B);
        
    if i > N - settings.batch_size
        J2(B * (i - 1) + 1: B * i, (B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = eye(B, B);
    end
    if ~ settings.batch_independent || (settings.batch_independent && i > N - settings.batch_size)
        F2(B * (i - 1) + 1: B * i) = beta_i - beta_j;
        J2(B * (i - 1) + 1: B * i, (B + T) * (j - 1) + 1:(B + T) * (j - 1) + B) = - eye(B, B);
    end
end

%% Assemble
F = [F1; sqrt(w2) * F2];
J = [J1; sqrt(w2) * J2];

