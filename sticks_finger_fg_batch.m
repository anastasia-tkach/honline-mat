function [F, J] = sticks_finger_fg_batch(X, segments0, joints, frames, N, D, batch_size, w2)

B = 3; T = 3;

%% Data term

F1 = zeros(D * N, 1);
J1 = zeros(D * N, N * (B + T));
for i = max(1, N - batch_size + 1):N    
    [f1, j1] = sticks_finger_fg_data(X((B + T) * (i - 1) + 1:(B + T) * i), segments0, joints, frames{i});    
    F1(D * (i - 1) + 1:D * i) = f1;
    J1(D * (i - 1) + 1:D * i, (B + T) * (i - 1) + 1:(B + T) * i) = j1;
end

%% Closeness term

F2 = zeros(2, 1);
J2 = zeros(2, N * (B + T));
count = 1;
for i = max(1, N - batch_size):N - 1
    j = i + 1;
    
    beta_i = X((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B);
    beta_j = X((B + T) * (j - 1) + 1:(B + T) * (j - 1) + B);
    
    F2(B * (count - 1) + 1: B * count) = beta_i - beta_j;
    if i > N - batch_size
        J2(B * (count - 1) + 1: B * count, (B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = eye(B, B);
    end
    J2(B * (count - 1) + 1: B * count, (B + T) * (j - 1) + 1:(B + T) * (j - 1) + B) = - eye(B, B);
    count = count + 1;
end

%% Assemble
F = [F1; sqrt(w2) * F2];
J = [J1; sqrt(w2) * J2];

