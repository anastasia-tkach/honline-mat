function [F, J] = sticks_finger_fg_batch(X, x0, segments0, joints, frames, N, D, settings, w2)

B = 3; T = 3;
L = min(N, settings.batch_size);

%% Data term

F1 = zeros(D * L, 1);
J1 = zeros(D * L, (B + T) * L);
for i = 1:L
    [f1, j1] = sticks_finger_fg_data(X((B + T) * (i - 1) + 1:(B + T) * i), segments0, joints, frames{N - L + i});
    F1(D * (i - 1) + 1:D * i) = f1;
    J1(D * (i - 1) + 1:D * i, (B + T) * (i - 1) + 1:(B + T) * i) = j1;
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

if settings.batch_online && N > settings.batch_size
    beta_1 = X(1:B);
    beta_0 = x0(1:B);
    F2(B * (L - 1) + 1: B * L) = beta_0 - beta_1;
    J2(B * (L - 1) + 1: B * L, 1:B) = - eye(B, B);
end

%% Robust online batch
if settings.batch_online_robust && N > settings.batch_size
    beta_1 = X(1:B);
    beta_0 = x0(1:B);
    r = beta_0 - beta_1;
    dr = - eye(B, B);
    [f, df] = german_mcclure_kernel(r, dr);
    F2(B * (L - 1) + 1: B * L) = f;
    J2(B * (L - 1) + 1: B * L, 1:B) = df;
end


%% Assemble
F = [F1; sqrt(w2) * F2];
J = [J1; sqrt(w2) * J2];

