function [F, J] = sticks_finger_fg_batch(X, segments0, joints, frames, N, D, batch_size, w2)

blocks = {[1, 2], [2, 3], [3, 4]};
B = 3; T = 3;

j1 = cell(N, 1);
f1 = cell(N, 1);

%% Update
betas = cell(N, 1);
thetas = cell(N, 1);
for i = 1:N
    betas{i} = X((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B);
    thetas{i} = X((B + T) * (i - 1) + B + 1:(B + T) * i);
end

for i = max(1, N - batch_size + 1):N
    data_points = frames{i};
    
    %% Initialize
    [segments] = shape_2D(segments0, betas{i});
    [segments] = pose_2D(segments, joints, thetas{i});
    
    %% Compute correspondences
    [segment_indices, model_points] = compute_correspondences_cpp_wrapper(segments, blocks, data_points);
    
    %% Compute Jacobians  
    beta = X((B + T) * (i - 1) + 1:(B + T) * (i - 1) + B);
    theta = X((B + T) * (i - 1) + B + 1:(B + T) * i);
    [f1_i, j1_i] = jacobian_shape_pose_cpp_wrapper(beta, theta, segments, joints, model_points, data_points, segment_indices, 'cpp');
    j1{i} =  j1_i;
    f1{i} = f1_i;
end

%% Build data jacobian
F1 = zeros(D * N, 1);
J1 = zeros(D * N, N * (B + T));
for i = max(1, N - batch_size + 1):N
    F1(D * (i - 1) + 1:D * i) = f1{i};
    J1(D * (i - 1) + 1:D * i, (B + T) * (i - 1) + 1:(B + T) * i) = j1{i};
end

%% Compute closeness
F2 = zeros(2, 1);
J2 = zeros(2, N * (B + T));
count = 1;
for i = max(1, N - batch_size):N - 1
    coefficients = [1; 1; 1];
    j = i + 1;
    F2(B * (count - 1) + 1: B * count) = diag(coefficients) * (betas{i} - betas{i + 1});
    if i > N - batch_size
        J2(B * (count - 1) + 1: B * count, (B + T) * (i - 1) + 1:(B + T) * (i - 1) + B) = diag(coefficients) * eye(B, B);
    end
    J2(B * (count - 1) + 1: B * count, (B + T) * (j - 1) + 1:(B + T) * (j - 1) + B) = - diag(coefficients) * eye(B, B);
    count = count + 1;
end

%% Assemble
F = [F1; sqrt(w2) * F2];
J = [J1; sqrt(w2) * J2];