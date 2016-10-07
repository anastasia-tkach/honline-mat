function [F, J, H] = sticks_finger_eg_single(x, segments0, joints, data_points)

%disp(x');

blocks = {[1, 2], [2, 3], [3, 4]};
B = 3; T = 3;

%% Update
betas = x(1:B);
thetas = x(B + 1:B + T);

%% Initialize
[segments] = shape_2D(segments0, betas);
[segments] = pose_2D(segments, joints, thetas);

%% Compute correspondences
[segment_indices, model_points] = compute_correspondences_cpp_wrapper(segments, blocks, data_points);

%% Compute vector F, G, H
[f, j, h] = jacobian_shape_pose_cpp_wrapper(betas, thetas, segments, joints, model_points, data_points, segment_indices);

%% Compute scalar F, G, H
F = f' * f;    
J = 2 * f' * j;    
H = 2 * j' * j;
for i = 1:size(f, 1)
    H = H + 2 * f(i) * squeeze(h(i, :, :));
end

