function [F, J] = sticks_finger_fg_data(x, segments0, joints, data_points)

%disp(x);

blocks = {[1, 2], [2, 3], [3, 4]};
B = 3; T = 3;

%% Initialize
betas = x(1:B);
thetas = x(B + 1:B + T);
[segments] = shape_2D(segments0, betas);
[segments] = pose_2D(segments, joints, thetas);

%% Compute correspondences
[segment_indices, model_points] = compute_correspondences_cpp_wrapper(segments, blocks, data_points);

%% Compute Jacobians
[F, J] = jacobian_shape_pose_cpp_wrapper(betas, thetas, segments, joints, model_points, data_points, segment_indices, 'cpp');
