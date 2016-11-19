function [F, J, H] = real_finger_data(x, segments0, joints, radii, blocks, data_points)

B = 3; T = 3;

%% initialize

beta = x(1:B);
theta = x(B + 1:B + T);
[segments] = shape_3D(segments0, beta);
[segments] = pose_3D(segments, joints, theta);
 
[model_points, segment_indices] = compute_correspondences_3D(segments, radii, blocks, data_points);

[F, J] = jacobian_cpp_wrapper(beta, theta, segments, joints, model_points, data_points, segment_indices);

