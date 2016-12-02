function [F, J, H] = real_finger_data(x, segments0, joints, radii, blocks, data_points, skeleton)

B = 3; T = 3;

%% initialize

beta = x(1:B);
theta = x(B + 1:B + T);
[segments] = shape_3D(segments0, beta);
[segments] = pose_3D(segments, joints, theta);

if skeleton
    [segment_indices, model_points] = compute_correspondences_3D_skeleton(segments, blocks, data_points);
    axis_projections = model_points;
else
    [model_points, axis_projections, segment_indices] = compute_correspondences_3D(segments, radii, blocks, data_points);
end

% figure; clf; hold on; axis off; axis equal; set(gcf,'color','w');
% display_finger_3D(segments, radii, blocks, data_points, model_points);
% drawnow; pause(0.05); 

[F, J] = jacobian_cpp_wrapper(beta, theta, radii, segments, joints, model_points, data_points, axis_projections, segment_indices);

