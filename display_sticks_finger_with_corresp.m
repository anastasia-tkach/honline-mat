function [] = display_sticks_finger_with_corresp(beta, theta, data_points, blocks)

figure('units', 'normalized', 'outerposition', [0.25, 0.275, 0.45, 0.7]);
axis off; axis equal; hold on;
[segments0, joints] = segments_and_joints_2D();
[segments0] = shape_2D(segments0, beta);
[segments] = pose_2D(segments0, joints, theta);
[~, model_points] = compute_correspondences_2D(segments, blocks, data_points);
clf; hold on; axis off; axis equal; set(gcf,'color','w');
display_sticks_finger(segments, data_points, model_points);
drawnow; pause(0.05);
