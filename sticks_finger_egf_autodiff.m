function [F, J, H] = sticks_finger_egf_autodiff(x, data)

num_points = 15;
num_segments = 4;
num_joints = 4;

%% Parse data
current_index = 0;
DataPoints = reshape(data(current_index + 1:current_index + num_points * 2), num_points, 2); 
current_index = current_index + num_points * 2;
%{
SegmentsLocal = reshape(data(current_index + 1:current_index + num_segments * 16), num_segments, 16); 
current_index = current_index + num_segments * 16;

SegmentsKinematicChain = reshape(data(current_index + 1:current_index + num_segments * 3), num_segments, 3); 
current_index = current_index + num_segments * 3;

SegmentsParentId = data(current_index + 1:current_index + num_segments); 
current_index = current_index + num_segments;

JointsSegmentId = data(current_index + 1:current_index + num_joints); 
current_index = current_index + num_joints;

JointsAxis = reshape(data(current_index + 1:current_index + num_joints * 3), num_joints, 3); 
%}
%% Create normal data structures
% data
data_points = cell(num_points, 1);
for j = 1:num_points
    data_points{j} = DataPoints(j, :)';
end
%{
% segments
segments = cell(num_segments, 1);
for j = 1:num_segments
   segments{j}.local =  reshape(SegmentsLocal(j, :), 4, 4);
end
for j = 1:num_segments - 1
    segments{j}.kinematic_chain = SegmentsKinematicChain(j, :);
    chain_length = find(segments{j}.kinematic_chain == -1, 1) - 1;
    if isempty(chain_length)
        chain_length = 3;
    end
    segments{j}.kinematic_chain = segments{j}.kinematic_chain(1:chain_length);
end
for j = 1:num_segments
    segments{j}.parent_id = SegmentsParentId(j);
end

% joints
joints = cell(num_joints, 1);
for j = 1:num_joints
    joints{j}.segment_id = JointsSegmentId(j);
end
for j = 1:num_joints
    joints{j}.axis = JointsAxis(j, :)';
    joints{j}.type = 'R';
end
%}
[segments0, joints] = segments_and_joints_2D();
blocks = {[1, 2], [2, 3], [3, 4]};
B = 3; T = 3;

%% Update
betas = x(1:B);
thetas = x(B + 1:B + T);

%% Initialize
[segments] = shape_2D_autodiff(segments0, betas);
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

