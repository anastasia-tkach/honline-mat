function [F, J] = jacobian_shape_pose_cpp_wrapper(segments, joints, model_points, data_points, segment_indices)

num_points = length(data_points);
num_segments = length(segments);
num_joints = length(joints);
max_kinematic_chain = 3;

DataPoints = zeros(num_points, 2);
for j = 1:length(data_points)
    DataPoints(j, :) = data_points{j}';
end
ModelPoints = zeros(num_points, 2);
for j = 1:length(model_points)
    ModelPoints(j, :) = model_points{j}';
end
SegmentsGlobal =  zeros(num_segments, 16);
for j = 1:length(segments)
    SegmentsGlobal(j, :) = segments{j}.global(:)';
end
SegmentsKinematicChain =  -1 * ones(num_segments, max_kinematic_chain);
for j = 1:length(segments) - 1
    SegmentsKinematicChain(j, 1:length(segments{j}.kinematic_chain)) = segments{j}.kinematic_chain';
end
JointsSegmentId = zeros(num_joints, 1);
for j = 1:num_joints
    JointsSegmentId(j) = joints{j}.segment_id;
end
JointsAxis = zeros(num_joints, 3);
for j = 1:num_joints
    JointsAxis(j, :) = joints{j}.axis';
end

[F, J] = jacobian_shape_pose_cpp([num_points, num_joints, num_segments, max_kinematic_chain], DataPoints, ModelPoints, segment_indices, SegmentsKinematicChain, SegmentsGlobal, JointsSegmentId, JointsAxis);

%[F, J] = jacobian_shape_pose_matlab([num_points, num_joints, num_segments, max_kinematic_chain], DataPoints, ModelPoints, segment_indices, SegmentsKinematicChain, SegmentsGlobal, JointsSegmentId, JointsAxis);



