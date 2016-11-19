function [F, J] = jacobian_cpp_wrapper(beta, theta, segments, joints, model_points, data_points, segment_indices)

num_segments = length(segments);
num_joints = length(joints);
max_kinematic_chain = 3;

DataPoints = zeros(0, 3);
ModelPoints = zeros(0, 3);
SegmentIndices = zeros(0, 1);
k = 1;
for j = 1:length(data_points)
    if isempty(model_points{j}), continue; end 
    if any(isinf(model_points{j})), continue; end  
    DataPoints(k, :) = data_points{j}';
    ModelPoints(k, :) = model_points{j}';
    SegmentIndices(k) = segment_indices(j);
    k = k + 1;
end
num_points = size(DataPoints, 1);

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
[F, J] = jacobian_analytical(beta, theta, [num_points, num_joints, num_segments, max_kinematic_chain], DataPoints, ModelPoints, SegmentIndices', SegmentsKinematicChain, SegmentsGlobal, JointsSegmentId, JointsAxis);

%[F, J] = jacobian_cpp(beta, theta, [num_points, num_joints, num_segments, max_kinematic_chain], DataPoints, ModelPoints, SegmentIndices', SegmentsKinematicChain, SegmentsGlobal, JointsSegmentId, JointsAxis);