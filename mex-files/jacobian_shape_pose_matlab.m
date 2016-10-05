function [F, J] = jacobian_shape_pose_matlab(sizes, DataPoints, ModelPoints, segment_indices, SegmentsKinematicChain, SegmentsGlobal, JointsSegmentId, JointsAxis)

num_points = sizes(1);
num_joints = sizes(2);
num_segments = sizes(3);
max_kinematic_chain = sizes(4);

J = zeros(num_points, num_joints - 1 + num_segments - 1);
F = zeros(num_points, 1);

%% Build the Jacobian matrix
for k = 1:num_points
    d = [DataPoints(k, :)'; 0];
    m = [ModelPoints(k, :)'; 0];
    
    if isempty(m) || norm(m - d) == 0, continue; end
    n = (d - m) / norm(m - d);
    
    j = zeros(3, num_joints - 1 + num_segments - 1);
    
    segment_kinematic_chain = SegmentsKinematicChain(segment_indices(k), :);
    for l = 1:length(segment_kinematic_chain)
        if segment_kinematic_chain(l) == -1, break; end
        
        joint_id = segment_kinematic_chain(l);
        segment_id = JointsSegmentId(joint_id);
        u = JointsAxis(joint_id, :)';
        
        segment_global = reshape(SegmentsGlobal(segment_id, :), 4, 4);
        
        p = segment_global(1:3, 4);
        T = segment_global;
                        
        %% shape
        v = T * [0; 1; 0; 1]; v = v(1:3) / v(4);
        v = v - p;        
        j(:, segment_id) = v;
        
        %% pose
        v = T * [u; 1]; v = v(1:3) / v(4);
        v = v - p;        
        j(:, num_segments - 1 + joint_id) = cross(v, m - p)';
    end
    
    %% accumulate sides
    J(k, :) = n' * j;
    F(k) = n' * (d - m);
end