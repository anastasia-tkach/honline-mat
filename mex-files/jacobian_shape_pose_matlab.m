function [F, J, H] = jacobian_shape_pose_matlab(beta, theta, sizes, DataPoints, ModelPoints, segment_indices, SegmentsKinematicChain, SegmentsGlobal, JointsSegmentId, JointsAxis)

num_points = sizes(1);
num_joints = sizes(2);
num_segments = sizes(3);
max_kinematic_chain = sizes(4);

F = zeros(num_points, 1);
J = zeros(num_points, num_joints - 1 + num_segments - 1);
H = zeros(num_points, num_joints - 1 + num_segments - 1, num_joints - 1 + num_segments - 1);

%% Build the Jacobian matrix
for k = 1:num_points
    d = [DataPoints(k, :)'; 0];
    m = [ModelPoints(k, :)'; 0];
    
    if isempty(m) || norm(m - d) == 0, continue; end
    n = (d - m) / norm(m - d);    
      
    segment_kinematic_chain = SegmentsKinematicChain(segment_indices(k), :);        
    
    %% compute hessian - function
    bt = [beta; theta];
    if segment_kinematic_chain(2) == -1
        k1 = norm(m - SegmentsGlobal(1, 13:15)') / bt(1);
        %k1 = 1;
        k2 = 0; k3 = 0;
    end
    if segment_kinematic_chain(3) == -1 && segment_kinematic_chain(2) > 0
        k1 = 1; k3 = 0;
        k2 = norm(m - SegmentsGlobal(2, 13:15)') / bt(2);
        %k2 = 1;
    end
    if segment_kinematic_chain(3) > 0
        k1 = 1; k2 = 1;
        k3 = norm(m - SegmentsGlobal(3, 13:15)') / bt(3);
        %k3 = 1;
    end
    m =  [k1 * bt(1) * cos(pi/2 + bt(4)) + k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) + k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6));...
        k1 * bt(1) * sin(pi/2 + bt(4)) + k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) + k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    %disp([m_, m(1:2)]);
    
    %% compute hessian - gradients
    
    dm_dt1 =  [- k1 * bt(1) * sin(pi/2 + bt(4)) - k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6)); ...
        k1 * bt(1) * cos(pi/2 + bt(4)) + k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) + k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6))];
    dm_dt2 =  [ - k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6)); ...
        k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) + k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6))];
    dm_dt3 =  [ - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6)); ...
        k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_db1 =  [k1 * cos(pi/2 + bt(4)); k1 * sin(pi/2 + bt(4))];
    dm_db2 =  [k2 * cos(pi/2 + bt(4) + bt(5)); k2 * sin(pi/2 + bt(4) + bt(5))];
    dm_db3 =  [k3 * cos(pi/2 + bt(4) + bt(5) + bt(6)); k3 * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    j =  [dm_db1, dm_db2, dm_db3, dm_dt1, dm_dt2, dm_dt3];
    
    %% compute hessian - beta-beta
    dm_db1_db1 =  [0; 0]; dm_db1_db2 =  [0; 0]; dm_db1_db3 =  [0; 0];
    dm_db2_db1 =  [0; 0]; dm_db2_db2 =  [0; 0]; dm_db2_db3 =  [0; 0];
    dm_db3_db1 =  [0; 0]; dm_db3_db2 =  [0; 0]; dm_db3_db3 =  [0; 0];
    
    %% compute hessian - beta-theta
    dm_db1_dt1 =  [- k1 * sin(pi/2 + bt(4)); k1 * cos(pi/2 + bt(4))];
    dm_db1_dt2 =  [0; 0];
    dm_db1_dt3 =  [0; 0];
    
    dm_db2_dt1 =  [- k2 * sin(pi/2 + bt(4) + bt(5)); k2 * cos(pi/2 + bt(4) + bt(5))];
    dm_db2_dt2 =  [- k2 * sin(pi/2 + bt(4) + bt(5)); k2 * cos(pi/2 + bt(4) + bt(5))];
    dm_db2_dt3 =  [0; 0];
    
    dm_db3_dt1 =  [- k3 * sin(pi/2 + bt(4) + bt(5) + bt(6)); k3 * cos(pi/2 + bt(4) + bt(5) + bt(6))];
    dm_db3_dt2 =  [- k3 * sin(pi/2 + bt(4) + bt(5) + bt(6)); k3 * cos(pi/2 + bt(4) + bt(5) + bt(6))];
    dm_db3_dt3 =  [- k3 * sin(pi/2 + bt(4) + bt(5) + bt(6)); k3 * cos(pi/2 + bt(4) + bt(5) + bt(6))];
    
    %% compute hessian - theta-theta
    
    % dt1_d..
    dm_dt1_dt1 =  [- k1 * bt(1) * cos(pi/2 + bt(4)) - k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k1 * bt(1) * sin(pi/2 + bt(4)) - k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_dt1_dt2 =  [- k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_dt1_dt3 =  [- k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    % dt3_d..
    dm_dt2_dt1 =  [ - k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_dt2_dt2 =  [ - k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_dt2_dt3 =  [ - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    % dt3_d..
    dm_dt3_dt1 =  [ - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_dt3_dt2 =  [ - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_dt3_dt3 =  [ - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_ddb1 =  [dm_db1_db1, dm_db1_db2, dm_db1_db3, dm_db1_dt1, dm_db1_dt2, dm_db1_dt3];
    dm_ddb2 =  [dm_db2_db1, dm_db2_db2, dm_db2_db3, dm_db2_dt1, dm_db2_dt2, dm_db2_dt3];
    dm_ddb3 =  [dm_db3_db1, dm_db3_db2, dm_db3_db3, dm_db3_dt1, dm_db3_dt2, dm_db3_dt3];
    dm_ddt1 =  [dm_db1_dt1, dm_db2_dt1, dm_db3_dt1, dm_dt1_dt1, dm_dt1_dt2, dm_dt1_dt3];
    dm_ddt2 =  [dm_db1_dt2, dm_db2_dt2, dm_db3_dt2, dm_dt2_dt1, dm_dt2_dt2, dm_dt2_dt3];
    dm_ddt3 =  [dm_db1_dt3, dm_db2_dt3, dm_db3_dt3, dm_dt3_dt1, dm_dt3_dt2, dm_dt3_dt3];    
    
    ddm =  [- n(1:2)' * dm_ddb1; - n(1:2)' * dm_ddb2; - n(1:2)' * dm_ddb3; - n(1:2)' * dm_ddt1; - n(1:2)' * dm_ddt2; - n(1:2)' * dm_ddt3];
    
    
    %% Sstore to matrices
    F(k) = n(1:2)' * (d(1:2) - m);
    J(k, :) = - n(1:2)' * j;
    H(k, :, :) = ddm;
end








