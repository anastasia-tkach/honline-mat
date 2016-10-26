function [F, J, H] = jacobian_shape_pose_analytical(beta, theta, sizes, DataPoints, ModelPoints, segment_indices, SegmentsKinematicChain, SegmentsGlobal, JointsSegmentId, JointsAxis)

num_points = sizes(1);
num_joints = sizes(2);
num_segments = sizes(3);
max_kinematic_chain = sizes(4);

F = zeros(num_points, 1);
J = zeros(num_points, num_joints - 1 + num_segments - 1);
H = zeros(num_points, num_joints - 1 + num_segments - 1, num_joints - 1 + num_segments - 1);

F_ = @(bt) [];
J_ = @(bt) [];
H_ = @(bt) [];

if num_points == 0, return; end

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
        
        c = 1;
        if l == 3 || segment_kinematic_chain(l + 1) == -1
           c = norm(m - p) / beta(l); 
        end

        j(:, segment_id) = c * v;
        
        %% pose
        v = T * [u; 1]; v = v(1:3) / v(4);
        v = v - p;        
        j(:, num_segments - 1 + joint_id) = cross(v, m - p)';
    end
    
    %%{
    %% compute hessian - function
    bt = [beta; theta];
    if segment_kinematic_chain(2) == -1    
        k1 = norm(m - SegmentsGlobal(1, 13:15)') / bt(1); 
        k2 = 0; k3 = 0;
    end    
    if segment_kinematic_chain(3) == -1 && segment_kinematic_chain(2) > 0
        k1 = 1; k3 = 0;
        k2 = norm(m - SegmentsGlobal(2, 13:15)') / bt(2);        
    end    
    if segment_kinematic_chain(3) > 0
        k1 = 1; k2 = 1;        
        k3 = norm(m - SegmentsGlobal(3, 13:15)') / bt(3);
    end
    m_ = @(bt) [k1 * bt(1) * cos(pi/2 + bt(4)) + k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) + k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6));...
        k1 * bt(1) * sin(pi/2 + bt(4)) + k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) +  k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    %disp([m_(bt), m(1:2)]);
    
    %% compute hessian - gradients
    
    dm_dt1 = @(bt) [- k1 * bt(1) * sin(pi/2 + bt(4)) - k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6)); ...
        k1 * bt(1) * cos(pi/2 + bt(4)) + k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) + k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6))];
    dm_dt2 = @ (bt)  [ - k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6)); ...
        k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) + k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6))];
    dm_dt3 = @ (bt)  [ - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6)); ...
        k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6))];    
   
    dm_db1 = @(bt) [k1 * cos(pi/2 + bt(4)); k1 * sin(pi/2 + bt(4))];
    dm_db2 = @(bt) [k2 * cos(pi/2 + bt(4) + bt(5)); k2 * sin(pi/2 + bt(4) + bt(5))];
    dm_db3 = @(bt) [k3 * cos(pi/2 + bt(4) + bt(5) + bt(6)); k3 * sin(pi/2 + bt(4) + bt(5) + bt(6))];   
    
    dm_ = @(bt) [dm_db1(bt), dm_db2(bt), dm_db3(bt), dm_dt1(bt), dm_dt2(bt), dm_dt3(bt)];
   
    %% compute hessian - beta-beta
    dm_db1_db1 = @(bt) [0; 0]; dm_db1_db2 = @(bt) [0; 0]; dm_db1_db3 = @(bt) [0; 0];
    dm_db2_db1 = @(bt) [0; 0]; dm_db2_db2 = @(bt) [0; 0]; dm_db2_db3 = @(bt) [0; 0];
    dm_db3_db1 = @(bt) [0; 0]; dm_db3_db2 = @(bt) [0; 0]; dm_db3_db3 = @(bt) [0; 0];
    
    %% compute hessian - beta-theta    
    dm_db1_dt1 = @(bt) [- k1 * sin(pi/2 + bt(4)); k1 * cos(pi/2 + bt(4))];
    dm_db1_dt2 = @(bt) [0; 0];
    dm_db1_dt3 = @(bt) [0; 0];
    
    dm_db2_dt1 = @(bt) [- k2 * sin(pi/2 + bt(4) + bt(5)); k2 * cos(pi/2 + bt(4) + bt(5))];
    dm_db2_dt2 = @(bt) [- k2 * sin(pi/2 + bt(4) + bt(5)); k2 * cos(pi/2 + bt(4) + bt(5))];
    dm_db2_dt3 = @(bt) [0; 0];    
    
    dm_db3_dt1 = @(bt) [- k3 * sin(pi/2 + bt(4) + bt(5) + bt(6)); k3 * cos(pi/2 + bt(4) + bt(5) + bt(6))];
    dm_db3_dt2 = @(bt) [- k3 * sin(pi/2 + bt(4) + bt(5) + bt(6)); k3 * cos(pi/2 + bt(4) + bt(5) + bt(6))];
    dm_db3_dt3 = @(bt) [- k3 * sin(pi/2 + bt(4) + bt(5) + bt(6)); k3 * cos(pi/2 + bt(4) + bt(5) + bt(6))];    
    
    %% compute hessian - theta-theta
    
    % dt1_d..
    dm_dt1_dt1 = @(bt) [- k1 * bt(1) * cos(pi/2 + bt(4)) - k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k1 * bt(1) * sin(pi/2 + bt(4)) - k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_dt1_dt2 = @(bt) [- k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_dt1_dt3 = @(bt) [- k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
     % dt3_d..    
    dm_dt2_dt1 = @ (bt)  [ - k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_dt2_dt2 = @ (bt)  [ - k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
         - k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_dt2_dt3 = @ (bt)  [ - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
        - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];

    % dt3_d..
    dm_dt3_dt1 = @ (bt)  [ - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
    - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_dt3_dt2 = @ (bt)  [ - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
    - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    dm_dt3_dt3 = @ (bt)  [ - k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6)); ...
    - k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];    

    dm_ddb1 = @ (bt) [dm_db1_db1(bt), dm_db1_db2(bt), dm_db1_db3(bt), dm_db1_dt1(bt), dm_db1_dt2(bt), dm_db1_dt3(bt)];
    dm_ddb2 = @ (bt) [dm_db2_db1(bt), dm_db2_db2(bt), dm_db2_db3(bt), dm_db2_dt1(bt), dm_db2_dt2(bt), dm_db2_dt3(bt)];
    dm_ddb3 = @ (bt) [dm_db3_db1(bt), dm_db3_db2(bt), dm_db3_db3(bt), dm_db3_dt1(bt), dm_db3_dt2(bt), dm_db3_dt3(bt)];
    dm_ddt1 = @ (bt) [dm_db1_dt1(bt), dm_db2_dt1(bt), dm_db3_dt1(bt), dm_dt1_dt1(bt), dm_dt1_dt2(bt), dm_dt1_dt3(bt)];
    dm_ddt2 = @ (bt) [dm_db1_dt2(bt), dm_db2_dt2(bt), dm_db3_dt2(bt), dm_dt2_dt1(bt), dm_dt2_dt2(bt), dm_dt2_dt3(bt)];
    dm_ddt3 = @ (bt) [dm_db1_dt3(bt), dm_db2_dt3(bt), dm_db3_dt3(bt), dm_dt3_dt1(bt), dm_dt3_dt2(bt), dm_dt3_dt3(bt)];   
    
    %nn = my_gradient(dm_dt3, bt);
    %aa = dm_ddt3(bt);
    %disp([nn(1, :); aa(1, :); nn(2, :); aa(2, :)]);
    
    ddm = @(bt) shiftdim([- n(1:2)' * dm_ddb1(bt); - n(1:2)' * dm_ddb2(bt); - n(1:2)' * dm_ddb3(bt); - n(1:2)' * dm_ddt1(bt); - n(1:2)' * dm_ddt2(bt); - n(1:2)' * dm_ddt3(bt)], -1);
    H(k, :, :) = ddm(bt); 
    %%}
    %% Sstore to matrices 
    %j = j./10;
    %j(:, 4:6) = j(:, 4:6)./50;
    
    F(k) = n' * (d - m);   
    J(k, :) = - n' * j;
    
    %for x = 1:3
    %    myline(m, m - j(:, x), 'r');
    %end
    %for x = 4:6
    %    myline(m, m - j(:, x), 'k');
    %end
    
    
    %%{
    F_ = @(bt) [F_(bt); n(1:2)' * (d(1:2) - m_(bt))];
    J_ = @(bt) [J_(bt); - n(1:2)' * dm_(bt)];
    H_ = @(bt) [H_(bt); ddm(bt)];
    %%}
end 

%{
V = my_gradient(F_, bt);
figure; imagesc(J_(bt) - V); axis equal; colorbar;
VV = my_gradient(J_, bt);
H__ = H_(bt);
for i = 1:6
    figure; imagesc(H__(:, :, i) - VV(:, :, i)); axis equal; colorbar;
end
%}

%% Compute scalar functions
%{
f = @(bt) F_(bt)' * F_(bt);    
j = @(bt) 2 * F_(bt)' * J_(bt);    
v = my_gradient(f, bt);
%disp([v; j(bt)]);

h = @(bt) hessian_for_scalar_objective(F_(bt), J_(bt), H_(bt));
vv = my_gradient(j, bt);
%disp([vv; h(bt)]);
%}





