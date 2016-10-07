function [F, J] = jacobian_ik_2D(beta, theta, segments, joints, model_points, data_points, segment_indices)
num_model_points = size(model_points, 1);
J = zeros(num_model_points, length(joints) - 1);
F = zeros(num_model_points, 1);
M_ = zeros(num_model_points, 2);

%% Build the Jacobian matrix
for k = 1:num_model_points
    d = [data_points{k}; 0];
    m = [model_points{k}; 0];
    
    if isempty(m) || norm(m - d) == 0, continue; end
    n = (d - m) / norm(m - d);
    
    j = zeros(3, length(joints) - 1);
    
    segment = segments{segment_indices(k)};
    for l = 1:length(segment.kinematic_chain)
        joint_id = segment.kinematic_chain(l);
        segment_id = joints{joint_id}.segment_id;
        u = joints{joint_id}.axis;
   
        p = segments{segment_id}.global(1:3, 4);
        T = segments{segment_id}.global;
        v = T * [u; 1]; v = v(1:3) / v(4);
        v = v - p;        
        
        switch joints{joint_id}.type
            case 'R'
                j(:, joint_id) = cross(v, m - p)';
            case 'T'
                j(:, joint_id) = v;
        end  
    end
    
    bt = [beta; theta];
    if length(segment.kinematic_chain) == 1        
        k1 = norm(m - segments{1}.global(1:3, 4)) / bt(1); 
        k2 = 0; k3 = 0;
    end    
    if length(segment.kinematic_chain) == 2
        k1 = 1; k3 = 0;
        k2 = norm(m - segments{2}.global(1:3, 4)) / bt(2);        
    end    
    if length(segment.kinematic_chain) == 3
        k1 = 1; k2 = 1;        
        k3 = norm(m - segments{3}.global(1:3, 4)) / bt(3);
    end
    m_ = @(bt) [k1 * bt(1) * cos(pi/2 + bt(4)) + k2 * bt(2) * cos(pi/2 + bt(4) + bt(5)) + k3 * bt(3) * cos(pi/2 + bt(4) + bt(5) + bt(6));...
        k1 * bt(1) * sin(pi/2 + bt(4)) + k2 * bt(2) * sin(pi/2 + bt(4) + bt(5)) +  k3 * bt(3) * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    disp([m_(bt), m(1:2)]);
       
    %% accumulate sides
    J(k, :) = n' * j;
    F(k) = n' * (d - m);
end
