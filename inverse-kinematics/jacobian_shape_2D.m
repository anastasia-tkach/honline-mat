function [F, J] = jacobian_shape_2D(beta, theta, segments, model_points, data_points, segment_indices)
num_model_points = size(model_points, 1);
J = zeros(num_model_points, length(segments) - 1);
F = zeros(num_model_points, 1);
u = [0; 1; 0];

%% Build the Jacobian matrix
for k = 1:num_model_points
    d = [data_points{k}; 0];
    m = [model_points{k}; 0];
    
    if isempty(m) || norm(m - d) == 0, continue; end
    n = (d - m) / norm(m - d);
    
    j = zeros(3, length(segments) - 1);
    
    segment = segments{segment_indices(k)};
    for l = 1:length(segment.shape_chain)
        segment_id = segment.shape_chain(l);
        p = segments{segment_id}.global(1:3, 4);
        T = segments{segment_id}.global;
        v = T * [u; 1]; v = v(1:3) / v(4);
        v = v - p;
        
        c = 1;
        if l == length(segment.shape_chain), c = norm(m - p) / beta(l); end
        
        j(:, segment_id) = c * v;
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
    
    %disp([m_(bt)'; m(1:2)']);
    
    dm_db1 = @(bt) [k1 * cos(pi/2 + bt(4)); k1 * sin(pi/2 + bt(4))];
    dm_db2 = @(bt) [k2 * cos(pi/2 + bt(4) + bt(5)); k2 * sin(pi/2 + bt(4) + bt(5))];
    dm_db3 = @(bt) [k3 * cos(pi/2 + bt(4) + bt(5) + bt(6)); k3 * sin(pi/2 + bt(4) + bt(5) + bt(6))];
    
    %v = my_gradient(m_, bt);
    %disp([j(1:2, :), v(:, 1:3)]);
    %disp([j(1:2, :), dm_db1(bt), dm_db2(bt) dm_db3(bt)]);
    
    dm_db1_db1 = @(bt) [0; 0]; dm_db1_db2 = @(bt) [0; 0]; dm_db1_db3 = @(bt) [0; 0];
    dm_db2_db1 = @(bt) [0; 0]; dm_db2_db2 = @(bt) [0; 0]; dm_db2_db3 = @(bt) [0; 0];
    dm_db3_db1 = @(bt) [0; 0]; dm_db3_db2 = @(bt) [0; 0]; dm_db3_db3 = @(bt) [0; 0];
    
    
    %% accumulate sides
    J(k, :) = n' * j;
    F(k) = n' * (d - m);
end

















