function [F, J, H] = jacobian_analytical(beta, theta, radii, sizes, DataPoints, ModelPoints, AxisProjections, segment_indices, SegmentsKinematicChain, SegmentsGlobal, JointsSegmentId, JointsAxis)

[segments, joints] = segments_and_joints_3D();
[segments] = shape_3D(segments, beta);

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
    d = DataPoints(k, :)';
    m = ModelPoints(k, :)';
    s = AxisProjections(k, :)';
    %disp(d');
    %disp(m');
    
    if isempty(m) || norm(m - d) == 0, continue; end
    n = (d - m) / norm(m - d);
    
    j = zeros(3, num_joints - 1 + num_segments - 1);
    
    segment_kinematic_chain = SegmentsKinematicChain(segment_indices(k), :);
    kinematic_chain_globals = cell(length(segment_kinematic_chain), 1);
    for l = 1:length(segment_kinematic_chain)
        if segment_kinematic_chain(l) == -1, break; end
        joint_id = segment_kinematic_chain(l);
        segment_id = JointsSegmentId(joint_id);
        u = JointsAxis(joint_id, :)';
        
        %% extract data
        segment_global = reshape(SegmentsGlobal(segment_id, :), 4, 4);
        t = segment_global(1:3, 4);
        T = segment_global;
        
        %% pose
        v = T * [u; 1];
        v = v(1:3) / v(4);
        v = v - t;
        j(:, num_segments - 1 + joint_id) = cross(v, s - t)';
        
        %% shape
        v = T * [0; 1; 0; 1]; v = v(1:3) / v(4);
        v = v - t;
        
        c = 1;
        if l == 3 || segment_kinematic_chain(l + 1) == -1
            c = norm(s - t) / beta(l);
        end
        
        j(:, segment_id) = c * v;
    end
    
    %% "anonynous" function
    l_last = 0;
    if segment_kinematic_chain(3) > 0, l_last = 3; end
    if segment_kinematic_chain(3) == -1, l_last = 2; end
    if segment_kinematic_chain(2) == -1, l_last = 1; end
    
    fraction = norm(s - t) / beta(l_last); % change
    x = [beta; theta];
    offset  = m - s;
    m_ = projection_function(x, d, segment_kinematic_chain, JointsAxis, JointsSegmentId, segments, fraction, offset);
    df = my_gradient(@(x) projection_function(x, d, segment_kinematic_chain, JointsAxis, JointsSegmentId, segments, fraction, offset), x);
    
    %% jacobian convsegment
    variables = {'c2'};
    if (l_last == 1), i1 = 1; i2 = 2; end
    if (l_last == 2), i1 = 2; i2 = 3; end
    if (l_last == 3), i1 = 3; i2 = 4; end
    T1 = reshape(SegmentsGlobal(i1, :), 4, 4);
    T2 = reshape(SegmentsGlobal(i2, :), 4, 4);
    c1 = T1(1:3, 4);
    c2 = T2(1:3, 4);
    r1 = radii(i1);
    r2 = radii(i2);
    [m_, dm] = jacobian_convsegment(m, c1, c2, r1, r2, variables);
    df_ = dm.dc2 * v;
    %myline(m, m + df_, 'b');
    %disp(df_./(m - s));

    %% two projections    
    v1 = c1 + r1 * (m - s) / norm(m - s);
    c3 = c2 + v;
    z = c1 + (c3 - c1) * r1 / (r1 - r2);
    
    t = intersect_line_line(s, m, v1, z);
    
    dd = t - m;
    %disp([df_, dd]);
    disp(df_./dd);
    
    [~, q1, s, ~] = projection_convsegment(m, c1, c2, r1, r2, 1, 2);
    [~, q2, s, ~] = projection_convsegment(m, c1, c2 + 1e-10 * v, r1, r2, 1, 2);
    ddd = (q2 - q1)./1e-10;
end



