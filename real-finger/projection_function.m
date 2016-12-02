function [m] = projection_function(x, p, segment_kinematic_chain, JointsAxis, JointsSegmentId, segments, fraction, offset)

beta = x(1:3);
theta = x(4:6);

[segments] = shape_3D(segments, beta);

kinematic_chain_globals = cell(length(segment_kinematic_chain), 1);
for l = 1:length(segment_kinematic_chain)
    if segment_kinematic_chain(l) == -1, break; end
    joint_id = segment_kinematic_chain(l);
    segment_id = JointsSegmentId(joint_id);
    u = JointsAxis(joint_id, :)';
    
    %% compute transformations
    segment = segments{segment_id};
    T = segment.local * makehgtform('axisrotate', u, theta(l));
    if l > 1
        segment.global = kinematic_chain_globals{l - 1} * T;
    else
        segment.global = T;
    end
    kinematic_chain_globals{l} = segment.global;
    
    if l == 3 || segment_kinematic_chain(l + 1) == -1
        t =  kinematic_chain_globals{l}(1:3, 4);
        direction =  kinematic_chain_globals{l} * [0; 1; 0; 1];
        direction = direction(1:3)/direction(4) - t;
        s = t + fraction * beta(l) * direction(1:3);
        m = s + offset;
    end
    
end
