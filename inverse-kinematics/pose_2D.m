function [segments] = pose_2D(segments, joints, theta)

theta = [theta; 0];

%% Pose segments
for i = 1:length(joints)
    segment = segments{joints{i}.segment_id};
    T = [];
    switch joints{i}.type
        case 'R'
            %T = segment.local * makehgtform('axisrotate', joints{i}.axis, theta(i));
            T = segment.local * my_axisrotate_in_plane(theta(i));            
        case 'T'
            T = segment.local * makehgtform('translate', joints{i}.axis * theta(i));
    end    
    if segment.parent_id > 0
        segment.local = T;
        segment.global = segments{segment.parent_id}.global * T;
    else
        segment.local = T;
        segment.global = T;
    end
    segments{joints{i}.segment_id} = segment;
end