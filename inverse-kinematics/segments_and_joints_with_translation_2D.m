function [segments, joints] = segments_and_joints_with_translation_2D()

%% Segments
segments{1}.kinematic_chain = 1:3;
segments{1}.length = 5;
segments{1}.parent_id = 0;
segments{1}.global = eye(4, 4);
segments{1}.local = eye(4, 4);

segments{2}.kinematic_chain = 1:4;
segments{2}.length = 4;
segments{2}.parent_id = 1;
segments{2}.local = eye(4, 4);
segments{2}.local(2, 4) = segments{segments{2}.parent_id}.length;

segments{3}.kinematic_chain = 1:5;
segments{3}.length = 3;
segments{3}.parent_id = 2;
segments{3}.local = eye(4, 4);
segments{3}.local(2, 4) = segments{segments{3}.parent_id}.length;

segments{4}.kimenatic_chain = 1:6;
segments{4}.length = 0;
segments{4}.parent_id = 3;
segments{4}.local = eye(4, 4);
segments{4}.local(2, 4) = segments{segments{4}.parent_id}.length;

%% Joints

joints{1}.segment_id = 1; 
joints{1}.axis = [1; 0; 0]; 
joints{1}.type = 'T';

joints{2}.segment_id = 1;
joints{2}.axis = [0; 1; 0]; 
joints{2}.type = 'T';

joints{3}.segment_id = 1;
joints{3}.axis = [0; 0; 1]; 
joints{3}.type = 'R';

joints{4}.segment_id = 2;
joints{4}.axis = [0; 0; 1];
joints{4}.type = 'R';

joints{5}.segment_id = 3;
joints{5}.axis = [0; 0; 1]; 
joints{5}.type = 'R';

joints{6}.segment_id = 4;
joints{6}.axis = [0; 0; 1]; 
joints{6}.type = 'R';
