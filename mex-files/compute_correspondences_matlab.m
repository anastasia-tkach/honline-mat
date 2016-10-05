function [I, M] = compute_correspondences_matlab(C, B, D)

RAND_MAX = 32767;
num_points = size(D, 1);

M = zeros(num_points, 2);
I = zeros(num_points, 1);

for i = 1:num_points    
    p = D(i, :)';   
    min_distance = RAND_MAX;
    for j = 1:size(B, 1)
        
        block = B(j, :)';  
        c1 = C(block(1), :)'; 
        c2 = C(block(2), :)'; 
        q = projection_segment_cpp(p, c1, c2);
       
        if norm(p - q) < min_distance
            min_distance = norm(p - q);
            M(i, :) = q';
            I(i) = j;
        end
    end
    
end


