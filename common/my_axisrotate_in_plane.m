function [R] = my_axisrotate_in_plane(theta)

c = cos(theta);
s = sin(theta);

R = eye(4, 4);

R(1, 1) = c;
R(1, 2) = -s;

R(2, 1) = s;
R(2, 2) = c;
return;

R(1, 1) = t*x*x + c;
R(1, 2) = t*x*y - z*s;
R(1, 3) = t*x*z + y*s;

R(2, 1) = t*x*y + z*s;
R(2, 2) = t*y*y + c;
R(2, 3) = t*y*z - x*s;

R(3, 1) = t*x*z - y*s;
R(3, 2) = t*y*z + x*s;
R(3, 3) = t*z*z + c;

