function [E, dE, Q, dQ, e1, de1, e2, de2] = simple_problem_fg_quadratic_two(xx, x_0, y, t, E_, dE_, ddE_, w2)

Q = E_ + dE_* (xx(1) - x_0) + (xx(1) - x_0)' * ddE_ * (xx(1) - x_0);
dQ_dx1 = dE_ +  2 * xx(1)' * ddE_ - 2 * x_0' * ddE_;
dQ_dx2 = 0;

f1 = exp(xx(2) * t)^2 - y;
df1_dx1 = 0;
df1_dx2 = 2 * t * exp(xx(2) * t)^2;

e1 = f1' * f1;
de1_dx1 = 2 * df1_dx1' * f1;
de1_dx2 = 2 * df1_dx2' * f1;

f2 = sqrt(w2) * (xx(2) - xx(1));
df2_dx1 = - sqrt(w2);
df2_dx2 = sqrt(w2);

e2 = f2' * f2;
de2_dx1 = 2 * df2_dx1' * f2;
de2_dx2 = 2 * df2_dx2' * f2;

E = e1 + e2 + Q;
dE_dx1 = de1_dx1 + de2_dx1 + dQ_dx1;
dE_dx2 = de1_dx2 + de2_dx2 + dQ_dx2;
dE = [dE_dx1, dE_dx2];

dQ = [dQ_dx1, dQ_dx2];
de1 = [de1_dx1, de1_dx2];
de2 = [de2_dx1, de2_dx2];

% disp(xx); disp(dE);