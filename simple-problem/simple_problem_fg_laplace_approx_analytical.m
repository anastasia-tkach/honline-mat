function [F, J, h] = simple_problem_fg_quadratic_two_lsqnonlin_analytical(xx, x_0, y, t, h_, w2)

%% Data term
F1 = @(xx) exp(xx(2) * t)^2 - y;
dF1_dx1 = @(xx) 0;
dF1_dx2 = @(xx) 2 * t * exp(xx(2) * t)^2;

ddF1_dx1_dx1 = @(xx) 0;
ddF1_dx1_dx2 = @(xx) 0;
ddF1_dx2_dx2 = @(xx) 4 * t.^2 * exp(xx(2) * t).^2;

ddF1 = @(xx) shiftdim([ddF1_dx1_dx1(xx), ddF1_dx1_dx2(xx); ddF1_dx1_dx2(xx), ddF1_dx2_dx2(xx)], -1);

%v = my_gradient(F1, xx);
%disp([v, dF1(xx)]);
%vv = my_gradient(dF1, xx);
%disp([vv, shiftdim(ddF1(xx), 1)]);

%% Closeness term
F2 = @(xx) sqrt(w2) * (xx(2) - xx(1));
dF2_dx1 = @(xx)  - sqrt(w2);
dF2_dx2 = @(xx) sqrt(w2);

ddF2 = @(xx) shiftdim([0, 0; 0, 0], -1);

%% Quadratic approx
if isempty(x_0)
    Q = @(xx) 0;
    dQ_dx1 = @(xx) 0;
    dQ_dx2 = @(xx) 0;
    ddQ = @(xx) shiftdim([0, 0; 0, 0], -1);
else
    A = h_(1, 1);
    B = h_(1, 2);
    C = h_(2, 1);
    D = h_(2, 2);
    %ddQ_dx2_dx2 = D - C * inv(A) * B;
    ddQ_dx2_dx2 = D;
    ddQ_dx2_dx2_sqrt = sqrt(ddQ_dx2_dx2);
    Q = @(xx) ddQ_dx2_dx2_sqrt * (xx(1) - x_0);
    dQ_dx1 = @(xx) ddQ_dx2_dx2_sqrt;
    dQ_dx2 = @(xx) 0;
    ddQ = @(xx) shiftdim([0, 0; 0, 0], -1);
end

F = @(xx) [Q(xx); F1(xx); F2(xx)];
J = @(xx) [dQ_dx1(xx), dQ_dx2(xx); dF1_dx1(xx), dF1_dx2(xx); dF2_dx1(xx), dF2_dx2(xx)];
H = @(xx) [ddQ(xx); ddF1(xx); ddF2(xx)];

%% Hessian for scalar objective
f = @(xx) F(xx)' * F(xx);    
j = @(xx) 2 * F(xx)' * J(xx);    
h = @(xx) hessian_for_scalar_objective(F(xx), J(xx), H(xx));

%v = my_gradient(f, xx);
%disp([v; j(xx)]);
%vv = my_gradient(j, xx);
%disp([vv, h(xx)]);

F = F(xx);
J = J(xx);
H = H(xx);

f = f(xx);
j = j(xx);
h = h(xx);
