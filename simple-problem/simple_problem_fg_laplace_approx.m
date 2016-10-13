function [F, J, h] = simple_problem_fg_laplace_approx(xx, x_0, y, t, h_, w2)

%%{
%% Data term
F1 = exp(xx(2) * t)^2 - y;
dF1_dx1 = 0;
dF1_dx2 = 2 * t * exp(xx(2) * t)^2;

ddF1_dx1_dx1 = 0;
ddF1_dx1_dx2 = 0;
ddF1_dx2_dx2 = 4 * t.^2 * exp(xx(2) * t).^2;

ddF1 = shiftdim([ddF1_dx1_dx1, ddF1_dx1_dx2; ddF1_dx1_dx2, ddF1_dx2_dx2], -1);

%% Closeness term
F2 = sqrt(w2) * (xx(2) - xx(1));
dF2_dx1 =  - sqrt(w2);
dF2_dx2 = sqrt(w2);

ddF2 = shiftdim([0, 0; 0, 0], -1);

%% Quadratic approx
if isempty(x_0)
  Q = 0;
  dQ_dx1 = 0;
  dQ_dx2 = 0;
  ddQ = shiftdim([0, 0; 0, 0], -1);
else
  a = h_(1, 1);
  b = h_(1, 2);
  c = h_(2, 1);
  d = h_(2, 2);
  ddQ_dx2_dx2 = d - c * inv(a) * b;
  %ddQ_dx2_dx2 = D;
  ddQ_dx2_dx2_sqrt = sqrt(ddQ_dx2_dx2);
  Q = ddQ_dx2_dx2_sqrt * (xx(1) - x_0);
  dQ_dx1 = ddQ_dx2_dx2_sqrt;
  dQ_dx2 = 0;
  ddQ = shiftdim([0, 0; 0, 0], -1);
end

F = [Q; F1; F2];
J = [dQ_dx1, dQ_dx2; dF1_dx1, dF1_dx2; dF2_dx1, dF2_dx2];
H = [ddQ; ddF1; ddF2];

%% Hessian for scalar objective
f = F' * F;  
j = 2 * F' * J;  
h = hessian_for_scalar_objective(F, J, H);

