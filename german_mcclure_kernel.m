function [f, df] = german_mcclure_kernel(r, dr)

e = 1.5;
B = 3;

rho = (r/e)'*(r/e)/(1 + (r/e)' * (r/e));
drho = 2 * r'/e * dr/e / (1 + (r/e)' * (r/e))^2;

a = sqrt(rho);
da = 0.5 * drho / sqrt(rho);
b = norm(r);
db = r' * dr / norm(r);
c = a/b;
dc = (da * b - a * db) / (b' * b);

f = c * r;
f = sqrt(rho)/norm(r) * r;
df = r * dc + c * dr;

%% Verify analytically
return
r = @(x) x_prev - x; %r = X(i) - X(i + 1);
dr = @(x) - eye(B, B);
rho = @(x) (r(x)/e)'*(r(x)/e)/(1 + (r(x)/e)' * (r(x)/e));
drho = @(x) 2 * r(x)'/e * dr(x)/e / (1 + (r(x)/e)' * (r(x)/e))^2;

a = @(x) sqrt(rho(x));
da = @(x) 0.5 * drho(x) / sqrt(rho(x));
b = @(x) norm(r(x));
db = @(x) r(x)' * dr(x) / norm(r(x));
c = @(x) a(x)/b(x);
dc = @(x) (da(x) * b(x) - a(x) * db(x)) / (b(x)' * b(x));

f = @(x) c(x) * r(x);
f = @(x) sqrt(rho(x))/norm(r(x)) * r(x);
df = @(x) r(x) * dc(x) + c(x) * dr(x);

v = my_gradient(f, x); disp([v; df(x)]);

