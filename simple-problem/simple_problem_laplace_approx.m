function [xx_opt, J, h] = simple_problem_laplace_approx(X0, X_prev, h_, Y, T, N, w2, num_iters)

y_ = Y(N - 1);
t_ = T(N - 1);
y = Y(N);
t = T(N);
x_0 = X_prev(N - 1);
x__0 = X_prev(N - 2);

%{
f1 =  exp(x_0 * t_)^2 - y_;
j1 = 2 * t_ * exp(x_0 * t_)^2;
h1 = 4 * t_.^2 * exp(x_0 * t_).^2;

f2 = sqrt(w2) * (x_0 - x__0);
j2 = sqrt(w2) * 1;
h2 = 0;

E1_ = f1' * f1; E2_ = f2' * f2;
dE1_ = 2 * f1' * j1; dE2_ = 2 * f2' * j2;
ddE1_ = 2 * j1' * j1 + 2 * f1 * h1; 
ddE2_ = 2 * j2' * j2 + 2 * f2 * h2;
E_ = E1_ + E2_; dE_ = dE1_ + dE2_; ddE_ = ddE1_ + ddE2_;
h_ = [0, 0; 0, ddE_];
%}

%%{
if N == 3
    [F_, J_, h_] = simple_problem_fg_laplace_approx(X_prev(N - 2: N - 1), [], y_, t_, h_, w2);
end
%%}

%% Quadratic approximation
xx0 = X0(N-1:N);

[xx_opt] = my_lsqnonlin(@(xx) simple_problem_fg_laplace_approx(xx, x_0, y, t, h_, w2), xx0, num_iters);
[F, J, h] = simple_problem_fg_laplace_approx(xx_opt, x_0, y, t, h_, w2);

%% Plot energy and laplace approx
f = F' * F;  
j = 2 * F' * J;  
sigma = inv(h);
K = 1/(2*pi)/sqrt(det(sigma));

x1 = xx_opt(1) - 0.3:0.02:xx_opt(1) + 0.3;
x2 = xx_opt(2) - 0.3:0.02:xx_opt(2) + 0.3;
[X1, X2] = meshgrid(x1, x2);
E = zeros(length(x1), length(x2));
E_approx = zeros(length(x1), length(x2));
E_mygauss = zeros(length(x1), length(x2));
p_mygauss = zeros(length(x1), length(x2));
p_approx = zeros(length(x1), length(x2));
for u = 1:length(x1)
    for v = 1:length(x2)
        xx = [X1(u, v); X2(u, v)];
        [F, ~, ~] = simple_problem_fg_laplace_approx(xx, x_0, y, t, h_, w2);
        E(u, v) = F' * F;
        
        A = 0.5 * (xx - xx_opt)' * h *  (xx - xx_opt);        
        
        p_approx(u, v) = exp(f) * exp(A);
        E_approx(u, v) = log(p_approx(u, v));
        
        E_mygauss(u, v) = -log(K) + A;
        p_mygauss(u, v) = K * exp(-A);
    end
        
end

mu = xx_opt;

p_gauss = mvnpdf([X1(:) X2(:)], mu', sigma);
p_gauss = reshape(p_gauss, length(x2), length(x1));
E_gauss = -log(p_gauss);

E_mygauss = E_mygauss + log(K) - f;
E_gauss = E_gauss + log(K) - f;

%% 
%p_mygauss = p_mygauss .* exp(f) ./ K; 

figure; hold on;
%mesh(X1, X2, p_approx);
% surf(X1, X2, p_gauss);
% surf(X1, X2, p_mygauss);

% mesh(X1, X2, E);
surf(X1, X2, E_approx);
mesh(X1, X2, E_gauss);
mesh(X1, X2, E_mygauss);


%mypoint(xx_opt, 'm', 20);
view([-45, 30]);
disp(' ');

