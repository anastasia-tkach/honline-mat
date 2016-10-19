function [] = display_posterior_and_laplace_approx(xx_opt, F, J, h, x_0, y, t, h_, w2)

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
        
        p_approx(u, v) = exp(-f) * exp(-A);
        %E_approx(u, v) = - log(p_approx(u, v));
        E_approx(u, v) = f + A;
        
        p_mygauss(u, v) = K * exp(-A);
        %E_mygauss(u, v) = -log(p_mygauss(u, v));
        E_mygauss(u, v) = -log(K) + A;
        
    end
        
end

p_gauss = mvnpdf([X1(:) X2(:)], xx_opt', sigma);
p_gauss = reshape(p_gauss, length(x2), length(x1));
E_gauss = -log(p_gauss);

E_mygauss = E_mygauss + log(K) + f;
E_gauss = E_gauss + log(K) + f;

figure; hold on;

surf(X1, X2, E);
mesh(X1, X2, E_approx); 
mesh(X1, X2, E_gauss); 
mesh(X1, X2, E_mygauss);

mypoint([xx_opt; f], 'm', 20);
view([-45, 30]);
