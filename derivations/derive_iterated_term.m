%close all;
clc; clear;
X = linspace(0, 1.5);
% f_ = @(x) exp(-1.3 * x).^2;
% j_ = @(x) -1.3 * exp(-1.3 * x);
% F = exp(-1.3 * X).^2;

f_ = @(x) (x - 2).^4;
j_ = @(x) 4 * (x - 2)^3;
F = (X - 2).^4;


x0 = 0.1;
x = x0;

line_colors = {[0.3, 0.8, 1.0], [1, 0.6, 0.1]};
%point_colors = {[0.7, 0.1, 0.6], [1, 0.4, 0.1]};
point_colors = {[0.3, 0.8, 1.0], [1, 0.6, 0.1]};
%%{
figure; hold on; plot(X, F, 'lineWidth', 3); mypoint([x0, f_(x0)], [0, 0.3, 0.7], 50);
lambdas = [0; 6];
for c = 1:length(lambdas)  
    lambda = lambdas(c);
    x = x0;
    for iter = 1:10
        f = f_(x);
        df = j_(x);
        
        %delta =  (x0 - x) + (df * (f - df * (x0 - x))) / (df * df + lambda);
        delta = (df * f) / (df * df + lambda);
        
        delta = - delta;
        myline([x, f], [x + delta, f + df * delta], line_colors{c});
        myline([x, f], [x, f + df * delta], line_colors{c});
        myline([x, f + df * delta], [x + delta, f + df * delta], line_colors{c});
        x = x + delta;
        myline([x, f + df * delta], [x, f_(x)], line_colors{c});
        mypoint([x, f_(x)], point_colors{c}, 50);
    end
end
%ylim([-0.3, 1]); xlim([0, 1.5]);
%%}

figure; hold on; plot(X, F, 'lineWidth', 3); mypoint([x0, f_(x0)], [0, 0.3, 0.7], 50); 
lambdas = [0, 1];
for c = 1:length(lambdas)  
    lambda = lambdas(c);
    x = x0;
    for iter = 1:4
        f = f_(x);
        df = j_(x);       
               
        f0 = f + df * (x0 - x);       
        x_plus = x0 - df * f0 / (df * df + lambda);         
        myline([x; f], [x; f + df * (x_plus - x)], line_colors{c});        
        myline([x0; f + df * (x0 - x)], [x; f], line_colors{c});        
        myline([x; f], [x_plus, f + df * (x_plus - x)], line_colors{c});
        myline([x_plus, f + df * (x_plus - x)], [x_plus, f_(x_plus)], line_colors{c});
        myline([x_plus, f + df * (x_plus - x)], [x0, f + df * (x_plus - x)],  line_colors{c});
        mypoint([x_plus, f_(x_plus)], point_colors{c}, 50);
        
        x = x_plus;
    end
end
%ylim([-0.3, 1]);
xlim([0, 2]);
