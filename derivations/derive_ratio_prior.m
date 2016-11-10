clc; clear;
n = 3;
x = rand(n, 1);

q = [2; 4; 1]; %vector of ratios of the model parameters, you should fill it in yourself, using the template model
S = ones(1, n);
A = inv(S * q) * q * S;
%{
the idea of this prior is to find a vector with the same sum, but with a
given ratios of the elements.
this is ok for the finger, but for the palm we might want something
different.

the question is: given the current hand model, how can we compute a scaling
factor for a best matching template model? What do you think? 
We may as well say that they should be the same heights
%}

%% the objective function 
F = @(x) (eye(n, n) -  A) * x;
options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt');
x_opt = lsqnonlin(F, x, [], [], options);

disp(x);
disp(x_opt);




