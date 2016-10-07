clear; clc;
m = 4;
x = rand(m, 1);
A = rand(m, m);
B = rand(1, m);
c = rand(1, 1);
D = [A(:); B'; c];
au_autodiff_generate(@test_function, x, D, 'E:/OneDrive/EPFL/Code/honline-mat/test_autodiff/au_autodiff_test.cpp', 'HESSIAN=1');

[V] = au_autodiff_test(x, D, 2);
j_ = V(1:m);
h_ = V(m + 1:end - 1);

J_ = j_';
H_ = zeros(m, m);
k = 0;
for  i = m:-1:1
    H_(1:i, i) = h_(k + 1:k + i);
    H_(i, 1:i-1) =  h_(k + 1:k + i - 1);
    k = k + i;
end

[F, J, H] = test_function(x, D);


