clc;
clear;
n = 3;
x = randn(n, 1);
A = randn(n, n);
B = randn(n, n);


%% Hessian of vector function
f = @(x) [x' * A * x; x' * B * x];

df = @(x) [x' * (A + A'); x' * (B + B')];

df_num = my_gradient(f, x);
%disp(df_num); disp(df(x));

ddf{1} = @(x) A + A';
ddf{2} = @(x) B + B';

% ddf_num = my_gradient(df, x);
% for i = 1:size(ddf_num, 1)
%     disp(['f(', num2str(i), ')'])
%     disp(squeeze(ddf_num(i, :, :))); disp(ddf{i}(x));
% end

%% Hessian of dot product
y = @(x) x' * A * x;
dy = @(x) x' * (A + A');
ddy = @(x) A + A';
f = @(x) y(x)' * y(x);
df = @(x) 2 * y(x)' * dy(x);
ddf = @(x) 2 * dy(x)' * dy(x) + 2 * y(x)' * ddy(x);

ddf_num = my_gradient(df, x);
disp(ddf_num);
disp(ddf(x));

%dE_dX_ = @(X, N) 2 * [F1(X, N); sqrt(w2) * F2(X, N)]' * [dF1(X, N) * eye(N, N - 1); sqrt(w2) * dF2(X, N) * eye(N, N - 1)];

            %dX_ = -  dE1_dX_(X_) \ E2(X_, x);
            %lambda = 1;   
            %X_init = X_;
            %X_opt = X_ + lambda * dX_;
            %while (E2(X_, x) <  E2(X_opt, x))
            %    lambda = lambda / 2;                
            %    X_opt = X_ + lambda * dX_;
            %end
            %disp([E2(X_, x), E2(X_opt, x)]);
            %energy_lin = E2(X_, x) +  dE1_dX_(X_) * (lambda * dX_);
            %dX_value = dX_;