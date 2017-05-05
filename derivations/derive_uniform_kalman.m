
T = 10;
x_hat = 0;
sigma_hat = 0;
for i = 1:T
    x = rand();
    sigma =  1/T;
    x_hat = sigma / (sigma_hat + sigma) * x_hat + sigma_hat / (sigma_hat + sigma) * x;
    sigma_hat = sigma_hat / (sigma_hat + sigma) * sigma;
end