%% Estimate Moments Using Independent Metropolis-Hastings Sampling
%

% Copyright 2015 The MathWorks, Inc.


%%
% Use Independent Metropolis-Hastings sampling to estimate the second
% order moment of a Gamma distribution.
rng default;  
alpha = 2.43;
beta = 1;
pdf = @(x)gampdf(x,alpha,beta); % Target distribution
proppdf = @(x,y)gampdf(x,floor(alpha),floor(alpha)/alpha);
proprnd = @(x)sum(exprnd(floor(alpha)/alpha,floor(alpha),1));
nsamples = 5000;
smpl = mhsample(1,nsamples,'pdf',pdf,'proprnd',proprnd, 'proppdf',proppdf);
            
%%
% Plot the results.
xxhat = cumsum(smpl.^2)./(1:nsamples)';
figure;
plot(1:nsamples,xxhat)

