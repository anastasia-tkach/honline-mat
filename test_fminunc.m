n = 100;
xstart = -ones(n,1);
xstart(2:2:n,1) = 1;
options = optimoptions(@fminunc,'Algorithm','trust-region', 'SpecifyObjectiveGradient',true);
[x,fval,exitflag,output] = fminunc(@brownfgh,xstart,options);