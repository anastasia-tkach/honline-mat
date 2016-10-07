function str = au_autodiff_generate(function_handle, example_arguments, example_data, filename, varargin)
% AU_AUTODIFF_GENERATE  Generate code for function and derivatives
%      To generate C code for a function and its jacobian, use
%       STR = AU_AUTODIFF_GENERATE(FUNCTION_HANDLE, EXAMPLE_ARG, EXAMPLE_DATA)
%      or 
%       AU_AUTODIFF_GENERATE(FUNCTION_HANDLE, EXAMPLE_ARG, EXAMPLE_DATA, FILE)
%      The first form returns the core code in a string, the second wraps
%      it in a MEX wrapper, and stores it in file FILE.
%      The EXAMPLE_ARG and EXAMPLE_DATA are used to determine the size, and
%      to perform finite-difference checks just in case...
%
%      FUNCTION_HANDLE must take a column vector and return a column vector.  
%      It can also take an optional argument DATA.

% awf, dec13

if nargin == 0
    %%
    f = @(x,data) [sin(x(3)/norm(x(1:2))); 14*cos((6*x(1)-6*x(3))+x(2)^2)];
    p = [1 2 3]';
    f(p)
    
    au_test_begin au_autodiff_generate
    au_autodiff_generate(f, p, [], 'C:/Users/tkach/Desktop/au_autodiff_generate_test_mex.cpp', 'HESSIAN=1');
 
    % Check it exists and is callable
    for do_jacobian = [0 1 2]
      out = au_autodiff_generate_test_mex([p p p], zeros(0,3), do_jacobian);
      out_val = out(1,:);
      out_jac = out(2:end,:);
    end
    
    au_test_should_fail o=au_autodiff_generate_test_mex(rand(4,1),zeros(0,1),0)
    au_test_should_fail o=au_autodiff_generate_test_mex(rand(3,1),zeros(2,1),0)

    au_test_end
    
    return
end

if nargin < 4
    filename = [];
end

opts = au_opts('CSE=1;JACOBIAN=1;HESSIAN=1;SIMPLIFY_TMAX=0', varargin{:});

%% Determine sizes etc
%au_assert size(example_arguments,2)==1
out = function_handle(example_arguments, example_data);
%au_assert size(out,2)==1

m = size(example_arguments,1);
n = size(out,1);
md = size(example_data,1);

%% Make symbolic variables and push them through.
in = sym('x', [m 1]);
assume(in, 'real');
data = sym('data', [md 1]);
assume(data, 'real');

tic
fprintf('au_autodiff_generate: making code for f:R^%d->R^%d ', m, n);
out_val = function_handle(in, data);
fprintf(' [%.1fsec]\n', toc);
if opts.JACOBIAN
  fprintf('au_autodiff_generate: computing jacobian ... ');
  times=[];
  for k=numel(in):-1:1
    tic
    out_jac(:,k) = diff(out_val, in(k));
    times(k) = toc;
    fprintf(' %d [%.1fsec]', k, times(k));
  end
  fprintf(' %.1fsec\n', sum(times));
  fprintf('au_autodiff_generate: transpose ... ');
  tic
  out_jac = [out_jac];
  fprintf(' %.1fsec\n', toc);
end

if opts.HESSIAN
  fprintf('au_autodiff_generate: computing Hessian ... ');
  times=[];
  nh=1;
  for k=numel(in):-1:1
    tic
    out_hess(:,nh:nh+k-1) = diff(out_jac(:,1:k), in(k));
    times(k) = toc;
    fprintf(' %d [%.1fsec]', k, times(k));
    nh=nh+k;
  end
  fprintf(' %.1fsec\n', sum(times));
  fprintf('au_autodiff_generate: transpose ... ');
  tic
  out_hess = [out_jac out_hess];
  fprintf(' %.1fsec\n', toc);
end


% awf xxfixme: handle symbolic sparsity find(~(out_all == 0))

%% Now wrap the computation in au_autodiff_generate_template.cpp
decls = '';
for k=1:m
    decls= [decls sprintf('    double x%d = in_ptr[%d];\n', k, k-1)];
end
for k=1:md
    decls= [decls sprintf('    double data%d = data_ptr[%d];\n', k, k-1)];
end
fprintf('au_autodiff_generate: computing ccode for do_jacobian=0\n');
BodyNoJac = [decls au_ccode(out_val.', [], opts.CSE, opts.SIMPLIFY_TMAX)];
BodyNoJac = strrep(BodyNoJac, '[0 * out_rows + ', '[');
str = BodyNoJac;
do_jacobian = 0;

if opts.JACOBIAN
  tic
  fprintf('au_autodiff_generate: computing ccode for do_jacobian=1\n');
  BodyJacobian = [decls au_ccode(out_jac.', [], opts.CSE, opts.SIMPLIFY_TMAX)];
  BodyJacobian = strrep(BodyJacobian, '[0 * out_rows + ', '[');
  codegen_time = toc;
  str = BodyJacobian;
  do_jacobian = 1;
end

if opts.HESSIAN
  tic
  fprintf('au_autodiff_generate: computing ccode for do_jacobian=2\n');
  BodyHessian = [decls au_ccode(out_hess.', [], opts.CSE, opts.SIMPLIFY_TMAX)];
  BodyHessian = strrep(BodyHessian, '[0 * out_rows + ', '[');
  codegen_time = toc;
  str = BodyHessian;
  do_jacobian = 2;
end

if isempty(filename)
    return
end

%% If there is a filename, make it a mexFunction

% File descriptor or name?
if ischar(filename)
    fprintf('au_autodiff_generate: writing file [%s]\n', filename);
    fd = fopen(filename, 'w');
else
    fd = filename;
    fprintf('au_autodiff_generate: writing file [fd=%d]\n', fd);
end

body0 = sprintf('  /* inner loop do_jac=0 */\n%s', BodyNoJac);
if opts.JACOBIAN
  body1 = sprintf('  /* inner loop do_jac=1 */\n%s', BodyJacobian);
else
  body1 = 'mexErrMsgTxt("nojac");';
end
if opts.HESSIAN
  body2 = sprintf('  /* inner loop do_jac=2 */\n%s', BodyHessian);
else
  body2 = 'mexErrMsgTxt("nohess");';
end

% Get Template Text
tfd = fopen([au_root_dir '/au_autodiff_generate_template.cpp'], 'r');
template = fread(tfd, inf, 'char');
fclose(tfd);

varname = inputname(1);
if isempty(varname)
    varname = '[Anonymous expression]';
end

template = char(template');
template = strrep(template, '@VarName', varname);
template = strrep(template, '@BodyNoJac', body0);
template = strrep(template, '@BodyJacobian', body1);
template = strrep(template, '@BodyHessian', body2);
template = strrep(template, '@InRows', num2str(m));
template = strrep(template, '@OutDim', num2str(n));
template = strrep(template, '@DataRows', num2str(md));
template = strrep(template, '@OutRows', num2str(size(out,1)));

fprintf(fd, '%s', template);
if ischar(filename)
    fclose(fd);
end

%% And mex it and test it...
if ischar(filename)
    fprintf('au_autodiff_generate: mexing... [%s]\n', filename)
    mex(['-I' au_root_dir], filename);

    [~,fn,~] = fileparts(filename);
        
    % Check we're not calling an old one on the path...
    d = dir(which(fn)); 
    mex_file_age_in_seconds = (now - datenum(d.date))*24*60*60;
    au_test_assert mex_file_age_in_seconds<5
    
    fprintf('au_autodiff_generate: testing... [%s]\n', fn)
    val_f = function_handle(example_arguments, example_data);
    val_mex = feval(fn, example_arguments(:), example_data(:), do_jacobian);
    
    val_f = val_f(:)';
    %au_test_equal(val_f, val_mex(1,:), 1e-8);
    disp([val_f, val_mex(1,:)]);

    if opts.JACOBIAN
      J = val_mex(2:(1+m),:)';
      f_dc = @(x) feval(fn, x, example_data(:), 0)';
      timeout = max(2, codegen_time/10);
      au_check_derivatives(f_dc, example_arguments(:), J, 'delta=1e-5;tol=1e-3', 'timeout', timeout);
    end
end

