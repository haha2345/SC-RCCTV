function [X, info] = run_new_method(m, Y, o)
% Baselines added in the revision.  o: struct with (optional) oracle noise parameters
%   o.sigma (scalar), o.sp_rate, o.st_rate  -- used only by S3TTV (constrained model needs them)
%   o.tctv_mult  -- lambda multiplier of t-CTV (chosen on the dev cubes)
info = struct('iters', NaN);
[M,N,p] = size(Y);
switch m
  case 'FastHyMix'                      % parameter-free; k_subspace = 8 (authors' default)
    X = FastHyMix(Y, 8);
  case 'L1HyMixDe'                      % authors' defaults: k_subspace = 8, p = 0.1
    X = L1HyMixDe(Y, 8, 0.1);
  case 'S3TTV'                          % authors' settings (JSTARS 2025, Sec. IV): rho = 0.95, 0.98 in cases 1 and 4,
    Nn = numel(Y); rho = 0.95;          % epsilon fixed to 0.01 in case 4; blocksize 10; stopcri 1e-5
    c = gf(o,'case',0); if any(c == [1 4]), rho = 0.98; end
    prm.epsilon = rho * o.sigma * sqrt(Nn * (1 - o.sp_rate));
    if c == 4, prm.epsilon = 0.01; end
    prm.alpha   = rho * (0.5 * Nn * o.sp_rate);
    prm.beta    = rho * Nn * (1 - o.sp_rate) * o.st_rate * 0.5 / 2;
    prm.blocksize = [10 10]; prm.stopcri = 1e-5; prm.maxiter = gf(o,'maxiter',20000);
    if isfield(o,'trace_clean'), prm.trace_clean = o.trace_clean; prm.trace_every = 100; end   % logging only
    prm.use_gpu = gf(o,'use_gpu',has_gpu());                 % env S3TTV_GPU=0/1 overrides (sim_run_new)
    [X, it] = s3ttv_cpu(Y, prm); info.iters = it;
  case 'tCTVd'                          % t-CTV with dev-tuned lambda
    randn('state',0); rand('state',0);
    [X,~,~,~,it] = tctv_trpca_lam(Y, struct('directions',[1 2],'detail',1,'lambda_mult',o.tctv_mult)); info.iters = it;
  otherwise
    error('unknown method %s', m);
end
X = reshape(double(X), [M,N,p]);
end

function g = has_gpu()
try, g = gpuDeviceCount("available") > 0; catch, g = false; end
end

function v = gf(o,f,d)
if isfield(o,f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
