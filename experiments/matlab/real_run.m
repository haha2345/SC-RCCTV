% Real-data run: env REAL_DS in {indian220, urban200}, REAL_METH in
% {RCCTV, SC, SC_noT, SC_noG, SC_scalar, E3DTV, LRTDTV, tCTV}.  Frozen parameters, nothing tuned on real data.
R = setup_paths();
ds = getenv('REAL_DS'); meth = getenv('REAL_METH');
out_dir = fullfile(res_dirs(R),'real'); if ~exist(out_dir,'dir'), mkdir(out_dir); end
Y = load_real(R, ds); [M,N,p] = size(Y);
rk = 8; lam = 6/sqrt(M*N);
t0 = tic; info = struct();
switch meth
  case 'RCCTV',     X = RCCTV_RPCA(Y, rk, lam);
  case 'SC',        [X,info] = sc_rcctv(Y, struct('rk',rk));
  case 'SC_noT',    [X,info] = sc_rcctv(Y, struct('rk',rk,'T_mode','none'));
  case 'SC_noG',    [X,info] = sc_rcctv(Y, struct('rk',rk,'G_mode','none','sig_online',false));
  case 'SC_scalar', [X,info] = sc_rcctv(Y, struct('rk',rk,'scalar_sigma',true));
  case 'SC_fix',    [X,info] = sc_rcctv(Y, struct('rk',rk,'sig_online',false));   % sigma_k from the input only (two-stage)
  case 'E3DTV',     rng(0); X = EnhancedTV(Y, 0.004*sqrt(M*N), [13 13 13]);
  case 'LRTDTV',    rng(0); X = LRTDTV(Y, 1, 10, [round(0.8*M) round(0.8*N) 10]);
  case 'tCTV',      randn('state',0); rand('state',0); opts = struct('directions',[1 2]); X = TCTV_TRPCA(Y, opts);
  otherwise, error('unknown method');
end
sec = toc(t0); X = reshape(X,[M,N,p]);
save(fullfile(out_dir,[ds '_' meth '.mat']), 'X', 'sec', 'info', '-v7');
fprintf('REAL %s %s done in %.1fs\n', ds, meth, sec);

function Y = load_real(R, ds)
switch ds
  case 'indian220'
    S = load(fullfile(R,'data','real','IndianPines','Indian_pines.mat')); Y = double(S.indian_pines);
  case 'urban200'
    S = load(fullfile(R,'data','real','Urban','Urban_F210.mat')); Y = double(S.Urban_F210(54:253, 54:253, :));
end
for b = 1:size(Y,3)                                   % per-band min-max to [0,1] (as in the old R protocol)
    a = Y(:,:,b); mn = min(a(:)); mx = max(a(:)); Y(:,:,b) = (a-mn)/max(mx-mn, eps);
end
end
