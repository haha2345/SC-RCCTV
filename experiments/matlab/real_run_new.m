% Real-data runs of the revision baselines (same loading/normalisation as real_run.m; nothing tuned on real data).
% env REAL_DS in {indian220, urban200}, REAL_METH in {FastHyMix, L1HyMixDe, tCTVd}, TCTV_MULT.
R = setup_paths_new();
warning('off','MATLAB:mpath:nameNonexistentOrNotADirectory');
ds = getenv('REAL_DS'); meth = getenv('REAL_METH');
tm = str2double(getenv('TCTV_MULT')); if isnan(tm), tm = 1; end
out_dir = fullfile(res_dirs(R),'real'); if ~exist(out_dir,'dir'), mkdir(out_dir); end
switch ds
  case 'indian220', S = load(fullfile(R,'data','real','IndianPines','Indian_pines.mat')); Y = double(S.indian_pines);
  case 'urban200',  S = load(fullfile(R,'data','real','Urban','Urban_F210.mat')); Y = double(S.Urban_F210(54:253, 54:253, :));
end
for b = 1:size(Y,3), a = Y(:,:,b); mn = min(a(:)); mx = max(a(:)); Y(:,:,b) = (a-mn)/max(mx-mn, eps); end
t0 = tic; [X, info] = run_new_method(meth, Y, struct('tctv_mult', tm)); sec = toc(t0);
save(fullfile(out_dir,[ds '_' meth '.mat']), 'X', 'sec', 'info', '-v7');
fprintf('REAL %s %s done in %.1fs\n', ds, meth, sec);
