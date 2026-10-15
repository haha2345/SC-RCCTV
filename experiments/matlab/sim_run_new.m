% Native-size simulation driver for the revision baselines (same inputs/scoring as sim_run.m).
% env: SIM_CUBE, SIM_METHS, SIM_SEEDS, SIM_CASES, SIM_TAG, [TCTV_MULT]
% S3TTV gets the oracle noise parameters its constrained model needs (true sigma; rms of sigma_k in
% cases 9-10; true impulse and stripe rates) -- an advantage the blind methods do not have.
R = setup_paths_new();
warning('off','MATLAB:mpath:nameNonexistentOrNotADirectory');
cube  = getenv('SIM_CUBE');
meths = strsplit(getenv('SIM_METHS'), ',');
seeds = str2double(strsplit(getenv('SIM_SEEDS'), ','));
cases = str2double(strsplit(getenv('SIM_CASES'), ','));
tag   = getenv('SIM_TAG');
tm    = str2double(getenv('TCTV_MULT')); if isnan(tm), tm = 1; end
SP = [0 .05 .05 0 0 0 .05 .05 0 .05];  ST = [0 0 0 1 1 1 1 1 0 1]*0.05;
[base, inb] = res_dirs(R);
csv   = fullfile(base, ['sim_' tag '.csv']);
C = load(fullfile(R,'data','locked_nhsi',cube,'clean.mat')); clean = C.clean; [M,N,p] = size(clean);
done = {};
if exist(csv,'file'), done = strsplit(fileread(csv), newline); else, fid = fopen(csv,'w'); fprintf(fid,'cube,case,seed,meth,mpsnr,mssim,sam,ergas,sec,iters\n'); fclose(fid); end
for c = cases
  for s = seeds
    f = fullfile(inb,'inputs',sprintf('%s_c%d_s%d.mat',cube,c,s));
    if ~exist(f,'file'), continue; end
    I = load(f); Y = I.Nhsi;
    o = struct('sigma', sqrt(mean(I.sigma_true(:).^2)), 'sp_rate', SP(c), 'st_rate', ST(c), 'tctv_mult', tm, 'case', c);
    s3max = str2double(getenv('S3TTV_MAXIT')); if ~isnan(s3max), o.maxiter = s3max; o.trace_clean = clean; end
    s3gpu = str2double(getenv('S3TTV_GPU'));   if ~isnan(s3gpu), o.use_gpu = s3gpu > 0; end
    for mi = 1:numel(meths)
      m = meths{mi};
      key = sprintf('%s,%d,%d,%s,', cube, c, s, m);
      if any(strncmp(done, key, numel(key))), continue; end
      t0 = tic;
      try
        [X, info] = run_new_method(m, Y, o);
      catch ME
        fprintf('ERROR %s c%d s%d %s : %s\n', cube, c, s, m, ME.message); continue;
      end
      sec = toc(t0);
      [ps,ss,sam,erg] = score_full(X, clean);
      fid = fopen(csv,'a'); fprintf(fid,'%s,%d,%d,%s,%.6f,%.6f,%.6f,%.6f,%.2f,%g\n', cube,c,s,m,ps,ss,sam,erg,sec,info.iters); fclose(fid);
      fprintf('%s c%d s%d %-9s mpsnr %.3f  %.1fs\n', cube,c,s,m,ps,sec);
    end
  end
end
fprintf('SIM DONE %s %s\n', cube, tag);
