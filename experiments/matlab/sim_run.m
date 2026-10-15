% Native-size simulation driver.  env: SIM_CUBE, SIM_METHS (comma list), SIM_SEEDS (e.g. 0,1,2,3,4,5), SIM_CASES, SIM_TAG
% Appends rows to results/sim_<tag>.csv ; resumes (skips finished rows).
R = setup_paths();
cube  = getenv('SIM_CUBE');
meths = strsplit(getenv('SIM_METHS'), ',');
seeds = str2double(strsplit(getenv('SIM_SEEDS'), ','));
cases = str2double(strsplit(getenv('SIM_CASES'), ','));
tag   = getenv('SIM_TAG');
[base, inb] = res_dirs(R);
csv   = fullfile(base, ['sim_' tag '.csv']);
C = load(fullfile(R,'data','locked_nhsi',cube,'clean.mat')); clean = C.clean; [M,N,p] = size(clean);
done = {};
if exist(csv,'file'), done = strsplit(fileread(csv), newline); else, fid = fopen(csv,'w'); fprintf(fid,'cube,case,seed,meth,mpsnr,mssim,sam,ergas,sec,iters\n'); fclose(fid); end
for c = cases
  for s = seeds
    f = fullfile(inb,'inputs',sprintf('%s_c%d_s%d.mat',cube,c,s));
    if ~exist(f,'file'), continue; end
    I = load(f); Y = I.Nhsi; sg = I.sigma_true; beta = I.beta_true;
    for mi = 1:numel(meths)
      m = meths{mi};
      key = sprintf('%s,%d,%d,%s,', cube, c, s, m);
      if any(strncmp(done, key, numel(key))), continue; end
      t0 = tic; iters = NaN;
      try
      switch m
        case 'RCCTV',     X = RCCTV_RPCA(Y, 8, 6/sqrt(M*N));
        case 'OLD',       [X,i] = sc_rcctv(Y, struct('T_mode','ball2','G_mode','peng','beta',beta,'sigma',sqrt(mean(sg.^2)),'sig_online',false)); iters = i.iters;
        case 'OLDband',   [X,i] = sc_rcctv(Y, struct('T_mode','ball2','G_mode','peng','beta',beta,'sigma',sg,'sig_online',false)); iters = i.iters;
        case 'SC',        [X,i] = sc_rcctv(Y, struct()); iters = i.iters;
        case 'SCv',       [X,i] = sc_rcctv(Y, struct('stripe_dirs','v')); iters = i.iters;
        case 'SC_scalar', [X,i] = sc_rcctv(Y, struct('scalar_sigma',true)); iters = i.iters;
        case 'SC_noT',    [X,i] = sc_rcctv(Y, struct('T_mode','none')); iters = i.iters;
        case 'SC_noG',    [X,i] = sc_rcctv(Y, struct('G_mode','none','sig_online',false)); iters = i.iters;
        case 'SC_fix',    [X,i] = sc_rcctv(Y, struct('sig_online',false)); iters = i.iters;   % sigma_k from the input only (two-stage)
        case 'E3DTV',     rng(0); X = EnhancedTV(Y, 0.004*sqrt(M*N), [13 13 13]);
        case 'LRTDTV',    rng(0); X = LRTDTV(Y, 1, 10, [round(0.8*M) round(0.8*N) 10]);
        case 'tCTV',      randn('state',0); rand('state',0); X = TCTV_TRPCA(Y, struct('directions',[1 2]));
        otherwise, error('unknown method %s', m);
      end
      catch ME
        fprintf('ERROR %s c%d s%d %s : %s\n', cube, c, s, m, ME.message); continue;
      end
      sec = toc(t0); X = reshape(X,[M,N,p]);
      [ps,ss,sam,erg] = score_full(X, clean);
      fid = fopen(csv,'a'); fprintf(fid,'%s,%d,%d,%s,%.6f,%.6f,%.6f,%.6f,%.2f,%g\n', cube,c,s,m,ps,ss,sam,erg,sec,iters); fclose(fid);
      fprintf('%s c%d s%d %-9s mpsnr %.3f  %.1fs\n', cube,c,s,m,ps,sec);
    end
  end
end
fprintf('SIM DONE %s %s\n', cube, tag);
