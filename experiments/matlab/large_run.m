% Large-scale run (accuracy + wall-clock on an otherwise idle machine; one method at a time).
% env LARGE_JOBS (comma list of paviaU_c8, paviaU_c10, urban307), LARGE_METHS (comma list), TCTV_MULT.
% Appends results/large.csv: job,meth,mpsnr,mssim,sam,ergas,stripe_top20,sec  (urban: no reference -> NaN scores)
R = setup_paths_new();
warning('off','MATLAB:mpath:nameNonexistentOrNotADirectory');
[base, inb] = res_dirs(R);
csv = fullfile(base,'large.csv');
jobs = strsplit(getenv('LARGE_JOBS'), ','); meths = strsplit(getenv('LARGE_METHS'), ',');
tm = str2double(getenv('TCTV_MULT')); if isnan(tm), tm = 1; end
done = {};
if exist(csv,'file'), done = strsplit(fileread(csv), newline); else, fid = fopen(csv,'w'); fprintf(fid,'job,meth,mpsnr,mssim,sam,ergas,stripe_top20,sec,M,N,p\n'); fclose(fid); end
for ji = 1:numel(jobs)
  job = jobs{ji};
  if strcmp(job,'urban307')
    S = load(fullfile(R,'data','real','Urban','Urban_F210.mat')); Y = double(S.Urban_F210); clean = [];
    for b = 1:size(Y,3), a = Y(:,:,b); mn = min(a(:)); mx = max(a(:)); Y(:,:,b) = (a-mn)/max(mx-mn,eps); end
    o = struct('tctv_mult', tm);
  else
    I = load(fullfile(inb,'inputs_large',[job '_s1.mat'])); Y = I.Nhsi; clean = I.clean;
    c = sscanf(job,'paviaU_c%d'); SP = [0 .05 .05 0 0 0 .05 .05 0 .05]; ST = [0 0 0 1 1 1 1 1 0 1]*0.05;
    o = struct('sigma', sqrt(mean(I.sigma_true(:).^2)), 'sp_rate', SP(c), 'st_rate', ST(c), 'tctv_mult', tm);
  end
  [M,N,p] = size(Y);
  SIy = stripe_index(Y); [~,ordI] = sort(SIy,'descend'); top = ordI(1:20);
  for mi = 1:numel(meths)
    m = meths{mi}; key = sprintf('%s,%s,', job, m);
    if any(strncmp(done, key, numel(key))), continue; end
    t0 = tic;
    try
      switch m
        case 'RCCTV',  X = RCCTV_RPCA(Y, 8, 6/sqrt(M*N));
        case 'SC',     X = sc_rcctv(Y, struct());
        case 'E3DTV',  rng(0); X = EnhancedTV(Y, 0.004*sqrt(M*N), [13 13 13]);
        case 'LRTDTV', rng(0); X = LRTDTV(Y, 1, 10, [round(0.8*M) round(0.8*N) 10]);
        otherwise,     X = run_new_method(m, Y, o);
      end
    catch ME
      fprintf('ERROR %s %s : %s\n', job, m, ME.message);
      fid = fopen(csv,'a'); fprintf(fid,'%s,%s,NaN,NaN,NaN,NaN,NaN,NaN,%d,%d,%d\n', job, m, M, N, p); fclose(fid); continue;
    end
    sec = toc(t0); X = reshape(X,[M,N,p]);
    si = stripe_index(X); st = mean(si(top));
    if isempty(clean), ps = NaN; ss = NaN; sam = NaN; erg = NaN; else, [ps,ss,sam,erg] = score_full(X, clean); end
    fid = fopen(csv,'a'); fprintf(fid,'%s,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.2f,%d,%d,%d\n', job, m, ps, ss, sam, erg, st, sec, M, N, p); fclose(fid);
    fprintf('LARGE %s %-9s mpsnr %.3f stripe %.4f %.1fs\n', job, m, ps, st, sec);
  end
end
fprintf('LARGE DONE\n');

function si = stripe_index(Y)
[~,~,p] = size(Y); c = squeeze(mean(Y,1)); r = squeeze(mean(Y,2));
si = zeros(1,p);
for b = 1:p, v = c(:,b); w = r(:,b); si(b) = mean(abs(v - movmedian(v,7))) + mean(abs(w - movmedian(w,7))); end
end
