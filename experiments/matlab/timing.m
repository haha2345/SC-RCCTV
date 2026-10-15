% Clean wall-clock benchmark (run on an otherwise idle machine, default MATLAB threads, one run per cell).
% Resumable: skips (ds,meth) rows already in results/timing.csv; errors are printed and skipped.
R = setup_paths(); [base, inb] = res_dirs(R); f = fullfile(base,'timing.csv');
done = '';
if exist(f,'file'), done = fileread(f); else, fid = fopen(f,'w'); fprintf(fid,'ds,meth,sec\n'); fclose(fid); end
jobs = {'pavia120_c8','belt100_c8','indian220'};
meths = {'RCCTV','SC','E3DTV','LRTDTV'};          % t-CTV (dev-tuned lambda) and the other new baselines: timing_new.m
for ji = 1:numel(jobs)
  ds = jobs{ji};
  if strncmp(ds,'indian',6)
    S = load(fullfile(R,'data','real','IndianPines','Indian_pines.mat')); Y = double(S.indian_pines);
    for b = 1:size(Y,3), a = Y(:,:,b); mn = min(a(:)); mx = max(a(:)); Y(:,:,b) = (a-mn)/max(mx-mn,eps); end
  else
    cube = ds(1:end-3); I = load(fullfile(inb,'inputs',[cube '_c8_s0.mat'])); Y = I.Nhsi;
  end
  [M,N,p] = size(Y);
  for mi = 1:numel(meths)
    m = meths{mi};
    if contains(done, [ds ',' m ',']), continue; end
    try
      t0 = tic;
      switch m
        case 'RCCTV',  X = RCCTV_RPCA(Y,8,6/sqrt(M*N));
        case 'SC',     X = sc_rcctv(Y,struct());
        case 'E3DTV',  rng(0); X = EnhancedTV(Y,0.004*sqrt(M*N),[13 13 13]);
        case 'LRTDTV', rng(0); X = LRTDTV(Y,1,10,[round(.8*M) round(.8*N) 10]);
        case 'tCTV',   randn('state',0); rand('state',0); X = TCTV_TRPCA(Y,struct('directions',[1 2]));
      end
      sec = toc(t0); fid = fopen(f,'a'); fprintf(fid,'%s,%s,%.2f\n',ds,m,sec); fclose(fid); fprintf('%s %s %.1fs\n',ds,m,sec);
    catch ME
      fprintf('ERROR %s %s : %s\n', ds, m, ME.message);
    end
  end
end
disp('TIMING DONE');
