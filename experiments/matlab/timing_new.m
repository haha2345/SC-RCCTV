% Idle-machine wall-clock of the revision baselines, same jobs/protocol as timing.m (default threads).
% S3TTV: per-iteration time from 20 iterations; total = per-iteration time x iterations of the
% evaluation run (results/sim_s3*_a.csv), because one full run takes hours on this CPU.
R = setup_paths_new();
warning('off','MATLAB:mpath:nameNonexistentOrNotADirectory');
[base, inb] = res_dirs(R);
f = fullfile(base,'timing_new.csv');
done = '';
if exist(f,'file'), done = fileread(f); else, fid = fopen(f,'w'); fprintf(fid,'ds,meth,sec,note\n'); fclose(fid); end
tm = str2double(getenv('TCTV_MULT')); if isnan(tm), tm = 0.5; end
jobs = {'pavia120_c8','belt100_c8','indian220'};
for ji = 1:numel(jobs)
  ds = jobs{ji};
  if strncmp(ds,'indian',6)
    S = load(fullfile(R,'data','real','IndianPines','Indian_pines.mat')); Y = double(S.indian_pines);
    for b = 1:size(Y,3), a = Y(:,:,b); Y(:,:,b) = (a-min(a(:)))/max(max(a(:))-min(a(:)),eps); end
    o = struct('tctv_mult', tm); meths = {'FastHyMix','L1HyMixDe','tCTVd'};
  else
    cb = strtok(ds,'_'); I = load(fullfile(inb,'inputs',[cb '_c8_s0.mat'])); Y = I.Nhsi;
    o = struct('sigma',0.1,'sp_rate',0.05,'st_rate',0.05,'tctv_mult',tm); meths = {'FastHyMix','L1HyMixDe','tCTVd','S3TTV'};
  end
  for mi = 1:numel(meths)
    m = meths{mi}; note = '';
    if contains(done, [ds ',' m ',']), continue; end
    % S3TTV on the CPU, as in its evaluation runs (sim_s3*.csv). Without use_gpu it picks the GPU whenever one is
    % visible, and where pagesvd does not accept gpuArray (R2022a) it falls back to the per-block prox, ~20x slower
    % (75 s/iter measured on 2026-10-14 vs 3.6 s/iter on the CPU).
    if strcmp(m,'S3TTV'), o.maxiter = 20; o.use_gpu = false; end
    t0 = tic; run_new_method(m, Y, o); sec = toc(t0);
    if strcmp(m,'S3TTV')
      C = readtable(fullfile(base, sprintf('sim_s3%s.csv', upper(cb(1)))));   % written by the S3TTV runs (tags s3P / s3B)
      vn = C.Properties.VariableNames; cc = vn{ismember(vn, {'case','case_','xCase'})};   % 'case' is renamed by readtable
      C = C(C.(cc)==8 & C.seed==0, :); it = C.iters(1);
      note = sprintf('%.2fs/iter x %d iters', sec/20, it); sec = sec/20*it;
    end
    fid = fopen(f,'a'); fprintf(fid,'%s,%s,%.2f,%s\n', ds, m, sec, note); fclose(fid);
    fprintf('TIMING %s %s %.1fs %s\n', ds, m, sec, note);
  end
end
