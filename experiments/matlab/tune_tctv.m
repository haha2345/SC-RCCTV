% Dev-set tuning of the t-CTV lambda multiplier (dev cubes only; never the evaluation cubes).
% env TUNE_CUBES (comma list), TUNE_MULTS (comma list).  Appends results/tune_tctv.csv (resumable).
R = setup_paths_new();
[base, inb] = res_dirs(R);
csv = fullfile(base,'tune_tctv.csv');
tag = getenv('TUNE_TAG'); if ~isempty(tag), csv = fullfile(base,['tune_tctv_' tag '.csv']); end   % one file per parallel job
cubes = strsplit(getenv('TUNE_CUBES'), ',');
mults = str2double(strsplit(getenv('TUNE_MULTS'), ','));
done = {};
if exist(csv,'file'), done = strsplit(fileread(csv), newline); else, fid = fopen(csv,'w'); fprintf(fid,'cube,case,mult,mpsnr,mssim,sam,ergas,sec,iters\n'); fclose(fid); end
for ci = 1:numel(cubes)
  for c = 1:10
    I = load(fullfile(inb,'inputs_dev',sprintf('%s_c%d_s1.mat',cubes{ci},c)));
    for mu = mults
      key = sprintf('%s,%d,%g,', cubes{ci}, c, mu);
      if any(strncmp(done, key, numel(key))), continue; end
      t0 = tic; [X,info] = run_new_method('tCTVd', I.Nhsi, struct('tctv_mult',mu)); sec = toc(t0);
      [ps,ss,sam,erg] = score_full(X, I.clean);
      fid = fopen(csv,'a'); fprintf(fid,'%s,%d,%g,%.6f,%.6f,%.6f,%.6f,%.2f,%g\n',cubes{ci},c,mu,ps,ss,sam,erg,sec,info.iters); fclose(fid);
      fprintf('TUNE %s c%d mult %g mpsnr %.3f %.0fs\n', cubes{ci}, c, mu, ps, sec);
    end
  end
end
fprintf('TUNE DONE\n');
