% Environment report and prerequisite check for the full re-run (run_all_main.ps1, phase check).
% Writes <RESULTS_DIR>/env_main.txt and stops with an error if a toolbox or a third-party code is missing.
R = setup_paths_new(); warning('off','MATLAB:mpath:nameNonexistentOrNotADirectory');
base = res_dirs(R);
f = fopen(fullfile(base,'env_main.txt'),'w');
fprintf(f,'date %s\nMATLAB %s\ncomputer %s\nmaxNumCompThreads %d\n', char(datetime('now')), version, computer, maxNumCompThreads);
fprintf(f,'BLAS %s\nLAPACK %s\n', version('-blas'), version('-lapack'));
v = ver; for i = 1:numel(v), fprintf(f,'toolbox %s %s\n', v(i).Name, v(i).Version); end
ok = true;
chk = {'psf2otf','Image Processing Toolbox'; 'fitcecoc','Statistics and Machine Learning Toolbox'; ...
       'pagesvd','MATLAB R2021b or newer (pagesvd)'};
for i = 1:size(chk,1)
  if exist(chk{i,1}) == 0, fprintf(f,'MISSING %s\n', chk{i,2}); ok = false; end %#ok<EXIST>
end
if ~license('test','Signal_Toolbox'), fprintf(f,'MISSING Signal Processing Toolbox\n'); ok = false; end
code = {'RCCTV_RPCA','official RCCTV'; 'EnhancedTV','E-3DTV'; 'LRTDTV','LRTDTV'; 'TCTV_TRPCA','t-CTV'; ...
        'FastHyMix','FastHyMix'; 'L1HyMixDe','L1HyMixDe'; 'bm3d_thr','BM3D mex'; 'Prox_S3TTV','S3TTV'};
for i = 1:size(code,1)
  if exist(code{i,1}) == 0, fprintf(f,'MISSING code %s (%s)\n', code{i,2}, code{i,1}); ok = false; end %#ok<EXIST>
end
st = 'OK'; if ~ok, st = 'FAILED'; end
fprintf(f,'STATUS %s\n', st); fclose(f);
type(fullfile(base,'env_main.txt'));
if ~ok, error('env_check:missing','environment check failed, see env_main.txt'); end
