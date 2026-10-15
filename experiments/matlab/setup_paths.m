function R = setup_paths()
% Adds SC-RCCTV, official RCCTV (no PROPACK; lansvd stand-in is norm-based), baselines. Returns repo root.
here = fileparts(mfilename('fullpath'));
R = fileparts(fileparts(here));                       % grsl_paper_experiments_2026-10-12
addpath(here);
addpath(fullfile(R,'code','official_rcctv'));
addpath(fullfile(R,'code','official_rcctv','TV_operator'));
addpath(genpath(fullfile(R,'code','official_tctv')));
addpath(fullfile(R,'code','baselines','LRTDTV'));
addpath(fullfile(R,'code','baselines','LRTDTV','prox_operators'));
addpath(fullfile(R,'code','baselines','LRTDTV','mylib'));
e3 = fullfile(fileparts(here),'external','Enhanced-3DTV','Denoising','Enhanced3DTV');
if exist(e3,'dir'), addpath(e3); end
end
