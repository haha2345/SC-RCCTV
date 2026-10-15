function R = setup_paths_new()
% setup_paths() + the baselines added in the revision (FastHyMix, L1HyMixDe, S3TTV helpers).
% BM3D mex files are byte-identical to the official TUNI package (SHA-256 checked); only one copy is added.
R = setup_paths();
X = fullfile(fileparts(fileparts(mfilename('fullpath'))),'external');
addpath(fullfile(X,'FastHyMix','scripts'));
addpath(fullfile(X,'FastHyMix','BM3D'));
addpath(fullfile(X,'L1HyMixDe','L1HyMixDe'));
addpath(fullfile(X,'S3TTV','sub_functions','prox'));
addpath(fullfile(X,'S3TTV','sub_functions','operator'));
end
