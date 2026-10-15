function [out, inp] = res_dirs(R)
% Output and input folders of all experiment scripts.
% env RESULTS_DIR redirects every output (CSV, real/*.mat, figs): a name relative to the experiment folder
% (e.g. results_main for the full re-run on the main machine) or an absolute path; default 'results'.
% The noisy inputs (inputs/, inputs_dev/, inputs_large/) are always read from <experiment folder>/results/.
nm = fileparts(fileparts(mfilename('fullpath')));   % experiment folder (parent of matlab/), whatever its name
d = getenv('RESULTS_DIR');
if isempty(d), d = 'results'; end
if numel(d) > 1 && (d(2) == ':' || d(1) == filesep || d(1) == '/'), out = d; else, out = fullfile(nm, d); end
inp = fullfile(nm, 'results');
if ~exist(out, 'dir'), mkdir(out); end
end
