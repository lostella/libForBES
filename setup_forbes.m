% Add ForBES directory to MATLAB's path
forbes_path = fileparts(mfilename('fullpath'));
addpath(forbes_path);
savepath;

% Compile necessary C source files
private_path = fullfile(forbes_path, 'private');
tests_path = fullfile(forbes_path, 'tests');
LBFGS_path = fullfile(forbes_path, 'private', 'LBFGS.c');
Riccati_path = fullfile(forbes_path, 'tests', 'RiccatiSolve.c');
mex('-outdir', private_path, LBFGS_path);
mex('-outdir', tests_path, Riccati_path);