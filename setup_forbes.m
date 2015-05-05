% Add ForBES directory to MATLAB's path
forbes_path = fileparts(mfilename('fullpath'));
display(['Adding ForBES directory to MATLAB path: ', forbes_path]);
addpath(forbes_path);
savepath;

% Compile necessary C source files
private_path = fullfile(forbes_path, 'private');
tests_path = fullfile(forbes_path, 'tests');
LBFGS_path = fullfile(forbes_path, 'private', 'LBFGS.c');
Riccati_path = fullfile(forbes_path, 'tests', 'RiccatiSolve.c');
error_msg = 'The C compiler could not succesfully compile ';
if mex('-outdir', private_path, LBFGS_path), error([error_msg, LBFGS_path]); end
if mex('-outdir', tests_path, Riccati_path), error([error_msg, Riccati_path]); end
display('ForBES was succesfully configured and installed');
