% Compile necessary C source files
mex -outdir private/ private/LBFGS.c
mex -outdir tests/ tests/RiccatiSolve.c

% Add ForBES directory to MATLAB's path
this_path = mfilename('fullpath');
split_path = strsplit(this_path, '/');
forbes_path = split_path(:,1:end-1);
addpath(forbes_path);
savepath;