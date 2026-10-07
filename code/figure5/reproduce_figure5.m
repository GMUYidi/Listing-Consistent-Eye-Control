function reproduce_figure5(python_executable)
% Reproduce Figure 5 and its numerical checks from this directory.
% Requires MATLAB and Python with pypdf. Default command: reproduce_figure5
% Optional Python path: reproduce_figure5('C:\path\to\python.exe')
if nargin < 1 || isempty(python_executable)
    python_executable = 'python';
end
code_dir = fileparts(mfilename('fullpath'));
original_path = path;
cleanup = onCleanup(@() path(original_path)); %#ok<NASGU>
addpath(code_dir, '-begin');

verify_dynamic_model_integration;
render_figure5;
pad_script = fullfile(code_dir, 'pad_figure5_pdf.py');
command = sprintf('"%s" "%s"', char(python_executable), pad_script);
[status, output] = system(command);
if status ~= 0
    error('Figure5:PDFMargin', ...
        ['The numerical run and MATLAB export completed, but PDF margin ' ...
         'adjustment failed. Use a Python environment with pypdf installed.\n%s'], output);
end
fprintf('%s\n', output);
fprintf('Figure 5 and numerical results are in: %s\n', fullfile(code_dir, '..', '..', 'outputs', 'figure5'));
end
