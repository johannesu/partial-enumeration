% Wrapper function for TRW-S solver
function [solution, energy, lower_bound, iterations] =  ...
			trws(type, unary, connectivity, pairwise, patch_size, labels, options )
      
labels = int32(labels);
patch_size = int32(patch_size);
assert(patch_size > 0);

if (nargin == 5)
    options = [];
else
    if isfield(options,'max_iter')
       options.max_iter = int32(options.max_iter); 
    end
    
    if isfield(options,'verbose')
        options.verbose = logical(options.verbose);
    end
end
assert(min(connectivity(:) > 0));
assert( max(connectivity(:) <= numel(unary)))

% Compile if need be
my_name = mfilename('fullpath');
my_path = fileparts(my_name);

cpp_file = ['trws_mex.cpp'];
out_file = ['trws_mex'];

extra_arguments = {['-I"' my_path '"']};
sources = {['trw-s' filesep 'MRFEnergy.cpp'], ...
		['trw-s' filesep 'minimize.cpp'], ...
		['trw-s' filesep 'ordering.cpp'], ...
		['trw-s' filesep 'treeProbabilities.cpp']};

compile(cpp_file, out_file, sources, extra_arguments);

% Solve
[solution, energy, lower_bound, iterations] = ...
 trws_mex(type, patch_size, unary, uint32(connectivity-1), pairwise, labels, options);
