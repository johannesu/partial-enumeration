function [solution, energy, lower_bound, iterations] =  ...
		trws_general(unary, connectivity, pairwise, options)

	type = 'general';
	patch_size = 1;
	labels = 0;

	[solution, energy, lower_bound, iterations] =  ...
	trws(type, unary, connectivity, pairwise, patch_size, labels,  options);
      