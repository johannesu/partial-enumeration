function [solution, energy, lower_bound, iterations] =  trws_partial_enumeration(unary, connectivity, pairwise, patch_size, labels, options)

	type = 'partial_enumeration';

	[solution, energy, lower_bound, iterations] = trws(type, unary, connectivity, pairwise, patch_size, labels, options);
      
