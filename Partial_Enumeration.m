% Johannes Ulén and Carl Olsson 2013
%
classdef Partial_Enumeration < handle
	
	properties
		max_iter = 10000;
		max_relgap = 1e-10;
		unary_patches;
		verbose = false;
	end
	
	properties (SetAccess = protected)
		patch_size = 2;
		labels;
		problem_size;
	end
	
	methods
		
		% Unary (required): width x height x 2.  Where unary(x,y,i) is cost to set label
		% (x,y) to i.
		% --
		% Patch_size (optional, default: 2): Size of one side of the patch. E.g 3 gives 3x3 patches
		% --
		% Labels (optional, default: all): Which of the 2^(patch_size^2) should be
		% used? Input should be a vector, e.g. 0:511 for 3x3 patches.
		
		function self = Partial_Enumeration(unary, patch_size, labels)
			
			self.problem_size = [size(unary,1) size(unary,2)];
			
			if nargin > 1
				self.patch_size = patch_size;
			end
			
			if nargin > 2
				self.labels = labels;
			else
				self.labels = 0:2^(self.patch_size^2)-1;
			end

			if strcmp(self.labels,'all')
				self.labels = 0:2^(self.patch_size^2)-1; 
			end
			
			assert(numel(self.labels) <= 2^(self.patch_size^2));
			
			% All zero
			if nnz(unary) == 0
				w = (self.problem_size(1) - self.patch_size +1);
				h = (self.problem_size(2) - self.patch_size +1);
				self.unary_patches = zeros(numel(self.labels), w*h);
			else
				
				% Setup patch costs
				self.unary_patches = setup_patch_cost(self, unary);
			end
		end
		
		function patch_cost = setup_patch_cost(self, unary)
			ps = self.patch_size;
			patch_variables = self.patch_size*self.patch_size;
			patch_cost = zeros(length(self.labels),(self.problem_size(1)-ps+1)*(self.problem_size(2)-ps+1));
			
			% Each variable is not covered by the same number of patches.
			% We re-weigh the unary term to make up for this.
			weights = zeros(self.problem_size);
			for row = 1:ps
				for col = 1:ps
					weights([0:end-ps]+row,[0:end-ps]+col) =  weights([0:end-ps]+row,[0:end-ps]+col)+1;
				end
            end
			
            % This weighting was used in the paper:
            % weights = patch_variables;  
            
			% Re-weighting
			unary(:,:,1) = unary(:,:,1)./weights;
			unary(:,:,2) = unary(:,:,2)./weights;

			parfor i = 1:length(self.labels);
				binlabels = dec2bin(self.labels(i),patch_variables);
				binlabels = binlabels(end:-1:1);
				patch = reshape(binlabels,[ps ps])';
				datatermi = zeros(self.problem_size-[ps-1,ps-1]);
				for row = 1:ps
					for col = 1:ps
						datatermi = datatermi + ...
							(patch(row,col)=='1').*unary([0:end-ps]+row,[0:end-ps]+col,2) + ...
							(patch(row,col)=='0').*unary([0:end-ps]+row,[0:end-ps]+col,1);
					end
				end
				
				patch_cost(i,:) = datatermi(:);
			end
		end
		
		function self = add_patch_cost(self, patch_cost)
			% Add a certain penalty to every n x n patch.
			if (numel(patch_cost) ~= numel(self.labels));
				error('You must define a cost for each configuration in the patch.');
			end
			
			self.unary_patches =  self.unary_patches + ...
				repmat(patch_cost(:),[1 size(self.unary_patches,2)]);
		end
		
		function [L,e,lb] = solve(self)
			
			% Output
			% L: Labelling
			% e: Energy
			% lb: lower bound
			settings.max_iter = self.max_iter;
			settings.max_relgap = self.max_relgap;
			settings.verbose = self.verbose;
			
			[conn, conn_type] = patch_consistency(self);

			[L, e,lb] = trws_partial_enumeration( ...
				self.unary_patches, ...
				conn, ...
				conn_type, ...
				self.patch_size, ...
				self.labels, ...
				settings);

			L = patch_to_pixel(self, L);

		end
		function [connectivity, pairwise] = patch_consistency(self)
			
			nodenr = zeros(self.problem_size - self.patch_size + 1);
			nodenr(:) = 1:length(nodenr(:));
			
			%Neighborhood
			start = nodenr(1:end-1,:);
			finish = nodenr(2:end,:);
			neighborhood.vertedge = [start(:) finish(:)];
			
			start = nodenr(:,1:end-1);
			finish = nodenr(:,2:end);
			neighborhood.horedge = [start(:) finish(:)];
			
			connectivity = [neighborhood.vertedge(:,1)' neighborhood.horedge(:,1)'; ...
				neighborhood.vertedge(:,2)' neighborhood.horedge(:,2)'];
			
			pairwise = [ones(size(neighborhood.horedge,1),1);
				2*ones(size(neighborhood.vertedge,1),1)]';
		end
		
		
		function pix_assignments = patch_to_pixel(self, L)
			assert(min(L(:)) > 0);
		
			bin_ind = dec2bin(self.labels(L),self.patch_size^2);
			bin_ind = bin_ind(:,end:-1:1);
		
			pix_assignments = zeros(self.problem_size(1:2));
	
			i = 0;
			for col = 1:self.patch_size
				for row = 1:self.patch_size
					i=i+1;
					pix_assignments(row-1+[1:end-self.patch_size+1],col-1+[1:end-self.patch_size+1]) ...
						= reshape(bin_ind(:,i)=='1', self.problem_size -self.patch_size + 1);
				end
			end
		end

		% Manually set unary patch cost
		function self = set.unary_patches(self,value)
			
			% Patch size
			w = (self.problem_size(1) - self.patch_size +1);
			h = (self.problem_size(2) - self.patch_size +1);
			unary_patch_size = [numel(self.labels) w*h];
			
			if ~( all(unary_patch_size == size(value)) )
				error('Input must be of same size OBJ.unary_patches');
			end
			
			self.unary_patches = value;
		end
	end
end