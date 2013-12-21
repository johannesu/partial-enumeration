classdef Binary_Deconvolution < handle
    properties (SetAccess = protected)
        L;
        e;
        lb;
        im;
        fullsz;
        
        % Pseudo binary functions
        pbf;
			
				% Partial enumeration
				pe;

        % Solver used to create solution
        solver;
		end
		
		properties
        % Settings
        max_relgap = 1e-14;
        max_iter = 300;
				verbose = false;
		end
    
    methods
        function self = Binary_Deconvolution(im, fullsz)
            
            % Add appropriate path
            my_name = mfilename('fullpath');
            my_path = fileparts(my_name);
            addpath([my_path filesep 'include']);
            
            
            
            self.im = im;
            self.fullsz = fullsz;
            
            % Setup cost functions
            setup_pbf_cost(self)
            setup_partial_enumeration_cost(self);
            
            self.L = zeros(self.fullsz);
            self.e = nan;
            self.lb = nan;
        end
        
        function [L,e,lb] = solve(self, solver, settings)
            
            % Default partial enumeration
            if nargin == 1
                solver = 'pe';
            end
            
            if nargin == 2
                settings = [];
            end
            
            if strcmp(solver,'RD')
                [L, e, lb, unlabelled] = rd(self.pbf.unary(1,:)', ...
                    self.pbf.unary(2,:)', ...
                    self.pbf.pairwise(1,:), ...
                    self.pbf.pairwise(2,:), ...
                    self.pbf.pairwise(3,:), ...
                    self.pbf.pairwise(4,:), ...
                    self.pbf.connectivity);
                
                % Unlaballed set to 0.5
                L(L < 0) = 0.5;
                
                % Store
                self.solver = 'Roof duality (RD)';
                
            elseif (strcmp(solver,'trws'))
               [L,e,lb] = trws_general( ...
                                self.pbf.unary, ...
                                self.pbf.connectivity, ...
                                self.pbf.pairwise, ...
                                settings);
  
                self.solver = 'TRW-S';
                
            elseif (strcmp(solver,'Partial Enumeration'))
							
							% Alternative way to to use partial enumeration set the patch
							% costs directly. This is needed if you want irregular patch
							% costs
							
							% Zero unary costs
							unary = sparse(self.fullsz(1),self.fullsz(2));
							patch_size = 3;
							PE = Partial_Enumeration(unary, patch_size);
							
							PE.verbose = self.verbose;
							PE.max_iter = self.max_iter;
							PE.max_relgap = self.max_relgap;
     
							% Update all patch costs.
							PE.unary_patches = self.pe.unary;	
							[L,e,lb] = PE.solve();
							
							self.solver = 'Partial enumeration optimized by TRW-S';  
								                    
            else
                error('Unkown solver');
            end
            
            self.L = reshape(L,self.fullsz); 
            self.e = e; 
            self.lb = lb;
        end
        
        function display(self)
           imagesc(self.L);
           colormap gray; axis equal;
           gap = abs(self.e - self.lb);
           
           title(sprintf('Solver: %s', self.solver));
           xlabel(sprintf('(Lower bound, Energy): (%g,%g) Duality gap: %g', ...
                           self.e,self.lb, gap));
        end
        % Second order pbf.
        function setup_pbf_cost(self)
            nodenr = zeros(self.fullsz);
            nodenr(:) = 1:length(nodenr(:));
            
            %Unary term
            datacost = zeros(self.fullsz);
            datacost(1:end-2,1:end-2) = datacost(1:end-2,1:end-2) + 2*self.im/9;
            datacost(1:end-2,2:end-1) = datacost(1:end-2,2:end-1) + 2*self.im/9;
            datacost(1:end-2,3:end) = datacost(1:end-2,3:end) + 2*self.im/9;
            datacost(2:end-1,1:end-2) = datacost(2:end-1,1:end-2) + 2*self.im/9;
            datacost(2:end-1,2:end-1) = datacost(2:end-1,2:end-1) + 2*self.im/9;
            datacost(2:end-1,3:end) = datacost(2:end-1,3:end) + 2*self.im/9;
            datacost(3:end,1:end-2) = datacost(3:end,1:end-2) + 2*self.im/9;
            datacost(3:end,2:end-1) = datacost(3:end,2:end-1) + 2*self.im/9;
            datacost(3:end,3:end) = datacost(3:end,3:end) + 2*self.im/9;
            datacost = -datacost(:);
            
            %Pairwise term
            edges = [];
            weights = [];
            for i = 1:9
                for j = 1:9
                    
                    if i==j
                        row = floor(((i-1)/3))+1;
                        col = mod(i-1,3)+1;
                        
                        start = nodenr([1:end-2]+row-1,[1:end-2]+col-1);
                        datacost(start(:)) = datacost(start(:))+1/9^2;
                    else
                        row = floor(((i-1)/3))+1;
                        col = mod(i-1,3)+1;
                        
                        start = nodenr([1:end-2]+row-1,[1:end-2]+col-1);
                        
                        row = floor(((j-1)/3))+1;
                        col = mod(j-1,3)+1;
                        
                        finish = nodenr([1:end-2]+row-1,[1:end-2]+col-1);
                        
                        newedges = sort([start(:)';finish(:)'],1);
                        
                        [mem,loc] = ismember(newedges',edges','rows','legacy');
                        weights(loc(mem)) = weights(loc(mem)) + 1/9^2;
                        
                        edges = [edges newedges(:,~mem)];
                        weights = [weights 1/9^2*ones(1,sum(~mem))];
                    end
                end
            end
            
            U0 = zeros(size(datacost));
            U1 = datacost;
            
            % Adding constant cost
            U0(1) = sum(self.im(:).^2);
            U1(1) = sum(self.im(:).^2);
            
            E00 = zeros(size(weights));
            E01 = zeros(size(weights));
            E10 = zeros(size(weights));
            E11 = weights;
            
            self.pbf.unary = [U0 U1]';
            self.pbf.pairwise = [E00; E01; E10; E11];
            self.pbf.connectivity = edges;
        end
        
        % 3x3 patches
        function setup_partial_enumeration_cost(self)
            patch_cost = zeros(2^9, numel(self.im));
            for i = 1:2^9
                ibin = dec2bin(i-1);
                ibin = ibin(end:-1:1);
                imat = (ibin=='1');
                
                nrfgpix = sum(imat);
                costi = (self.im - nrfgpix/9).^2;
                patch_cost(i,:) = costi(:)';
            end

            self.pe.unary = patch_cost;
        end
    end
end