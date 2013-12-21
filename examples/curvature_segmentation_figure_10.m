% This example file corresponds to Figure 10 in:
%
% Partial Enumeration and Curvature Regularization
% Carl Olsson, Johannes Ulén, Yuri Boykov and Vladimir Kolmogorov
% International Conference on Computer Vision 2013
%
% In this version the unary costs are correctly weighted leading
% results sligthly differing from the ones presented in the paper.
%
% N.B. early feasible solutions are usually pretty decent.
% set max_relgap = 0.5, or similar, to abort early. 

clear all;
close all;
addpath('../');

% Creating the patches can be done in paralllell.
if matlabpool('size') == 0
	matlabpool('open');
end

im = double(imread('cameraman.tif')) / 255;

unary(:,:,1) = (0.5-im).^2;
unary(:,:,2) = im.^2;

% Regularization strength
lambda = 1;

% Setup object
C = Curvature_Segmentation(unary);

% Settings
C.verbose = false; % 
C.max_iter = 1000; % Stop when maximum number of iterations are reached.
C.max_relgap = 1e-10; % Stops when (e-lb)/e < max_relgap.
C.patch_size = 2;
C.lambda = 1;

%% 2x2
figure(); tic
[L,e,lb] = C.solve();
C.display(); xlabel(sprintf('Time: %g', toc));

%% 3x3
C.patch_size = 3;
figure(); tic;
[L,e,lb] = C.solve();
C.display(); xlabel(sprintf('Time: %g', toc));

%% 5x5
C.max_relgap = 1e-5;
C.patch_size = 5;
figure(); tic;
[L,e,lb] = C.solve();
C.display(); xlabel(sprintf('Time: %g', toc));
