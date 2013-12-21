% This script performs the experiment of Figure 12 in [1].
%
% [1] Partial Enumeration and Curvature Regularization
% Carl Olsson, Johannes Ulén, Yuri Boykov and Vladimir Kolmogorov
% International Conference on Computer Vision 2013
%
clear all; close all;

% Load image
addpath('../')
im = double(rgb2gray(imread('../data/iccv.png')))/255;

% Add noise
im_noise = conv2(single(im),ones(3,3)./9,'valid')+0.05*randn(size(im)-2);

figure(), imagesc(im_noise); colormap gray;
title('Noisy image'); axis equal;

% Create problem instance
B = Binary_Deconvolution(im_noise, size(im));

%% Solve using Roof duality (12 b)
B.solve('RD');
figure(); B.display();

%% Solve using TRW-S (12 d)
B.verbose = false;
B.solve('trws');
figure(); B.display();

%% Solve by reformulating the problem via Partial Enumeration
%  and then using TRW-S on the new problem instance (12 e)
B.solve('Partial Enumeration');
figure(); B.display();