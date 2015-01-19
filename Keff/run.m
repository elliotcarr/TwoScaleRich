% ------------------------------------------------------------------------------
% run.m
%
% Written by Dr Elliot Carr (2013-2014)
% Ecole Centrale Paris and Queensland University of Technology
% 
% This code is part of TwoScalRich.
%
% This MATLAB scripts computes the effective conductivity over the micro-cell
% and stores the appropriate results to a .txt file.
%
% ------------------------------------------------------------------------------

commandwindow
clc
clear all
close all

%-------------------------------------------------------------------------------
% User changeable parameters 

% Image-based configuration
visualise      = true;
Keff_file      = 'Keff.txt';
mesh_name      = 'mesh_Keff.mat'; % see micro_mesh_unstructured.m

% Solve periodic cell problems at a geometric sequence of discrete h values
n = 100; % number of discrete values
h_range = zeros(n,1);
h_range(1) = -2.01e1;
h_range(n) = -1e-2;
r = exp((log(-h_range(n))-log(-h_range(1)))/(n-1));
for i = 2:n
    h_range(i) = h_range(i-1)*r;
end

% Hydaulic properties of soils
% Sub-domain A
soilA.Ksat   = 0.044 / 3600;
soilA.thetar = 0.058;
soilA.thetas = 0.41;
soilA.alpha  = 7.3;
soilA.n      = 1.89;
soilA.m      = 1 - 1/soilA.n;

% Sub-domain B (inclusions)
soilB.Ksat   = 0.044 / 3600 * 1e-3;
soilB.thetar = 0.058;
soilB.thetas = 0.41;
soilB.alpha  = 7.3;
soilB.n      = 1.89;
soilB.m      = 1 - 1/soilB.n;
%-------------------------------------------------------------------------------

% Compile mex functions
mex Ffunc1.c
mex Ffunc2.c
mex compute_Keff.c

user_data.soilA = soilA;
user_data.soilB = soilB;

% Mesh micro-cell geometry for periodic cell problem and generate properties 
% (CV areas, edges lengths, etc.)
user_data.micro_mesh = micro_mesh(mesh_name, visualise);

% Solve periodic cell problems
effective_conductivity(Keff_file, h_range, user_data, visualise);