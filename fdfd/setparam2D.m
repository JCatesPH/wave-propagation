%% Set the default values for the fdfdparam structure
% Author: Jalen Cates
% Date: 2/4/2021
% -------------------------
%% Mesh specifications
% Number of cells to split vacuum wavelengths into
param2D.lx = 50; 
param2D.ly = 50;
% Number of cells for PML
param2D.Lx = 40; 
% Enforce distance of PML from TF region
param2D.bufsize = 80; 
% Total number of cells for x and y axis
param2D.Nx = 500 + param2D.Lx + param2D.bufsize;
param2D.Ny = 501;
% Set the index where reflectance and transmittance are calculated
param2D.Rpind = param2D.Lx + 5;
param2D.Tpind = param2D.Nx - param2D.Lx - 5;

%% PML parameters
% Set the polynomial order. (1 <= p < 5)
param2D.sx_m = 3;
% Desired reflectivity coefficient
param2D.sx_R = 1e-8;
% Set \sigma_max
param2D.sx_sigmax = 1;
% Set s_max. (0 <= s_max <= 5)
param2D.sx_smax = 3;

%% Space parameters
% Set the boundaries 
% -- Currently only for region split into two halves (ADD FUNCTIONALITY LATER)
% param2D.xbounds = [];
% Set the relative permittivities and permeabilities in each region
%   Only for second region in current implementation.
param2D.epsr = 2.25;
param2D.mur = 1;
param2D.Floquet = 1;



