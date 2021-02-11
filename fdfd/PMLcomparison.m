%% Compares the UPML in Rumpf and the SCPML in Shin
%   Uses new function for FDFD and param structure
%% Mesh specifications
% Number of cells to split vacuum wavelengths into
param2D.lx = 40; 
param2D.ly = 40;
% Number of cells for PML
param2D.Lx = 30; 
% Enforce distance of PML from TF region
param2D.bufsize = 210; 
% Total number of cells for x and y axis
param2D.Nx = 800 + param2D.Lx + param2D.bufsize;
param2D.Ny = 201;
% Set the index where reflectance and transmittance are calculated
%param2D.Rpind = param2D.Lx + 5;
%param2D.Tpind = param2D.Nx - param2D.Lx - 5;
param2D.Rpind = param2D.lx;
param2D.Tpind = param2D.Nx - param2D.lx;

%% PML parameters
% Set the polynomial order. (1 <= p < 5)
param2D.sx_m = 3;
% Desired reflectivity coefficient
param2D.sx_R = 1e-9;
% Set \sigma_max
param2D.sx_sigmax = 0.004;
% Set s_max. (0 <= s_max <= 5)
param2D.sx_smax = 3;

%% Space parameters
% Set the boundaries 
% -- Currently only for region split into two halves (ADD FUNCTIONALITY LATER)
% param2D.xbounds = [];
% Set the relative permittivities and permeabilities in each region
%   Only for second region in current implementation.
param2D.epsr = 1;
param2D.mur = 1;
param2D.Floquet = 1;

%%
paramsR = param2D
paramsS = param2D;

%% Set the omeg sampling and and angle of incidence
Nf = 2;
Na = 2;
%omegs = linspace(430e12, 750e12, Nf);
%thetai = linspace(0.0, 0.35*pi, Na)';
omegs = [430e12; 550e12];
thetai = [0.0; 0.1];

R1 = zeros(Na, Nf);
T1 = zeros(Na, Nf);
R2 = zeros(Na, Nf);
T2 = zeros(Na, Nf);
%% Test
for m = 1:Na
    for n = 1:Nf
        pathhead = sprintf("PMLtest/UPML_om%02dth%02d_", n, m);
        [R1(m,n),T1(m,n)] = fdfd2D(omegs(n), thetai(m), paramsR, pathhead)
        pathhead = sprintf("PMLtest/SCPML_om%02dth%02d_", n, m);
        [R2(m,n),T2(m,n)] = fdfd2D_SCPML(omegs(n), thetai(m), paramsR, pathhead)
	sprintf("Both BC with omeg = %f and thetai = %f complete. (%d,%d)", omegs(n), thetai(m), m, n)
    end
end

%% Find accuracy of reflection, transmittance, and energy cons
% Angle of refraction from Snell's law
n1 = 1;
n2 = sqrt(param2D.epsr);
thetat = asin(n1 / n2 * sin(thetai));
% Compute Fresnel relations
R_Fres = ((n1*cos(thetai)-n2*cos(thetat))./(n1*cos(thetai)+n2*cos(thetat))).^2
T_Fres = 1 - R_Fres; % Assumes energy conserved

% Compute mean and std dev for calculated R and T
R1_mean = mean(R1,2)
R1_std = std(R1, 0, 2)
T1_mean = mean(T1,2)
T1_std = std(T1, 0, 2)
Econs1 = R1_mean + T1_mean

R2_mean = mean(R2,2)
R2_std = std(R2, 0, 2)
T2_mean = mean(T2,2)
T2_std = std(T2, 0, 2)
Econs2 = R2_mean + T2_mean

%% Compute the relative errors
R1_relerr = abs(R1_mean - R_Fres) ./ R_Fres 
T1_relerr = abs(T1_mean - T_Fres) ./ T_Fres 

R2_relerr = abs(R2_mean - R_Fres) ./ R_Fres 
T2_relerr = abs(T2_mean - T_Fres) ./ T_Fres 

Ec1_relerr = abs(Econs1 - 1)
Ec2_relerr = abs(Econs2 - 1)

%% Save the relevant data
finaldat = table(thetai, R1_mean, R2_mean, R_Fres, Ec1_relerr, Ec2_relerr, R1_relerr, R2_relerr, T1_relerr, T2_relerr)
writetable(finaldat, "bigtest/RefDatPMLComp.txt");
