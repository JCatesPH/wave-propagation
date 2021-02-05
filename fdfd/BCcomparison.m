%% Compares the Floquet and basic periodic boundary conditions
%   Uses new function for FDFD and param structure
setparam2D

%% Adjust parameters for case of Floquet and basic periodic BC
paramsF = param2D
paramsP = param2D;
paramsP.Floquet = 0

%% Set the omeg sampling and and angle of incidence
Nf = 5;
Na = 15;
omegs = linspace(430e12, 750e12, Nf);
thetai = linspace(0.0, pi/3, Na)';

RF = zeros(Na, Nf);
TF = zeros(Na, Nf);
RP = zeros(Na, Nf);
TP = zeros(Na, Nf);
%% Test
for m = 1:Na
    for n = 1:Nf
        pathhead = sprintf("testfun/Floq_om%02d%02d_", m, n);
        [RF(m,n),TF(m,n)] = fdfd2D(omegs(n), thetai(m), paramsF, pathhead);
        pathhead = sprintf("testfun/Peri_om%02d%02d_", m, n);
        [RP(m,n),TP(m,n)] = fdfd2D(omegs(n), thetai(m), paramsP, pathhead);
    end
end

%% Find accuracy of reflection, transmittance, and energy cons
% Angle of refraction from Snell's law
n1 = 1;
n2 = sqrt(paramsF.epsr);
thetat = asin(n1 / n2 * sin(thetai));
% Compute Fresnel relations
R_Fres = ((n1*cos(thetai)-n2*cos(thetat))./(n1*cos(thetai)+n2*cos(thetat))).^2
T_Fres = 1 - R_Fres; % Assumes energy conserved

% Compute mean and std dev for calculated R and T
RF_mean = mean(RF,2)
RF_std = std(RF, 0, 2)
TF_mean = mean(TF,2)
TF_std = std(TF, 0, 2)
EconsF = RF_mean + TF_mean

RP_mean = mean(RP,2)
RP_std = std(RP, 0, 2)
TP_mean = mean(TP,2)
TP_std = std(TF, 0, 2)
EconsP = RP_mean + TP_mean

%% Compute the relative errors
RF_relerr = abs(RF_mean - R_Fres) ./ R_Fres 
TF_relerr = abs(TF_mean - T_Fres) ./ T_Fres 

RP_relerr = abs(RP_mean - R_Fres) ./ R_Fres 
TP_relerr = abs(TP_mean - T_Fres) ./ T_Fres 

EcF_relerr = abs(EconsF - 1)
EcP_relerr = abs(EconsP - 1)

%% Save the relevant data
finaldat = table(thetai, RF_mean, RP_mean, R_Fres, EcF_relerr, EcP_relerr, RF_relerr, RP_relerr, TF_relerr, TP_relerr)
writetable(finaldat, "testfun/RefDatFloqComp.txt");