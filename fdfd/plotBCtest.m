%%
reftab=readtable("bigtest/RefDatFloqComp.txt")

%%
figure(1);
plot(reftab.thetai, reftab.RF_relerr, reftab.thetai, reftab.RP_relerr)
legend('Floquet','Periodic')
title('Relative error in reflectance')
xlabel('Angle of incidence')
saveas(gcf,"figs/FloqPeriComp_EnergyCons.png");

figure(2);
plot(reftab.thetai, reftab.EcF_relerr, reftab.thetai, reftab.EcP_relerr)
legend('Floquet','Periodic')
title('Relative error in conservation of energy')
xlabel('Angle of incidence')
saveas(gcf,"figs/FloqPeriComp_RefError.png");

%%
figure(3);
clf;
plot(reftab.thetai, reftab.RF_mean, reftab.thetai, reftab.RP_mean)
hold on;
plot(reftab.thetai, reftab.R_Fres, '--k');
legend('Floquet','Periodic','Fresnel')
title('Calculated reflectance compared to Fresnel calculation')
xlabel('Angle of incidence')
saveas(gcf,"figs/FloqPeriComp_Reflectance.png");

%%
figure(4);
clf;
plot(reftab.thetai(1:end-1), reftab.RF_mean(1:end-1))
hold on;
plot(reftab.thetai(1:end-1), reftab.R_Fres(1:end-1), '--k');
legend('Floquet','Fresnel')
title('Calculated reflectance compared to Fresnel calculation')
xlabel('Angle of incidence')
saveas(gcf,"figs/Floq_Reflectance.png");

%%
figure(5);
semilogy(reftab.thetai, reftab.RF_relerr, reftab.thetai, reftab.TF_relerr, reftab.thetai, reftab.EcF_relerr)
legend('Reflectance', 'Transmittance', 'Energy')
title('Relative error for Floquet boundary condition')
xlabel('Angle of incidence')
saveas(gcf,"figs/Floq_RefandEnergyCons.png");
