%%
reftab=readtable("testfun/RefDatFloqComp.txt")

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