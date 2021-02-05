%%
R = readmatrix("data/reflectivity.csv");
omegs = R(:,1);

lambs0 = 6e8 .* pi ./ omegs;
lambsx = lambs0 .* sin(thetai);
lambsy = lambs0 .* cos(thetai);

%f = figure('visible','off');

for j = 1:Nf
    omeg = omegs(j);
    lambx = lambsx(j);
    lamby = lambsy(j);
    
    xfile = sprintf("data/X_%02d.csv", j);
    yfile = sprintf("data/Y_%02d.csv", j);
    efile = sprintf("data/Ez_%02d.csv", j);
    x = 1e6*readmatrix(xfile);
    y = 1e6*readmatrix(yfile);
    Ez = readmatrix(efile);
    Ez = reshape(Ez, length(x), length(y));

    f = figure(j)
    [M,c] = contourf(x, y, real(Ez'), 64);
    set(c,'LineColor','none')
    xlabel('x [{\mu}m]');
    ylabel('y [{\mu}m]');
    colorbar("EastOutside");
    ylim([y(1) y(end)]);
    xlim([x(1) x(end)]);
    hold on;
    %lims = max(Ez(i,:));
    %lims = 1;
    plot([x(sfr),x(sfr)], [y(1),y(end)], '--r', "linewidth", 2.5)
    plot([x(p),x(p)], [y(1),y(end)], '-k', "linewidth", 2.5)
    %plot([x(q),x(q)], [y(1),y(end)], '-k', "linewidth", 2.5)
    %plot([z(p),z(p)], [-lims,lims], '-k', "linewidth", 2.5)
    %plot([z(q),z(q)], [-lims,lims], '-k', "linewidth", 2.5)
    %plot(z,Q(1:floor(end/2),1:floor(end/2)) * real(fsrc(z))', "-r;Re[f(z)];");
    %legend("hide");
    %plot([x(sfr+20), x(sfr+20)+2*pi/ki_x*1e6], [y(20), y(20)+2*pi/ki_x*1e6], '-r', "linewidth", 2.5)
    quiver(x(sfr+20), y(20), lamby*1e6, lambx*1e6, '-r', "linewidth", 2.5)

    freq_text = sprintf("%.2e", omeg);
    title(strcat('E_z(x,y) for \omega = ', freq_text));
    %text(0.025, 0.05, ['\omega = ' freqtext], "units", "normalized");
    %text(z(sfr+10), -0.9, '\epsilon = \epsilon_0');
    %text(z(p+50), -0.9, '\epsilon = 5\epsilon_0')
    %text(z(q+10), -0.9, '\epsilon = \epsilon_0')
    
    figfile = sprintf("data/Ez_%02d.png", j);
    saveas(f, figfile)
    %system ("wslview figs/test.png")
end