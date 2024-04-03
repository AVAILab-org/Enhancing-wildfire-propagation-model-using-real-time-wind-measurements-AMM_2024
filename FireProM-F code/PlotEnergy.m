function [] = PlotEnergy(X, Y, E, FigE, time, Model, EE)

% Author: Mauro S. Innocente
% Date: 11/01/2018
% This function plots the grey colormap for energy field.

Emin = min(min(min(EE)));
Emax = max(max(max(EE)));

figure(FigE)
set(gca,'FontSize', 14)
% E(4,3)
% E(4,4)
% E(5,3)
pcolor(X,Y,E);
colormap gray
%colormap jet
% colormap parula
shading interp
% shading faceted
% shading flat
temp = ['Fuel Energy [J/m2] - ' Model ' Fire Model: t = ' num2str(time) ' seconds'];
title(temp,'fontsize',14)
xlabel('x [m]','fontsize',14)
ylabel('y [m]','fontsize',14)
axis equal
axis([min(X(1,:)), max(X(1,:)), min(Y(:,1)), max(Y(1,:))])
caxis([Emin Emax])
colorbar('EastOutside')
set(gca,'FontSize', 14)


end