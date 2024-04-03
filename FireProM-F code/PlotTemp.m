function [] = PlotTemp(X, Y, T, FigT, time, Model, TT)

% Author: Mauro S. Innocente
% Date: 11/01/2018
% This function plots the colormap for temperature field.

Tmin = min(min(min(TT)));
Tmax = max(max(max(TT)));

figure(FigT)

pcolor(X,Y,T);

colormap jet
% colormap parula
% shading interp
% shading faceted
shading interp
temp = ['Temperature [K] - ' Model ' Fire Model: t = ' num2str(time) ' seconds'];
title(temp,'fontsize',14)
xlabel('x [m]','fontsize',14)
ylabel('y [m]','fontsize',14)
axis equal
axis([min(X(1,:)), max(X(1,:)), min(Y(:,1)), max(Y(1,:))])
caxis([Tmin Tmax])
% caxis([295 1295])                                                        % Paolo's limits.
colorbar('EastOutside')
set(gca,'FontSize', 14)
hold on


end