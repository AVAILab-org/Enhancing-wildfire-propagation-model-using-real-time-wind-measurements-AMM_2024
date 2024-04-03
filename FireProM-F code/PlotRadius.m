function [] = PlotRadius(AB, Dt, SimulationLength, ign_time, Model, FigR)

% Author: Mauro S. Innocente
% Date: 11/01/2018
% This function plots the grey colormap for energy field.

figure(FigR)
set(gca,'FontSize', 16)
tspace = (ign_time:Dt:SimulationLength*Dt)-ign_time;
ll = length(tspace);
Radius = (AB/pi).^(1/2);
plot(tspace/60,Radius(end-ll+1:end),'o-b','LineWidth',2);
temp = ['Radius of the fire front - ' Model ' Fire Model'];
title(temp,'fontsize',16)
xlabel('time [min]','fontsize',16)
ylabel('radius [m]','fontsize',16)
grid on
pbaspect([1 1 1])
axis([0, tspace(end)/60, 0, ceil(max(Radius))])
set(gca,'FontSize', 16)

end