function [] = PlotBurnt(AB, Dt, SimulationLength, ign_time, Model, FigAb)

% Author: Mauro S. Innocente
% Date: 11/01/2018
% This function plots the grey colormap for energy field.

Amax = max(max(max(AB)));

figure(FigAb)
set(gca,'FontSize', 16)
tspace = (ign_time:Dt:SimulationLength*Dt)-ign_time;
ll = length(tspace);
plot(tspace/60,AB(end-ll+1:end),'o-b','LineWidth',2);
temp = ['Burnt Area [m^2] - ' Model ' Fire Model'];
title(temp,'fontsize',16)
xlabel('time [min]','fontsize',16)
ylabel('burnt area [m^2]','fontsize',16)
grid on
pbaspect([1 1 1])
axis([0, tspace(end)/60, 0, ceil(Amax)])
set(gca,'FontSize', 16)

end