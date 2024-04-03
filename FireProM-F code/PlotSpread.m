% Plot comparison radius of spread from central ignition point

time =     [0   5   10     15     20    25    30];

radius_FAR = [0  5.91 12.31 19.67 27.71 34.49 40.50]; % from FARSITE on tall grass model 3
Spread = figure('units','normalized','outerposition',[0 0 1 1]);
figure(Spread)
plot(time,radius_FAR,'-*r')
xlabel('simulation time [min]','fontsize',16)
ylabel('radius of spread [m]','fontsize',16)

time = [-3 time]; % 3 min ignition time before comparison with FARSITE
radius_MOD = [0 0.98 6.09 12.17 19.13 26.52 33.26 38.91];
figure(Spread)
hold on
plot(time,radius_MOD,'-ob')
legend('FARSITE','2D Model');
set(gca,'FontSize', 16)
