% MainTestFire
% Function used to test the FireProM-F_02 model

% Authors: Mauro S. Innocente     Paolo Grasso
% Date:    06/12/2017             05/04/2022
% Wrapper function.

% CONTENTS:
% M.1  Settings
% M.2  Initialisation of FireProM-F
% M.3  Initialisation of figures
% M.4  Simulation
% M.5  Display
% M.6  Save

clear all
close all
clc
rng('Default')


divide = 60; % divisions in 1 min    !!! CHANGE HERE (if needed) !!!
% It can be as low as 37 which results in a faster simulation,
% but it is advisable to use 60 (i.e. 1 s time step)
Dt = 60/divide; % timestep in seconds

%% (M.1) SETTINGS %% SETTINGS %% SETTINGS %% SETTINGS %% SETTINGS %% SETTINGS %% 

%%% !!! CHANGE FROM HERE ... %%%
% filename formats examples:
% Temperature field: Temper_#ms_###deg_##min
% Fuel Energy field: Energy_#ms_###deg_##min
filename1 = 'Temper_6ms_90deg_5min.gif';      % !!! CHANGE HERE !!!
filename2 = 'Energy_6ms_90deg_5min.gif';      % !!! CHANGE HERE !!!



% % @@@ define the following as a nx by ny matrix for a frozen wind
% simualtion
    U_w=ones(100,100,4)*4.6;  % [m/s] reliable max is 6 m/s, altough it can work also for stronger winds but the solution deteriorates
    V_w=ones(100,100,4)*0;
%  

% importing velocity components from source (FDS here) for dynamic wind
% simuations-25 second update interval here
% u1=xlsread('C64-0s','C:D');   U_w(:,:,1) = reshape(u1(:,1),[100,100]);    V_w(:,:,1) = reshape(u1(:,2),[100,100]);
% u2=xlsread('C64-25s','C:D');  U_w(:,:,2) = reshape(u2(:,1),[100,100]);    V_w(:,:,2) = reshape(u2(:,2),[100,100]);
% u3=xlsread('C64-50s','C:D');  U_w(:,:,3) = reshape(u3(:,1),[100,100]);     V_w(:,:,3) = reshape(u3(:,2),[100,100]);
% u4=xlsread('C64-75s','C:D');  U_w(:,:,4) = reshape(u4(:,1),[100,100]);     V_w(:,:,4) = reshape(u4(:,2),[100,100]);

%
%
%
tic
timeofsimulation =2; % [min]


ignition_time = 0*Dt; % 

ignition_steps = ceil(ignition_time/Dt);
SimulationLength = timeofsimulation*divide + ignition_steps;   % number of time steps
fprintf('ignition steps = %i    sim. steps = %i\n',ignition_steps,SimulationLength);
%==========================================================================
% BOUNDS OF SEARCH-SPACE.
Width_Space = 100;                  % [m]  !!! CHANGE HERE !!!
xmin = zeros(1,2);
xmax = ones(1,2)*Width_Space;       % [m];
%==========================================================================
% MESH.
nx = 100;     % Number of nodes in x.  !!! CHANGE HERE !!!
ny = 100;     % Number of node sin y.  !!! CHANGE HERE !!!

Dx = (xmax(1,1) - xmin(1,1))/(nx);   % [m];
Dy = (xmax(1,2) - xmin(1,2))/(ny);   % [m].
x = linspace(xmin(1,1), xmax(1,1), nx);
y = linspace(xmin(1,2), xmax(1,2), ny);
[X,Y] = meshgrid(x,y);
Y = flipud(Y); % Changes the order of Y so that is increasing to the north.
%==========================================================================
% WIND auto setting
% orientation = deg2rad(angles); % it can be a matrix % conversion from [deg] to [rad]
reduce_wind = 0.137;  % wind reduction coefficient- use as a calibration coefficient


%adding Gaussian noise to the velocity components----
 Rn=0.0;    %noise ratio-default to zero
 U_w = rot90(((Rn*randn(nx,nx).*U_w)+U_w).*reduce_wind,3);
 V_w = rot90(-((Rn*randn(nx,nx).*V_w)+V_w).*reduce_wind,3);
%%%%%%%%%%    wind matrix from FDS should rotate 270 degrees to become aligned with Temperature matrix in FirePromF

%%% end SETTINGS %%% end SETTINGS %%% end SETTINGS %%% end SETTINGS %%% end SETTINGS 

%% (M.2) INITIALISATION OF FIRE-SPREAD MODEL
model = 2;       % ENTER the code for the Fire Model desired. Suggested to use model 2 = FireProM-F
FlagF = true;    % Auxiliary variable for initialisation of fire and fuel .gif. automatically turns to false once the initialisation is completed.
disp('###################################################################')
fprintf('Computations\n\n')
switch model
    case 1 % not sure it still works (suggested not to use)
        %-----------------------------
        fprintf('Simplified Fire-Spread Model\n\n')
        Model = 'Simplified';
        Tamb = 21 + 273.15;  % [K] Ambient temperature
        Tig = 220 + 273.15;  % [K] Ignition temperature
        kappa = 100;                         % [W/m/K].
        rhoGas = 1.853;                      % [kg/m3].
        ThicknessGas = 1;                    % [m].
        MassGas = rhoGas*Dx*Dy*ThicknessGas; % [m3].
        cp = 1670;                           % [J/kg/K].
        alpha = kappa/cp/rhoGas;             % [m2/s].
        % alpha = 1.11 * 10^-4;              % [m2/s].
        %-----------------------------
        % Stability check (I need to make sure the check is correct).
        Dt = min(Dt,(Dx^2+Dy^2)/8/alpha);
        %-----------------------------
        qdotref = 10^5;                      % [W/m2].
        h = 40;                              % [W/m2/K] (convection coefficient).
        %-----------------------------
        % RADIATION (IGNORED at present!!)
        epsilon = 0;                         % [] (emissivity factor).
        sigma = 5.67*10^-8;                  % [W/m2/K4] (Stefan-Botlzmann constant).
        %-----------------------------
        [T, E] = FireModel1(FlagF, nx, ny, Dx, Dy, Dt, Tamb, Tig, kappa, alpha, qdotref, h, epsilon, sigma, [], []);
        %-----------------------------
    case 2 % modified 06/02/2018 Paolo, modified 05/04/2022
        %-----------------------------
        fprintf('  FireProM-F v.02\n')
        fprintf('  Progress:   0 %%')
        Model = 'Efficient';
        rhoGas = 1.853;         % [kg/m3] estimated gas mixture density
        ThicknessGas = 0.3;     % [m] ought to the the representation of the fuel available in gaseous form
        % the energy is function of the gas volume (so thickness is needed)
        
        MassGas = rhoGas*Dx*Dy*ThicknessGas;   % [m3]
        
        FireState.FLPerc = 100*ones(nx,ny); % percentage distribution of the fuel, 0% means no fuel, can use intermediate values
        % example random FLPerc distribution:
%       FireState.FLPerc = 100*rand(nx,ny); % uncomment this if want to try random distribution of one fuel

        FireCoeff.Tamb = 21 + 273.15; % [K]
        
        %%% DEFINE HERE THE IGNITION DESIRED !!! CHANGE HERE !!!
        % case single ignition in the centre
        FireCoeff.Xig = [0];
        FireCoeff.Yig = [Width_Space/2];
        % parabolic ingnition (many points distributed along a parabola) just for fun
%         FireCoeff.Xig = linspace(40,60,80);
%         FireCoeff.Yig = 50-((FireCoeff.Xig-50).^2)/10;


        FireCoeff.nx = nx;
        FireCoeff.ny = ny;
        FireCoeff.Dx = Dx;
        FireCoeff.Dy = Dy;
        FireCoeff.Dt = Dt;
        
        % !!! CHANGE FROM HERE ... %%%
        FireCoeff.SingleFM = 'Yes'; % 'Yes' if only one Fuel Model; 'No' if distribution on different fuel types. Please keep Yes.
        FireCoeff.FM = 1; % Fuel Model #ID (=1,2,6) or matrix distribution of #IDs if chose 'No' at the previous line
        FireCoeff.posSpy = [52 50];    % position where to plot the 1D dynamics (e.g. fuel consumption rate)
        % ... TO HERE !!! %%%%%%%%%%%%
        
        %FireState.wind = [0,0];        % initialisation of x and y components of the wind (0 m/s magnitude at the beginning)
        % @@@ if matrix of wind (you still need the two components, so you will have two matrices defined as follows:
        FireState.wind(:,:,1) = zeros(nx,ny); % x components on the mesh,initialised to zero
        FireState.wind(:,:,2) = zeros(nx,ny); % y components on the mesh,initialised to zero
        FireCoeff.alpha = zeros(nx,ny);
        % @@@ it can work also if matrix of different orientations of wind
        
        % Calling FireProM-F for the first time for initialisation (i.e.
        % with FlagF = true. Then it will turn FlagF = false, so that the
        % fire model will know that it has been already initialised)
        [FireState, FireCoeff, FireSpy] = FireProMF_02(FlagF, ThicknessGas, FireState, FireCoeff);
        store_Tspy  = FireSpy.T;
        store_Xfspy = FireSpy.Xf;
        store_XCspy = FireSpy.XC;
        
        % case mixed fuel (i.e. when a mixture of fuels is considered, FireCoeff.SingleFM = 'No';)
        % might not work !!! NEED TO CHECK !!!
%         [T, E, Ab, Pmax, store_Tspy, store_Xfspy, store_XCspy] = FireProMF01_mix(FlagF, nx, ny, Dx, Dy, Dt, Tamb, ThicknessGas, sigma, [], 0*u, 0*v,posx,posy);
        %-----------------------------
        % to here (*Paolo)
    case 3
        fprintf('High-Fidelity Fire-Spread Model\n\n')
        %-----------------------------
end
FlagF = false; % from now on, when calling FireProM-F, the function will know that it is no more initialisation mode

% Initialisation of the history storage variables
TT = FireCoeff.Tamb * ones(ny, nx, SimulationLength + 1);  % history of Temperature.
EE = TT;                        % history of Energy.
TT(:,:,1) = FireState.Temp;     % initial T state
EE(:,:,1) = FireState.Enrg;     % initial E state
AB(1)     = FireState.Ab;       % 1st element history of burnt area


%% (M.3) INITIALISATION OF FIGURES & GIFS
%%% DISPLAY %%% DISPLAY %%% DISPLAY %%% DISPLAY %%% DISPLAY %%% DISPLAY %%%

FigT = figure('units','normalized','outerposition',[0.25 0.1 0.45 0.85]);
FigE = figure('units','normalized','outerposition',[0.25 0.1 0.45 0.85]);
FigAb = figure('units','normalized','outerposition',[0.25 0.1 0.45 0.85]);
FigR = figure('units','normalized','outerposition',[0.25 0.1 0.45 0.85]);
% FigT = figure();   % if want to use your default Matlab settings
% FigE = figure();   % if want to use your default Matlab settings
% FigAb = figure();  % if want to use your default Matlab settings
% FigR = figure();   % if want to use your default Matlab settings

% calling function that plots the temperature field
PlotTemp(X, Y, FireState.Temp, FigT, 0-ignition_time, Model, TT)

% calling function that plots the fuel energy field
PlotEnergy(X, Y, FireState.Enrg, FigE, 0-ignition_time, Model, EE)

% formatting the temperature figure
figure(FigT)
set(gcf,'color','w');
% saving the frame of the temperature in an animated gif
frame = getframe(gcf);
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,filename1,'gif','LoopCount',1,'DelayTime',0);

% formatting the energy figure
figure(FigE)
set(gcf,'color','w');
% saving the frame of the energy in an animated gif
frame = getframe(gcf);
im = frame2im(frame);
[A,map] = rgb2ind(im,256);
imwrite(A,map,filename2,'gif','LoopCount',1,'DelayTime',0);
%%% DISPLAY %%% DISPLAY %%% DISPLAY %%% DISPLAY %%% DISPLAY %%% DISPLAY %%%



%% (M.4) SIMULATION %%% SIMULATION %%% SIMULATION %%% SIMULATION %%% SIMULATION 
for ts=1:SimulationLength
    if model==1

      [T, E] = FireModel1(FlagF, nx, ny, Dx, Dy, Dt, Tamb, Tig, kappa, alpha, qdotref, h, epsilon, sigma, T, E);
    else
        
        % modified (*) 06/02/2018 Paolo, 05/04/2022 Paolo
        if ts<2 % ignition_steps + 1
         % addtional 1st step with no wind
         FireState.Temp(nx/2+1,1)=550;
         FireState.Temp(nx/2-1,1)=550;
         [FireState, FireCoeff, FireSpy] = FireProMF01_02(FlagF, ThicknessGas, FireState, FireCoeff);
         store_Tspy(ts+1)  = FireSpy.T;
         store_Xfspy(ts+1) = FireSpy.Xf;
         store_XCspy(ts+1) = FireSpy.XC;
         % case mixed fuel
%          [T, E, Ab, Pmax, store_Tspy(ts+1), store_Xfspy(ts+1), store_XCspy(ts+1)] = FireProMF01_mix(FlagF, nx, ny, Dx, Dy, Dt, Tamb, ThicknessGas, sigma, T, 0*u, 0*v,posx,posy);

        else
            
        % wind field definition  ----here for 25 second update interval  
        if ts<=25 %0-25s
           FireState.wind(:,:,1)= U_w(:,:,1); % matrix of x components of wind [m/s]
           FireState.wind(:,:,2)= V_w(:,:,1); % matrix of y component of wind [m/s]
         elseif ts<=50   
           FireState.wind(:,:,1)= U_w(:,:,2); % matrix of x components of wind [m/s]
           FireState.wind(:,:,2)= V_w(:,:,2); % matrix of y component of wind [m/s]
        elseif ts<=75  
           FireState.wind(:,:,1)= U_w(:,:,3); % matrix of x components of wind [m/s]
           FireState.wind(:,:,2)= V_w(:,:,3); % matrix of y component of wind [m/s]
        elseif ts<=100   
           FireState.wind(:,:,1)= U_w(:,:,4); % matrix of x components of wind [m/s]
           FireState.wind(:,:,2)= V_w(:,:,4); % matrix of y component of wind [m/s]
        end
            
        %Defining the ignition line-just increase the temperature and
        %ignition will ontinue
        if ts<26
            i=ts;
            FireState.Temp(i+50,1)=550;
            FireState.Temp(50-i,1)=550;           
        end
   
         [FireState, FireCoeff, FireSpy] = FireProMF01_02(FlagF, ThicknessGas, FireState, FireCoeff);
         store_Tspy(ts+1)  = FireSpy.T;
         store_Xfspy(ts+1) = FireSpy.Xf;
         store_XCspy(ts+1) = FireSpy.XC;
         % case mixed fuel
%          [T, E, Ab, Pmax, store_Tspy(ts+1), store_Xfspy(ts+1), store_XCspy(ts+1)] = FireProMF01_mix(FlagF, nx, ny, Dx, Dy, Dt, Tamb, ThicknessGas, sigma, T, u, v,posx,posy);
        end
        % progress counter update
        perc = round(ts/SimulationLength*100); % percentage completion
        dp = 1; % visualised percentage increment - e.g. 5 means 5% increment
        if (mod(perc,dp)==0)&&(perc>=dp)
          for jperc=0:1:(floor(log10(perc))+2)
              fprintf('\b'); % delete previous counter display
          end
          fprintf('%d %%', perc);
%           pause(.2); % allows time for display to update
        end
    end
    
    % History of Temperature, Energy and Positions.
    TT(:,:,ts+1) = FireState.Temp;
    EE(:,:,ts+1) = FireState.Enrg;
    AB(ts+1)     = FireState.Ab;
    
end 
fprintf('  Done\n')
% if exist('store_XCspy','var')==1
store_XCspy = store_XCspy/store_XCspy(end);
Fig_Dyn = figure('units','normalized','outerposition',[0.15 0.1 0.65 0.85]);
figure(Fig_Dyn)
plot(Dt*(0:SimulationLength)/60,store_Tspy,'-r','Linewidth',2);
hold on
plot(Dt*(0:SimulationLength)/60,store_Xfspy,'-b','Linewidth',2);
hold on
plot(Dt*(0:SimulationLength)/60,store_XCspy,'-.m','Linewidth',2);
hold on
mass_loss_rate = -(store_Xfspy(2:end)-store_Xfspy(1:end-1));
mass_loss_rate = mass_loss_rate/(max(mass_loss_rate));
plot(Dt*((1:SimulationLength)-0.5)/60,mass_loss_rate,':k','Linewidth',2);
legend('Temperature','fuel mass','CO_2 mass','fuel mass rate')
textitle = ['Transient at [' num2str(FireCoeff.posSpy(1)) ' ' num2str(FireCoeff.posSpy(2)) '] m'];
title(textitle,'fontsize',14)
xlabel('time [min]','fontsize',14)
ylabel('normalised variables','fontsize',14)
set(gca,'FontSize', 14)
grid on
% end

fprintf('\n  Results showcase:\n')
fprintf('    Total Burnt Area: %d [m^2]\n',FireState.Ab);
fprintf('    Maximum Power Output: %.2s [kW]\n',FireState.Pmax);
fprintf('    Maximum peak Temperature: %.2s [K]\n',max(max(max(TT))));
fprintf('\n');

%%% SIMULATION %%% SIMULATION %%% SIMULATION %%% SIMULATION %%% SIMULATION 
toc

%% (M.5) DISPLAY %%% DISPLAY %%% DISPLAY %%% DISPLAY %%% DISPLAY %%% DISPLAY %%%
disp('###################################################################')
fprintf('Plots of Temperature Field (Fire Only)\n')
tic
f=1;
for ts=0:3:SimulationLength
    time = fix(ts*Dt);
    T = TT(:,:,ts+1);
    figure(FigT)
    PlotTemp(X, Y, T, FigT, time-ignition_time, Model, TT)
    
    %---plotting wind vectors-----
%      if ts<=10*f
%      quiver(rot90(U_w(:,:,f)),rot90(-V_w(:,:,f)))  
%      else 
%      f=f+1;
%      quiver(rot90(U_w(:,:,f)),rot90(-V_w(:,:,f))) 
%      end
     %-------------------------------   
%     quiver(U_w(:,:,4),V_w(:,:,4))
    
    set(gcf,'color','w');
    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    imwrite(A,map,filename1,'gif','WriteMode','append','DelayTime',0);  
end; clear ts
toc
disp('###################################################################')
fprintf('Plots of Energy Field (Fire Only)\n')
tic
for ts=0:5:SimulationLength
    time = fix(ts*Dt);
    E = EE(:,:,ts+1);
    figure(FigE)
    PlotEnergy(X, Y, E, FigE, time-ignition_time, Model, EE)
    set(gcf,'color','w');
    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    imwrite(A,map,filename2,'gif','WriteMode','append','DelayTime',0);  
end; clear ts

toc
disp('###################################################################')
fprintf('Plot of Burnt Area (Fire Only)\n')
%%
PlotBurnt(AB, Dt, SimulationLength, ignition_time, Model, FigAb)
fprintf('Plot of Radius of the fire front (Fire Only)\n')
PlotRadius(AB, Dt, SimulationLength, ignition_time, Model, FigR)

%%% DISPLAY %%% DISPLAY %%% DISPLAY %%% DISPLAY %%% DISPLAY %%% DISPLAY %%%

%% (M.6) SAVE DATA %% SAVE DATA %% SAVE DATA %% SAVE DATA %% SAVE DATA %% SAVE
tspace = ((ignition_time*60:Dt:SimulationLength*Dt))/60-ignition_time;
ll = length(tspace);

filename='Radius.txt';
fID = fopen(filename,'w');
fprintf(fID,'Radius [m]\n');
Radius = (AB/pi).^(1/2);
n=length(AB);
for i=1:ll
fprintf(fID,'%7.2f \n',Radius(i+n-ll));
end
fclose(fID);
%
filename='Area.txt';
fID = fopen(filename,'w');
fprintf(fID,'Area [m^2]\n');
n=length(AB);
for i=1:ll
fprintf(fID,'%7.2f \n',AB(i+n-ll));
end
fclose(fID);
%
filename='Time.txt';
fID = fopen(filename,'w');
fprintf(fID,'Time [s]\n');
n=length(AB);
for i=1:ll
fprintf(fID,'%7.2f \n',tspace(i));
end
fclose(fID);
