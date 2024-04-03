function [FireState, FireCoeff, FireSpy] = FireProMF01_02(FlagF, ThicknessGas, FireState, FireCoeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% FireProMF01 is FireProM-F v.01 function
% written for Matlab R2017a
% (should work fine also for older Matlab versions)
%
% Coventry, 03 August 2020
% authors:
% Dr. Paolo Grasso       (paologk90@gmail.com)
% Dr. Mauro S. Innocente (Mauro.S.Innocente@coventry.ac.uk)
%
% Please refer to the following journal paper:
% "Physics-based model of wildfire propagation towards
%  faster-than-real-time simulations"
% https://doi.org/10.1016/j.camwa.2020.05.009
% Thanks
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [OUTPUT] = FireProMF01(INPUT)
%
% INPUT = (FlagF,FireState, FireCoeff, ) * * * * * * * * * * * * * * * * * * * * 
%
%  FlagF        : 0 if requires inizialisation                [1]
%
%  ThicknessGas : thickness of the mixing layer               [m]
%
%  FireState                                          [structure]  of: ...                                                       
%     (i.e. state variables -- Note: need to provide the initial condition    
%      - i.e. initial state) This variable is a structure so access to         
%      the inner variable through 'FireState.xxxx' - e.g. FireState.Temp.      
%     (variables starting with Capital letters are mostly matrices)                 
%   
%   .FLPerc     : Fuel Load Percentage (provide I.C)          [%]
%   .wind       : [u v] asymptotic components of wind (B.C.)  [m/s]
%   .Xs(:,:,id) : Mole fraction of chemical species           [1]
%                 id = 1,...,5 for CH4, O2, CO2, H20 and 'air' respectively
%                 where CH4 is the representation of the fuel
%   .M          : molar mass distribution                     [kg/mol]
%   .cp         : constant pressure coeficcient of mixture    [J/kg/K]
%   .Temp       : temperature matrix (provide I.C.)           [K]     
%   .Tmax       : maximum temperature (monitor variable)      [K]
%   .hc         : specific combustion enthalpy                [kJ/kg]
%   .Enrg       : energy density matrix                       [J/m^2]
%   .Ab         : total burnt area                            [m^2]
%   .Pmax       : maximum output power (monitor variable)     [kW]
%
%  FireCoeff (i.e. the coefficients)                  [structure]  of: ...
%   .nx, .ny    : number of cells centres in x and y dir      [1]
%   .Dx, .Dy    : cell size in x and y directions             [m]
%   .Dt         : time step size                              [s]
%   .Xig, .Yig  : vectors of ignition (ig) points [X(i),Y(i)] [m]
%   .SingleFM   : 'Yes' if only one fuel model, 'No' else     [string]
%   .FM         : distribution of fuel model IDs (or only one value)  [#ID]
%   .posSpy     : [posx posy] position on scope point         [m m]
%   .alpha      : orientation of asymptotic wind              [rad]
%   .R          : universal gas constant             = 8.3140 [J/mol/K]                      
%   .Tref       : reference temperature              = 298.15 [K]
%   .Pref       : reference pressure                 = 101325 [Pa]
%   .sigma      : Stefan-Botlzmann constant       = 5.6704e-8 [W/m^2/K^4]
%   .coef       : stoichiometric coefficient of CH4 comb.     [1]
%   .ms         : molar mass of each species (s)              [kg/mol]
%   .Hrefs      : molar formation enthalpies (s)              [J/mol]
%   .cps        : constant pressure coefficients (s)          [J/kg/K]
%   .Xe         : fuel extintion molar fraction               [1]
%   .Xo2e       : oxigen extintion molar fraction             [1]
%   .overT      : T peak of ignition gaussian distributions   [K]
%   .FLm        : flux limiter coefficient                    [1]
%   .FModels    : collection of calibration coefficients for various fuel
%                 types by rows = [Ch Ar Kth Ot Dz Ca Tig; ........ ; ...]
%   .Ch         : distribution or single value of enthalpy    [1]
%                 correction coefficient
%   .Ar         : Arrhenius pre-exponential coefficient       [~]
%   .Kth        : thermal conductivity constant               [W/m/K]
%   .Ot         : optical thickness in 2D direction           [m]
%   .Dz         : optical thickness in vertical direction     [m]
%   .Ca         : turbulent convenction coefficient
%   .Tig        : autoignition temperature                    [K]
%   .Ta         : activation temperature                      [K]
%   .epsilon    : emissivity coefficient                      [~]
%   .corr       : correction coefficient for power output     [1]
%   .fl_cases   : vector of typical fuel loads for each FM    [1]
%   .Fl         : fuel load distribution (or just one value)  [1]
%   .Xfref      : reference CH4 mole fraction of the mixture  [1]
%   .f0         : boundary percentages of forced 0 fuel       [%]
%   .Xf0        : initial distribution of fuel mole fraction  [1]
%   .rho0       : initial average density                     [kg/m^2]
%   .cp0        : initial constant pressure coefficient       [J/kg/K]
%
%
% OUTPUT: [FireState, FireCoeff, FireSpy] * * * * * * * * * * * * * * * * *
%
%  FireState  (the update for the following time step ~ t = t + dt)
%
%  Fire Coeff (the complete set of coefficients
%
%  FireSpy                                            [structure]  of: ...
%        (intended for scoping values on particular places in the domain)
%   .T          : temperature profile                         [K]
%   .Xf         : fuel consumption                            [1]
%   .XC         : CO2 formation                               [1]
%    !          : could add more monitor variables            [!]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rng(11503)
persistent Xf
% global Y0D
% global FlD
% global hc0 
global r

nx = FireCoeff.nx;
ny = FireCoeff.ny;
Dx = FireCoeff.Dx;
Dy = FireCoeff.Dy;
Dt = FireCoeff.Dt;

% (1) FIRST STEP INITALISATION
if FlagF 
    %==========================================================================
    % MISCELLANEA PROPERTIES
    FireCoeff.R = 8.314; % gas constant [J/mol/K]   
    FireCoeff.Tref  = 298.15; % [K] 
    FireCoeff.Pref  = 101325; % [Pa]
    FireCoeff.sigma = 5.6704e-8; % [W/m^2/K^4] Boltzman constant
    
    %                   CH4         O2          CO2          H20         air
    FireCoeff.coef =  [ -1          -2          1            2                   ];
    FireCoeff.ms   =  [ 16.043      32          44.011       18.016      28.9644 ]; % g/mol
    FireCoeff.ms = FireCoeff.ms*0.001;                                              % kg/mol
    FireCoeff.Hrefs = [ -74.873      0         -393.522     -241.826             ]; % kJ/mol
    FireCoeff.Hrefs = FireCoeff.Hrefs*1000;                                         % J/mol
    FireCoeff.cps  =  [ 2220        919         844          1996         1010   ]; % J/kg/K
    FireCoeff.Xe = 0.005;
    FireCoeff.Xo2e = 0.005;
    FireCoeff.overT = 300; % for gaussian distribution

    %%%%%%%%%%%%-------------------------------------------%%%%%%%%%%%%%%%%%%%%%%%%%
    FireCoeff.FLm =0.0; % put FM1 0.5     FM2 0.25     FM6  0  !!! CHANGE HERE !!! (it will be automatic in case of mixture of multiple FMs, but need to check this)
    %                           Ch      Ar     Kth   Ot
    FireCoeff.FModels      = [ 0.995 5.757e-5 1.467 0.0343 1.67 0.0492 367.9     %   FM1
                               0.988 4.324e-5 0.749 0.0365 1.49 0.0516 368.7];   %   FM2
    FireCoeff.FModels(6,:) = [ 0.803 4.154e-5 0.383 0.0298 1.91 0.0216 356.7];   %   FM6 (assigned to line 6 because models from 3 to 5 are missing)                     
    xx = FireCoeff.FModels(FireCoeff.FM,:);
    FireCoeff.Ch  = xx(1);
    FireCoeff.Ar  = xx(2); 
    FireCoeff.Kth = xx(3); % note: is has been corrected by a 1.5 factor larger than the original values
    FireCoeff.Ot  = xx(4); 
    FireCoeff.Dz  = xx(5); 
    FireCoeff.Ca  = xx(6);
    FireCoeff.Tig = xx(7); 

    FireCoeff.Ta = FireCoeff.Tig;
    FireCoeff.epsilon = 1; % for ideal black body; or change to values smaller than 1 if not ideal emittance (e.g. 0.85)
    FireCoeff.corr = 1.1; % combustion energy correction

    FireCoeff.fl_cases = [0.74 4]; % fuel loading: [FM1 FM2]
    FireCoeff.fl_cases(6) = 3.5; % extending the previous fl vector: [FM1 FM2 0 0 0 FM6]
    tuning_factor = 0.22417; % for the fuel loading
    FireCoeff.Fl = FireCoeff.fl_cases(FireCoeff.FM);
    FireCoeff.Fl = FireCoeff.Fl*tuning_factor;
    
    Xs0  = [ 0         0.232       0.008        0.001       0.759];
    FireState.Xs = zeros(nx,ny,5);
    FireCoeff.Xfref = 0.1;
    % force zero-fuel boundaries for stability of the solution when fire
    % approaches the domain boundaries
    FireCoeff.f0(1)=1+floor(nx/100); % minimum x index
    FireCoeff.f0(2)=ceil(99*nx/100); % maximum x index
    FireCoeff.f0(3)=1+floor(ny/100); % minimum y index
    FireCoeff.f0(4)=ceil(99*ny/100); % maximum y index
    for ii = FireCoeff.f0(1):1:FireCoeff.f0(2)
        for jj = FireCoeff.f0(3):1:FireCoeff.f0(4)
            FireState.Xs(ii,jj,1) = FireState.FLPerc(ii,jj)*FireCoeff.Xfref/100;
        end
    end
    FireCoeff.Xf0 = FireState.Xs(:,:,1);
    for tt=2:1:5
        FireState.Xs(:,:,tt) = Xs0(tt)*ones(nx,ny).*(1-FireState.Xs(:,:,1));
    end
    
    % mixture molar mass estimation
    FireState.M = zeros(nx,ny);
    for s=1:1:5
        FireState.M = FireState.M + FireState.Xs(:,:,s)*FireCoeff.ms(s);
    end
    Xf = FireState.Xs(:,:,1);
    % Density
    FireCoeff.rho0=(FireCoeff.Pref/FireCoeff.R)*FireState.M(1,1)*(FireCoeff.Tamb^(-1));
    % Heat capacity
    FireState.cp = zeros(nx,ny);
    for s=1:1:5
        FireState.cp = FireState.cp + FireState.Xs(:,:,s)*FireCoeff.ms(s)*FireCoeff.cps(s);
    end
    FireState.cp  = FireState.cp./FireState.M;
    FireCoeff.cp0 = FireState.cp(1,1);
    
    % INITIAL TEMPERATURE FIELD
    FireState.Temp = ones(ny,nx) * FireCoeff.Tamb;
    % Spread ignition points
    sigmaT = 0.015; % a value of 0.015 has been used for calibration, but this value can be changed if needed for wider ignition region
    sigmaT = sigmaT*nx*Dx;
    nxig = round(FireCoeff.Xig/Dx); % can be a vector of integers for multiple ignition points
    nyig = round(FireCoeff.Yig/Dx); % can be a vector of integers for multiple ignition points
    nyig = ny - nyig;
    for cx0 = 1:1:length(nxig) % iterates for every ignition point
        % summing up Gaussian ditributions of temperature with peak higher
        % than the ignition temperature by a FireCoeff.overT value
        for cx1=1:1:nx
            for cx2=1:1:ny
                dist2=(((cx1-nxig(cx0))*Dx)^2+((cx2-nyig(cx0))*Dy)^2);
                Tbox = FireCoeff.Tamb +(FireCoeff.Tig+FireCoeff.overT-FireCoeff.Tamb)*exp(-dist2/sigmaT^2);
                FireState.Temp(cx1,cx2) = max(FireState.Temp(cx1,cx2),Tbox);
            end
        end
    end
    FireState.Tmax = max(max(FireState.Temp));
    
    % Combustion enthalpy estimation (i.e. formation enthalpies of the
    % various species)
    FireState.hc = zeros(nx,ny);
    for s=1:1:4
        FireState.hc = FireState.hc - FireCoeff.coef(s)*(FireCoeff.Hrefs(s)+(FireCoeff.ms(s)*FireCoeff.cps(s))*(FireState.Temp-FireCoeff.Tref)).*(FireState.M.^(-1));
    end
    FireState.hc = FireState.hc*FireCoeff.Ch;
    r = zeros(nx,ny);
       
    %======================================================================
    % INITIAL FUEL ENERGY FIELD
    FireState.Enrg = FireState.hc.*(FireState.Xs(:,:,1))/(FireCoeff.Xfref)*FireCoeff.rho0*ThicknessGas; % [J/kg * kg/m3 * m] = [J/m2]
    FireState.Enrg = FireState.Enrg';
    count = 0;
    for i=1:1:nx
        for j=1:1:nx
            if (FireState.Xs(i,j,1)<FireCoeff.Xf0(i,j))
                count = count+1;
            end
        end
    end
    FireState.Ab = count*Dx*Dy;
    %======================================================================
    FireState.Pmax=0;
    FireSpy.T  = FireState.Temp(FireCoeff.posSpy(1),FireCoeff.posSpy(2))/FireCoeff.Tig;
    FireSpy.Xf = FireState.Xs(FireCoeff.posSpy(1),FireCoeff.posSpy(2),1)/FireCoeff.Xfref;
    FireSpy.XC = FireState.Xs(FireCoeff.posSpy(1),FireCoeff.posSpy(2),3);
    FireState.Temp=FireState.Temp';
else
    % (2) NEXT STEPS
    FireState.Temp = FireState.Temp';
    Tr = FireCoeff.rho0*FireState.Temp.*FireState.cp;
    % 4.1 "a" calculation
    [a(:,:,1), hc1] = F_Tr(FireState,FireCoeff, Tr,FireState.hc,r,FireState.cp,FireState.M); 
    a(:,:,1) = Dt*a(:,:,1);
    ra   = F_Xf(FireCoeff, Tr,FireState.cp,FireState.Xs(:,:,1),FireState.Xs(:,:,2),nx,ny);
    Xsa = FireState.Xs;
    for jj=1:1:ny
        for ii=1:1:nx
            if FireState.Xs(ii,jj,1)~=0 
            for s=1:1:4 
                if (FireState.Xs(ii,jj,1)~=FireCoeff.Xe)
                Xval = FireState.Xs(ii,jj,s)-0.5*(Dt*FireCoeff.coef(s)/FireCoeff.ms(1))*ra(ii,jj)*FireState.M(ii,jj);
                if (s==1)&&(Xval<FireCoeff.Xe); Xval=FireCoeff.Xe;%Xval=Xe; 
                    % combustion rate correction in case of extinction
                    ra(ii,jj)=(Xsa(ii,jj,s)-FireCoeff.Xe)*FireCoeff.ms(1)/(Dt*FireCoeff.coef(s)*FireState.M(ii,jj));
                end
                % record Zero Fuel condition
                Xsa(ii,jj,s) = Xval;
                end
            end
            end
        end
    end  
     a(:,:,2) = Dt*ra;
     Ma = zeros(nx,ny);
    for s=1:1:5
        Ma = Ma + Xsa(:,:,s)*FireCoeff.ms(s);
    end
    cpa = zeros(nx,ny);
    for s=1:1:5
        cpa = cpa + (FireCoeff.cps(s)*FireCoeff.ms(s))*Xsa(:,:,s)./Ma;
    end
    hca = zeros(nx,ny);
    for s=1:1:4
        hca = hca - FireCoeff.coef(s)*(FireCoeff.Hrefs(s)+FireCoeff.ms(s)*FireCoeff.cps(s)*((Tr+0.5*a(:,:,1))./cpa/FireCoeff.rho0-FireCoeff.Tref)).*(Ma.^(-1));
    end
    hca = hca*FireCoeff.Ch;

    % 4.2 "b" calculation   
    [b(:,:,1), hc2] = F_Tr(FireState,FireCoeff, Tr+0.5*a(:,:,1),hca,ra,cpa,Ma);
    b(:,:,1) = Dt*b(:,:,1);
    rb   = F_Xf(FireCoeff, Tr+0.5*a(:,:,1),FireState.cp,Xsa(:,:,1),Xsa(:,:,2),nx,ny);
    clear Xsa
    Xsb = FireState.Xs;
    for jj=1:1:ny
        for ii=1:1:nx
            if FireState.Xs(ii,jj,1)~=0
            for s=1:1:4 % only for fuel and oxigen (partial for RK steps)
                if (FireState.Xs(ii,jj,1)~=FireCoeff.Xe)
                Xval = FireState.Xs(ii,jj,s) -0.5*(Dt*FireCoeff.coef(s)/FireCoeff.ms(1))*rb(ii,jj)*FireState.M(ii,jj);
                if (s==1)&&(Xval<FireCoeff.Xe); Xval=FireCoeff.Xe;%Xval=Xe; 
                    % combustion rate correction in case of extinction
                    rb(ii,jj)=(Xsb(ii,jj,s)-FireCoeff.Xe)*FireCoeff.ms(1)/(Dt*FireCoeff.coef(s)*FireState.M(ii,jj));
                end
                % record Zero Fuel condition
                Xsb(ii,jj,s) = Xval;
                end
            end
            end
        end
    end
     b(:,:,2) = Dt*rb;
     Mb = zeros(nx,ny);
    for s=1:1:5
        Mb = Mb + Xsb(:,:,s)*FireCoeff.ms(s);
    end
    cpb = zeros(nx,ny);
    for s=1:1:5
        cpb = cpb + (FireCoeff.cps(s)*FireCoeff.ms(s))*Xsb(:,:,s)./Mb;
    end
    hcb = zeros(nx,ny);
    for s=1:1:4
        hcb = hcb - FireCoeff.coef(s)*(FireCoeff.Hrefs(s)+FireCoeff.ms(s)*FireCoeff.cps(s)*((Tr+0.5*b(:,:,1))./cpb/FireCoeff.rho0-FireCoeff.Tref)).*(Mb.^(-1));
    end
    hcb = hcb*FireCoeff.Ch;
    
    % 4.3 "c" calculation 
    [c(:,:,1), hc3] = F_Tr(FireState,FireCoeff, Tr+0.5*b(:,:,1),hcb,rb,cpb,Mb);
    c(:,:,1) = Dt*c(:,:,1);
    rc   = F_Xf(FireCoeff, Tr+0.5*b(:,:,1),FireState.cp,Xsb(:,:,1),Xsb(:,:,2),nx,ny);
    clear Xsb hcb rb cpb  
    
    Xsc = FireState.Xs;
    for jj=1:1:ny
        for ii=1:1:nx
            if FireState.Xs(ii,jj,1)~=0
            for s=1:1:4 % only for fuel and oxigen (partial for RK steps)
                if (FireState.Xs(ii,jj,1)~=FireCoeff.Xe)
                Xval = FireState.Xs(ii,jj,s) -(Dt*FireCoeff.coef(s)/FireCoeff.ms(1))*rc(ii,jj)*FireState.M(ii,jj);
                if (s==1)&&(Xval<FireCoeff.Xe); Xval=FireCoeff.Xe;%Xval=Xe; 
                    % combustion rate correction in case of extinction
                    rc(ii,jj)=(Xsc(ii,jj,s)-FireCoeff.Xe)*FireCoeff.ms(1)/(Dt*FireCoeff.coef(s)*FireState.M(ii,jj));
                end
                % record Zero Fuel condition
                Xsc(ii,jj,s) = Xval;
                end
            end
            end
        end
    end
    % uncomment the following only for debug purposes    
%     fprintf('3# r max = %.2d     r min = %.2d \n',max(max(r)),min(min(r)));
%     fprintf('3# Xsc   = %.2d     Xsc   = %.2d \n',max(max(Xsc(f0xmin:f0xmax,f0ymin:f0ymax,1))),min(min(Xsc(f0xmin:f0xmax,f0ymin:f0ymax,1))));
    c(:,:,2) = Dt*rc;
    Mc = zeros(nx,ny);
    for s=1:1:5
        Mc = Mc + Xsc(:,:,s)*FireCoeff.ms(s);
    end
    cpc = zeros(nx,ny);
    for s=1:1:5
        cpc = cpc + (FireCoeff.cps(s)*FireCoeff.ms(s))*Xsc(:,:,s)./Mc;
    end
    hcc = zeros(nx,ny);
    for s=1:1:4
        hcc = hcc - FireCoeff.coef(s)*(FireCoeff.Hrefs(s)+FireCoeff.ms(s)*FireCoeff.cps(s)*((Tr+0.5*c(:,:,1))./cpc/FireCoeff.rho0-FireCoeff.Tref)).*(Mc.^(-1));
    end
    hcc = hcc*FireCoeff.Ch;
    
    % 4.4 "d" calculation 
    [d(:,:,1), hc4] = F_Tr(FireState,FireCoeff, Tr+c(:,:,1),hcc,rc,cpc,Mc);
    d(:,:,1) = Dt*d(:,:,1);
    rd   = F_Xf(FireCoeff, Tr+c(:,:,1),FireState.cp,Xsc(:,:,1),Xsc(:,:,2),nx,ny);
    Xsd = FireState.Xs;
    for jj=1:1:ny
        for ii=1:1:nx
            if FireState.Xs(ii,jj,1)~=0
            for s=1:1:2 % only for fuel and oxigen (partial for RK steps)
                if (FireState.Xs(ii,jj,1)~=FireCoeff.Xe)
                Xval = FireState.Xs(ii,jj,s) -(Dt*FireCoeff.coef(s)/FireCoeff.ms(1))*rd(ii,jj)*FireState.M(ii,jj);
                if (s==1)&&(Xval<FireCoeff.Xe); Xval=FireCoeff.Xe;%Xval=Xe; 
                    % combustion rate correction in case of extinction
                    rd(ii,jj)=(Xsd(ii,jj,s)-FireCoeff.Xe)*FireCoeff.ms(1)/(Dt*FireCoeff.coef(s)*FireState.M(ii,jj));
                end
                % record Zero Fuel condition
                Xsd(ii,jj,s) = Xval;
                end
            end
            end
        end
    end
    % uncomment the following only for debug purposes
%     fprintf('4# r max = %.2d     r min = %.2d \n',max(max(r)),min(min(r)));
%     fprintf('4# Xsd   = %.2d     Xsd   = %.2d \n',max(max(Xsd(f0xmin:f0xmax,f0ymin:f0ymax,1))),min(min(Xsd(f0xmin:f0xmax,f0ymin:f0ymax,1))));
    d(:,:,2) = Dt*rd;
    clear Xsc hcc rc cpc
    
    % 4.5 Runge Kutta 4th order approximated solution
    Trn  = Tr  + (a(:,:,1)+2*b(:,:,1)+2*c(:,:,1)+d(:,:,1))/6;
    DXfn = (a(:,:,2)+2*b(:,:,2)+2*c(:,:,2)+d(:,:,2))/6;
    r = DXfn/Dt;

    Xfn  = Xf  + DXfn.*FireState.M/FireCoeff.ms(1);
    
    % 4.6 Final update other properties     
    for jj=1:1:ny
        for ii=1:1:nx
            for s=1:1:4
                if (FireState.Xs(ii,jj,s)~=0)
                Xval = FireState.Xs(ii,jj,s) -DXfn(ii,jj)*FireCoeff.coef(s)/FireCoeff.ms(1)*FireState.M(ii,jj);
                if ((s==1)||(s==2))&&(Xval<FireCoeff.Xe); Xval=FireCoeff.Xe;%Xval=Xe; 
                    % combustion rate correction in case of extinction
                    r(ii,jj)=(FireState.Xs(ii,jj,s)-FireCoeff.Xe)*FireCoeff.ms(1)/(Dt*FireCoeff.coef(s)*FireState.M(ii,jj));
                end
                FireState.Xs(ii,jj,s) = Xval;
                if s==1
                    Xfn(ii,jj) = Xval;
                end
                end
            end
        end
    end
    
    FireState.M = zeros(nx,ny);
    for s=1:1:5
        FireState.M = FireState.M + FireState.Xs(:,:,s)*FireCoeff.ms(s);
    end
    
    FireState.cp = zeros(nx,ny);
    for s=1:1:5
        FireState.cp = FireState.cp + FireState.Xs(:,:,s)*FireCoeff.ms(s)*FireCoeff.cps(s);
    end
    FireState.cp = FireState.cp./FireState.M;
    
    FireState.hc = zeros(nx,ny);
    for s=1:1:4
        FireState.hc = FireState.hc - FireCoeff.coef(s)*(FireCoeff.Hrefs(s)+(FireCoeff.ms(s)*FireCoeff.cps(s))*(Trn./FireState.cp/FireCoeff.rho0-FireCoeff.Tref)).*(FireState.M.^(-1));
    end
    FireState.hc = FireState.hc*FireCoeff.Ch;
    
    FireState.Pmax = (hc1+2*hc2+2*hc3+hc4)/6;
    FireState.Pmax = FireState.Pmax * FireCoeff.corr*(FireCoeff.Fl)*(1/FireCoeff.Xfref)*((600*FireCoeff.rho0)/(700))/1000;
    FireState.Enrg = (FireState.hc.*(FireState.Xs(:,:,1))/FireCoeff.Xfref*FireCoeff.rho0*ThicknessGas);
    FireState.Temp = (Trn./FireState.cp/FireCoeff.rho0);
    % Calculate area affected by fire
    count = 0;
    for i=1:1:nx
        for j=1:1:ny
            if (FireState.Xs(i,j,1)<FireCoeff.Xf0(i,j))
                count = count+1;
            end
        end
    end
    FireState.Ab = count*Dx*Dy;
    FireSpy.T  = FireState.Temp(FireCoeff.posSpy(1),FireCoeff.posSpy(2))/FireCoeff.Tig;
    FireSpy.Xf = FireState.Xs(FireCoeff.posSpy(1),FireCoeff.posSpy(2),1)/FireCoeff.Xfref; 
    FireSpy.XC = FireState.Xs(FireCoeff.posSpy(1),FireCoeff.posSpy(2),3);
    storess = FireState.Xs(3:end-2,3:end-2,1);
    [jjnx,iinx] = find(storess<FireCoeff.Xfref); %correggi con Xf0
    jj_min= find(jjnx==min(min(jjnx)));
    jj_max= find(jjnx==max(max(jjnx)));
    
    FireState.Tmax = max(max(FireState.Temp));
    FireState.Temp = FireState.Temp';
    FireState.Enrg = FireState.Enrg';
    % uncomment the following only for debug purposes
%     fprintf('   T max = %.2d     T min = %.2d \n',max(max(FireState.Temp)),min(min(FireState.Temp)));
    %======================================================================
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. SUB-ROUTINES

function [value, dhcomb] = F_Tr(FireState,FireCoeff, Tr,hc,r,cp,M)  %      

nx = FireCoeff.nx;
ny = FireCoeff.ny;
Dx = FireCoeff.Dx;
Dy = FireCoeff.Dy;


%   matrix wind define as follows:
U = FireState.wind(:,:,1); % selecting matrix of x components
V = FireState.wind(:,:,2); % selecting matrix of y components

% persistent dTcomb Tdiffx Tdiffy Hdiffx Hdiffy dTRadx dTRady dTRad3 dTCa
T =(1/FireCoeff.rho0)*Tr./cp;  %initializing T matrix

% Combustion enthalpy hc*r    
    dhcomb = max(max(-M.*hc.*r/FireCoeff.ms(1)));
    dTcomb = -FireCoeff.rho0*hc.*M.*r/FireCoeff.ms(1); 
    % uncomment following fprintf only for debuggin purposes
%     fprintf('hc: %s %s\n', max(max(hc)), min(min(hc)));
%     fprintf('r: %s %s\n', max(max(r)), min(min(r)));
%     fprintf('dTcomb: %s %s\n', max(max(dTcomb./cp/FireCoeff.rho0)), min(min(dTcomb./cp/rho0)));
    value = dTcomb;
    

% % T Wind transport ( first order upwind scheme with linear damping- finite difference)
%-------------% initialisations for wind loop---------

%     dTrdx=zeros(nx,ny);
%     dTrdy=zeros(nx,ny);

     a=1; % coefficient of upwindness ( 1 for U+ and 0 for U- )[is automaticly calculated in loop]

    for i=1:nx   %here starts the  "wind_loop" / calculates heat convection due to wind
         for j =1:ny
   %transport in x-direction- first order accurate upwind scheme
   Cd=0.01;   %linear damping coefficient to reduce false diffusion due to use of first order method

    if U(i,j)>0
          a=1;
        if i==1  %|| i==2                  
           dTrdx(i,j) = a*(T(i,j)-FireCoeff.Tamb)/Dx + (1-a)*(T(i+1,j)-T(i,j))/Dx; 
        elseif i==nx %|| i==nx-1
           dTrdx(i,j) = a*(T(i,j)-T(i-1,j))/Dx + (1-a)*(FireCoeff.Tamb-T(i,j))/Dx;
        else
%             dTrdx(i,j) = a*(3*T(i,j)-4*T(i-1,j)+T(i-2,j))/(2*Dx) + (1-a)*(-T(i+2,j)+4*T(i+1,j) -3*T(i,j))/(2*Dx);%+0.02*(T(i,j+1)-2*T(i,j)+T(i,j-1));
             dTrdx(i,j) = a*(T(i,j) - T(i-1,j))/(Dx) + (1-a)*(T(i+1,j)-T(i,j))/(Dx)+(Cd*(T(i+1,j)-2*T(i,j)+T(i-1,j)));
        end
    else
          a=0;
        if i==1 %|| i==2                  
           dTrdx(i,j) = a*(T(i,j)-FireCoeff.Tamb)/Dx + (1-a)*(T(i+1,j)-T(i,j))/Dx; 
        elseif i==nx %|| i==nx-1
           dTrdx(i,j) = a*(T(i,j)-T(i-1,j))/Dx + (1-a)*(FireCoeff.Tamb-T(i,j))/Dx;
        else
%             dTrdx(i,j) = a*(3*T(i,j)-4*T(i-1,j)+T(i-2,j))/(2*Dx) + (1-a)*(-T(i+2,j)+4*T(i+1,j) -3*T(i,j))/(2*Dx);%+0.02*(T(i,j+1)-2*T(i,j)+T(i,j-1));
             dTrdx(i,j) = a*(T(i,j) - T(i-1,j))/(Dx) + (1-a)*(T(i+1,j)-T(i,j))/(Dx)+(Cd*(T(i+1,j)-2*T(i,j)+T(i-1,j)));
        end
    end 
     
     
    %transport in y-direction- first order accurate upwind scheme   
    if V(i,j)>=0;
        a=1;
        if (j==1) %|| (j==2)
             dTrdy(i,j) = a*(T(i,j)-FireCoeff.Tamb)/Dy + (1-a)*( FireCoeff.Tamb-T(i,j))/Dy;
        elseif (j==ny) %|| (j==ny-1)
             dTrdy(i,j) = a*(T(i,j)-T(i,j-1))/Dy + (1-a)*(diff([T(i,j); FireCoeff.Tamb]))/Dy;
        else
 %           dTrdy(i,j) = a*(3*T(i,j)-4*T(i,j-1)+T(i,j-2))/(2*Dy) + (1-a)*(-T(i,j+2)+4*T(i,j+1) -3*T(i,j))/(2*Dy);%+0.02*(T(i,j+1)-2*T(i,j)+T(i,j-1));
             dTrdy(i,j) = a*(T(i,j)-T(i,j-1))/(Dy) + (1-a)*(T(i,j+1) -T(i,j))/(Dy)+(Cd*(T(i,j+1)-2*T(i,j)+T(i,j-1)));
        end
    else
        a=0;
        if (j==1) %|| (j==2)
             dTrdy(i,j) = a*(T(i,j)-FireCoeff.Tamb)/Dy + (1-a)*( FireCoeff.Tamb-T(i,j))/Dy;
        elseif (j==ny) %|| (j==ny-1)
             dTrdy(i,j) = a*(T(i,j)-T(i,j-1))/Dy + (1-a)*(diff([T(i,j); FireCoeff.Tamb]))/Dy;
        else
%             dTrdy(i,j) = a*(3*T(i,j)-4*T(i,j-1)+T(i,j-2))/(2*Dy) + (1-a)*(-T(i,j+2)+4*T(i,j+1) -3*T(i,j))/(2*Dy);%+0.02*(T(i,j+1)-2*T(i,j)+T(i,j-1));
              dTrdy(i,j) = a*(T(i,j)-T(i,j-1))/(Dy) + (1-a)*(T(i,j+1) -T(i,j))/(Dy)+(Cd*(T(i,j+1)-2*T(i,j)+T(i,j-1)));
         end
    end
        
        end
    end          %%%%%%%%%%  here the "wind_loop" ends
    dTrdx   = - U.*FireCoeff.rho0.*cp.*dTrdx;
    dTrdy   = - V.*FireCoeff.rho0.*cp.*dTrdy;
%----------------------%wind loop ended here---------------------   
 

% Uncomment the following fprintf only for inspection    
%     fprintf('U= %.2e, V= %.2e, Vel= %.2e\n',U(1,1),V(1,1),Vel);
%     fprintf('sw:%i se:%i ne:%i nw:%i\n',sw,se,ne,nw);
%     fprintf('dTrdx= %.2e, dTrdy= %.2e, dTrdiag= %.2e\n',dTrdx(51,51),dTrdy(51,51),dTrdiag(51,51));
%     fprintf('\n');
%    % sw=0; se=0; ne=0; nw=0; % variables only used for inspection, not used in computation
    inx = find(dTrdx<0);
    iny = find(dTrdy<0);
    dTrdx(inx)=dTrdx(inx)*FireCoeff.FLm;
    dTrdy(iny)=dTrdy(iny)*FireCoeff.FLm;
    
    value = value + dTrdx + dTrdy; 
    clear dTrdx dTrdy dTrdiag
%%
% T Diffusion                                                  
  % d2T/Dx2
    % West (Dirichlet)
    a1 = 2*(2./(cp(1,:)+cp(2,:))+1./cp(1,:));
    a2 = 2*(2./(cp(1,:)+cp(2,:))-3./cp(1,:)+16./(cp(1,:)+FireCoeff.cp0));
    a3 = 2*(-4./cp(1,:)+16./(cp(1,:)+FireCoeff.cp0));
    Tdiffx(1,:) = FireCoeff.Kth*(a1.*Tr(2,:)/FireCoeff.rho0-a2.*Tr(1,:)/FireCoeff.rho0+a3*FireCoeff.cp0*FireCoeff.Tamb)/(3*Dx^2);
    % East (Dirichlet)
    a1 = 2*(2./(cp(nx,:)+cp(nx-1,:))+1./cp(nx,:));
    a2 = 2*(2./(cp(nx,:)+cp(nx-1,:))-3./cp(nx,:)+16./(cp(nx,:)+FireCoeff.cp0));
    a3 = 2*(-4./cp(nx,:)+16./(cp(nx,:)+FireCoeff.cp0));
    Tdiffx(nx,:) = FireCoeff.Kth*(a1.*Tr(nx-1,:)/FireCoeff.rho0-a2.*Tr(nx,:)/FireCoeff.rho0+a3*FireCoeff.cp0*FireCoeff.Tamb)/(3*Dx^2);
    % Internal
    TrR = Tr/FireCoeff.rho0;
    cpm = 1./avg(cp);
    for i=2:1:nx-1
        Tdiffx(i,:) = FireCoeff.Kth*(TrR(i+1,:).*cpm(i,:)-TrR(i,:).*(cpm(i-1,:)+cpm(i,:))+TrR(i-1,:).*cpm(i-1,:))/(Dx^2);
    end
    % uncomment the following only for debug purposes    
%     fprintf('Tdiffx: %s %s\n', max(max(Tdiffx./cp/rho0)), min(min(Tdiffx./cp/rho0)));
    value = value + Tdiffx;
    clear cpm
  % d2T/Dy2
    % South (Dirichlet)
    a1 = 2*(2./(cp(:,1)+cp(:,2))+1./cp(:,1));
    a2 = 2*(2./(cp(:,1)+cp(:,2))-3./cp(:,1)+16./(cp(:,1)+FireCoeff.cp0));
    a3 = 2*(-4./cp(:,1)+16./(cp(:,1)+FireCoeff.cp0));
    Tdiffy(:,1) = FireCoeff.Kth*(a1.*Tr(:,2)/FireCoeff.rho0-a2.*Tr(:,1)/FireCoeff.rho0+a3*FireCoeff.cp0*FireCoeff.Tamb)/(3*Dx^2);
    % North (Dirichlet)
    a1 = 2*(2./(cp(:,ny)+cp(:,ny-1))+1./cp(:,ny));
    a2 = 2*(2./(cp(:,ny)+cp(:,ny-1))-3./cp(:,ny)+16./(cp(:,ny)+FireCoeff.cp0));
    a3 = 2*(-4./cp(:,ny)+16./(cp(:,ny)+FireCoeff.cp0));
    Tdiffy(:,ny) = FireCoeff.Kth*(a1.*Tr(:,ny-1)/FireCoeff.rho0-a2.*Tr(:,ny)/FireCoeff.rho0+a3*FireCoeff.cp0*FireCoeff.Tamb)/(3*Dx^2);
    % Internal
    TrR = Tr/FireCoeff.rho0;
    cpm = 1./avg((cp)')';
    for j=2:1:ny-1
        Tdiffy(:,j) = FireCoeff.Kth*(TrR(:,j+1).*cpm(:,j)-TrR(:,j).*(cpm(:,j-1)+cpm(:,j))+TrR(:,j-1).*cpm(:,j-1))/(Dy^2);
    end
    % uncomment the following only for debug purposes    
%     fprintf('Tdiffy: %s %s\n', max(max(Tdiffy./cp/rho0)), min(min(Tdiffy./cp/rho0)));
    value = value + Tdiffy;
    clear TrR cpm
    
% % hc Diffusion     
% d2hc/Dx2  ************************************************************
    Hdiffx = zeros(nx,ny);
    TcpM = T./(cp.*M);
    TcpMavg = avg(TcpM);     
    diffx = 0;
    for s=1:5
        diffx = diffx + FireCoeff.cps(s)*FireCoeff.ms(s)*(FireState.Xs(FireCoeff.f0(1)+1:FireCoeff.f0(2),:,s)-FireState.Xs(FireCoeff.f0(1):FireCoeff.f0(2)-1,:,s));
    end
    diffx = TcpMavg(FireCoeff.f0(1):FireCoeff.f0(2)-1,:).*diffx/Dx;
    Hdiffx(FireCoeff.f0(1)+1:FireCoeff.f0(2)-1,:) = FireCoeff.Kth*(diffx(2:end,:)-diffx(1:end-1,:))/Dx;
    value = value + Hdiffx;
    clear cpm
    
  % d2hc/Dx2  ************************************************************
     Hdiffy = zeros(nx,ny);
    TcpM = T./(cp.*M);
    TcpMavg = (avg(TcpM'))';
    diffy = 0;
    for s=1:5
        diffy = diffy + FireCoeff.cps(s)*FireCoeff.ms(s)*(FireState.Xs(:,FireCoeff.f0(3)+1:FireCoeff.f0(4),s)-FireState.Xs(:,FireCoeff.f0(3):FireCoeff.f0(4)-1,s));
    end
    diffy = TcpMavg(:,FireCoeff.f0(3):FireCoeff.f0(4)-1).*diffy/Dy;
    Hdiffy(:,FireCoeff.f0(3)+1:FireCoeff.f0(4)-1) = FireCoeff.Kth*(diffy(:,2:end)-diffy(:,1:end-1))/Dy;
    value = value + Hdiffy;
    clear cpm
    
% Radiation                                                    
  % dRad/Dx
    % West (Dirichlet)
    dTRadx(1,:) = (FireCoeff.Ot*FireCoeff.sigma*FireCoeff.epsilon/(3*Dx^2))*(8*(FireCoeff.Tamb^4)+16*(FireCoeff.Tamb^3)*T(1,:)-48*FireCoeff.Tamb*(T(1,:).^3)+15*(T(1,:).^4)+6*(T(1,:).^3).*T(2,:)+2*T(1,:).*(T(2,:).^3)+T(2,:).^4);
    % East (Dirichlet)
    dTRadx(nx,:) = (FireCoeff.Ot*FireCoeff.sigma*FireCoeff.epsilon/(3*Dx^2))*(8*(FireCoeff.Tamb^4)+16*(FireCoeff.Tamb^3)*T(nx,:)-48*FireCoeff.Tamb*(T(nx,:).^3)+15*(T(nx,:).^4)+6*(T(nx,:).^3).*T(nx-1,:)+2*T(nx,:).*(T(nx-1,:).^3)+T(nx-1,:).^4);    
    % Internal
    Tavg = avg(T);
    for i=2:1:nx-1
        dTRadx(i,:) = (FireCoeff.Ot*FireCoeff.sigma*FireCoeff.epsilon*4/Dx^2)*((T(i+1,:)-T(i,:)).*Tavg(i,:).^3+(T(i-1,:)-T(i,:)).*Tavg(i-1,:).^3);
    end
%     fprintf('dTRadx: %s %s\n', max(max(dTRadx)), min(min(dTRadx)));
    value = value + dTRadx;
    clear Tavg
  % dRad/Dy
    % West (Dirichlet)
    dTRady(:,1) = (FireCoeff.Ot*FireCoeff.sigma*FireCoeff.epsilon/(3*Dy^2))*(8*(FireCoeff.Tamb^4)+16*(FireCoeff.Tamb^3)*T(:,1)-48*FireCoeff.Tamb*(T(:,1).^3)+15*(T(:,1).^4)+6*(T(:,1).^3).*T(:,2)+2*T(:,1).*(T(:,2).^3)+T(:,2).^4);
    % East (Dirichlet)
    dTRady(:,ny) = (FireCoeff.Ot*FireCoeff.sigma*FireCoeff.epsilon/(3*Dy^2))*(8*(FireCoeff.Tamb^4)+16*(FireCoeff.Tamb^3)*T(:,ny)-48*FireCoeff.Tamb*(T(:,ny).^3)+15*(T(:,ny).^4)+6*(T(:,ny).^3).*T(:,ny-1)+2*T(:,ny).*(T(:,ny-1).^3)+T(:,ny-1).^4);
    % Internal
    Tavg = avg(T')';
    for j=2:1:ny-1
        dTRady(:,j) = (FireCoeff.Ot*FireCoeff.sigma*FireCoeff.epsilon*4/Dy^2)*((T(:,j+1)-T(:,j)).*Tavg(:,j).^3+(T(:,j-1)-T(:,j)).*Tavg(:,j-1).^3);
    end
%     fprintf('dTRady: %s %s\n', max(max(dTRady)), min(min(dTRady)));
    value = value + dTRady;
    clear Tavg
    
% 3D pseudo-Radiation (z-dir.)                                 
    dTRad3 = FireCoeff.sigma*FireCoeff.epsilon*(FireCoeff.Tamb^4 - T.^4)*(1/FireCoeff.Dz);
%     fprintf('dTRad3: %s %s\n', max(max(dTRad3)), min(min(dTRad3)));
    value = value + dTRad3;
    
% 3D pseudo-Convection (z-dir.)                                
    dTCa = FireCoeff.Ca*(FireCoeff.Tamb-T);
%     fprintf('dTCa: %s %s\n', max(max(dTCa)), min(min(dTCa)));
    value = value + dTCa;
    clear T

function value = F_Xf(FireCoeff, Tr,cp,Xf,Xo2,nx,ny)%                       
    Temp = (1/FireCoeff.rho0)*Tr./cp;
    value = zeros(nx,ny);
    for j=1:1:ny
        for i=1:1:nx
            kdFl = Spr(Temp(i,j),Xf(i,j),Xo2(i,j),FireCoeff); % switch Flame ignition
            value(i,j) = -kdFl*FireCoeff.Ar*Temp(i,j)*(FireCoeff.rho0^0.3)*(Xf(i,j).^0.5)*Xo2(i,j)*exp(-FireCoeff.Ta/Temp(i,j));
        end
    end

% Average matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = avg(A,k)
    if nargin<2, k = 1; end
    if size(A,1)==1, A = A'; end
    if k<2
        B = (A(2:end,:)+A(1:end-1,:))/2; 
      else
        B = avg(A,k-1);
    end
    if size(A,2)==1, B = B'; end
 
% Ignition Condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spread = Spr(T,X,Xo,FireCoeff)
    if (T>FireCoeff.Tig)&&(X>FireCoeff.Xe)&&(Xo>FireCoeff.Xo2e) 
         spread = 1;
    else spread = 0; 
    end