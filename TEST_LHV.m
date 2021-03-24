%% Test data base
clear all;
% close all;
%% This global is needed by the Cp, Cv etc functions
addpath('matlab/General/Nasa'); % Add directory of Nasa routines to Matlab-path
global Runiv Pref Tref
Runiv = 8.314;
Pref=1.01235e5; % Reference pressure, 1 atm!
Tref=298.15;    % Reference Temperature
%% Some convenient units and 
kJ=1e3;                                             
kmol=1e3;
dm=0.1;
bara=1e5;
kPa = 1000;
kN=1000;
mm  = 0.001;    %[m]
cm3 = 10^(-6);  %[m^2]

%% 
DBdir = 'General\Nasa';
DBname = 'NasaThermalDatabase';
load(fullfile(DBdir,DBname));
%% Constants
iElements   = myfind({Sp.Name},{'C2H5OH','Gasoline','H2O', 'CO2', 'N2', 'O2'});
Elements    = Sp(iElements);
NElements   = length(Elements);
Mi          = [Elements.Mass];

Xair        = [0 0 0 0 0.79 0.21];                                          % Molar compisition of air
Mair        = Xair*Mi';                                                     % Molar mass of air
Yair        = Xair.*Mi/Mair;                                                % Mass compisition of air

Vt          = 196*cm3;                                                      % Total volume cylinder [m^3]
r           = 8.5;                                                          % Compression ratio [-]
Vc          = Vt/r;                                                         % Clearance volume [m^3]
Vd          = Vt-Vc;                                                        % Displacement volume [m^3]

p1          = 0.5*bara;                                                     % Take one from coresponding list [Pa]
T1          = Tref;                                                         % Starting Temperature [K]
V1          = Vcyl(180, Vc, Vd)                                             % Starting Volume [m^3] 

E           = 0;                                                           % E number of fuel [%]
rhoC8H18    = 0.7;                                                          % Density of Gasoline [g/cm^3] 
rhoC2H5OH   = 0.79;                                                         % Density of Ethanol [g/cm^3] 
rhoFuel     = [rhoC2H5OH, rhoC8H18, 0, 0, 0, 0];
Zfuel       = [(E/100), (1-E/100), 0, 0, 0, 0];                             % Volume compisition 
Yfuel       = Zfuel.*rhoFuel/(sum(Zfuel.*rhoFuel));                         % Mass compisition
nC8H18_O2   = Elements(2).Elcomp(3)+0.25*Elements(2).Elcomp(2) - 0.5*Elements(2).Elcomp(1);  %Number of moles of O2 required for Gasoline
nC2H5OH_O2  = Elements(1).Elcomp(3)+0.25*Elements(1).Elcomp(2) - 0.5*Elements(1).Elcomp(1);  %Number of moles of O2 required for Ethanol   
mair        = ((nC8H18_O2*Mi(6)*Yfuel(2)/Mi(2)) +(nC2H5OH_O2*Mi(6)*Yfuel(1))/Mi(1))/Yair(6);
AF          = mair/1;                                                       % Stoichemetric (mass of) air to fuel ratio    

YmixIn      = (Yair*AF + Yfuel)/(AF+1);                                     % Mass compistion of initial mixture
RgMixIn     = Runiv/((YmixIn./Mi)/(sum(YmixIn./Mi))*Mi');                   % Specific gas constant for initial mixture            
mMixIn      = p1*V1/(T1*RgMixIn);                                           % Mass of air in cylinder [kg]                                           
mFuel       = mMixIn/(1+AF);                                                % Mass of fuel in cylinder [kg]
mGAS        = Yfuel(2)*mFuel;                                               % Mass of Gasoline in cylinder [kg]        
mETH        = Yfuel(1)*mFuel;                                               % Mass of Ethanol in cylinder [kg]
mAir        = mFuel*AF;                                                     % Mass of air in cylinder [kg]

%% Combustion computations
%%%=============================================
%
%GASOLINE - COMPLETE COMBUSTION
% C7.76H13.1    +11.035(O2)   +   3.76(N2)   -->7.76CO2   +13.1/2*H2O  + N2
% Mi(2)          11.035(Mi(6)     3.76(Mi(5))  7.76(Mi(4) 13.1/2(Mi(3))  N2
%
%%%=============================================

mGAS_fuel_in    = Mi(2)/Mi(2);                                              % 1 kg of fuel - Gasoline
mGAS_O2_in      = nC8H18_O2*Mi(6)/Mi(2);                                    % Mass of O2 when 1 kg of fuel is burned
mGAS_N2_in      = nC8H18_O2*3.76*Mi(5)/Mi(2);                               % Mass of N2 when 1 kg of fuel is burned
mGAS_rate_in    = mGAS_fuel_in + mGAS_O2_in + mGAS_N2_in;                   % Total mass of mixture coming in
YGAS_in         = [0, mGAS_fuel_in, 0, 0, mGAS_N2_in, mGAS_O2_in]./mGAS_rate_in;

mGAS_fuel_out   = 0;                                                        % Complete Combustion no fuel left
mGAS_CO2_out    = 7.76*Mi(4)/Mi(2);                                         % Mass of CO2 when 1 kg of fuel is burned
mGAS_H20_out    = 13.1/2*Mi(3)/Mi(2);                                       % Mass of H2O when 1 kg of fuel is burned
mGAS_N2_out     = mGAS_N2_in;                                               % Mass of N2 remains constant, no reaction during combustion
mGAS_rate_out   = mGAS_fuel_out + mGAS_CO2_out + mGAS_H20_out +mGAS_N2_in;
YGAS_out        = [0, 0, mGAS_H20_out, mGAS_CO2_out, mGAS_N2_out, 0]./mGAS_rate_out;


%%%=============================================
%
%ETHANOL - COMPLETE COMBUSTION
%C2H5OH     + 3(O2 + 3.76N2)   -->     3H2O    +    2CO2    + 3*3.76N2
%Mi(1)      + 3(Mi(6) 3.76(Mi(5)))       3(Mi(3))     2(Mi(4))  3*3.76Mi(5)
%
%%%=============================================

mETH_fuel_in    = Mi(1)/Mi(1);                                              % 1 kg of fuel - ETHANOL
mETH_O2_in      = nC2H5OH_O2*Mi(6)/Mi(1);                                   % Mass of O2 when 1 kg of fuel is burned
mETH_N2_in      = nC2H5OH_O2*3.76*Mi(5)/Mi(1);                              % Mass of N2 when 1 kg of fuel is burned
mETH_rate_in    = mETH_fuel_in + mETH_O2_in + mETH_N2_in;                   % Total mass of mixture coming in
YETH_in         = [mETH_fuel_in, 0, 0, 0, mETH_N2_in, mETH_O2_in]./mETH_rate_in;

mETH_fuel_out   = 0;                                                        % Complete Combustion no fuel left
mETH_CO2_out    = 2*Mi(4)/Mi(1);                                            % MASS of CO2 when 1 kg of fuel is combusted
mETH_H20_out    = 3*Mi(3)/Mi(1);                                            % Mass of H2O when 1 kg of fuel is combusted
mETH_N2_out     = mETH_N2_in;                                               % MASS of N2 Remains constant
mETH_rate_out   = mETH_fuel_out + mETH_CO2_out + mETH_H20_out +mETH_N2_in;  
YETH_out        = [0, 0, mETH_H20_out, mETH_CO2_out, mETH_N2_out, 0]./mETH_rate_out;

YMixOut = (1-E/100)*YGAS_out + (E/100)*YETH_out;

%%%=============================================
%
% Qlhv = m(exhaust)/m(fuel) * Cp(exhaust) * (T(exhaust) - T(ambient))
%
%%%=============================================
Texh = 478.15                                                               % Exhaust temperatuer from literature

for i = 1:NElements
    hf(i)   = HNasa(Tref, Elements(i));                                     % Enthalpies of formation    
%     Cp_mix(i)   = CpNasa(Texh, Elements(i));                                    % Specific heat capacity at Texh
%     Cv_mix(i)   = Cp
    hexh(i) = HNasa(Texh, Elements(i));                                     % Enthalpies at Texh
end

hf_GAS_in   = hf*YGAS_in';                                                  % Enthalpy of formation of initial compisition Gasoline 
hf_GAS_out  = hf*YGAS_out';                                                 % Enthalpy of formation of exhaust compisition Gasoline            
hf_ETH_in   = hf*YETH_in';                                                  % Enthalpy of formation of initial compisition Ethanol
hf_ETH_out  = hf*YETH_out';                                                 % Enthalpy of formation of exhaust compisition Gasoline        

hc_GAS = (-hf_GAS_out + hf_GAS_in)/YGAS_in(2);                               % Enthalpy of combustion of Gasoline
hc_ETH = (-hf_ETH_out + hf_ETH_in)/YETH_in(1);                               % Enthalpy of combustion of Ethanol

hc_MIX = [hc_ETH, hc_GAS, 0, 0, 0, 0];                                      
Qlhv_mix = Yfuel*hc_MIX'                                                    % Enthalpy of combustion of mixture - Lower Heating Value of mixture

% Cp_GAS = Cp_mix*YGAS_out';
% Cp_ETH = Cp_mix*YETH_out'; 
% Cp_mix = [Cp_ETH, Cp_GAS, 0,0,0,0]*Yfuel';
% 
% loss = mMixIn*Cp_mix*(Texh-Tref); 
% Qcomb = Qlhv_mix*mFuel -loss
% 
% Qlhv_GAS = mMixIn/mGAS*Cp_GAS*(Texh-Tref) ;
% Qlhv_ETH = mMixIn/mETH*Cp_ETH*(Texh-Tref) ;


%% Combustion Cycle
Ca(1)       = 180;                                                                % Initial crank angle [degree]
NCa         =360;                                                                    % Number of crank-angles 
dCa         =0.5;                                                                    % Stepsize [degree]
NSteps      =NCa/dCa;
m(1)        = mMixIn;
p(1)        = p1;
T(1)        = T1;
V(1)        = V1;
RgMixOut    = Runiv/((YMixOut./Mi)/(sum(YMixOut./Mi))*Mi');                 % Specific gas constant for initial mixture
Rg(1)       = RgMixIn
pr = 0;
Tr = 0;
Vr = 0;
x = 0;
for i=2:NSteps
    Ca(i)   =Ca(i-1)+dCa;
    V(i)    = Vcyl(Ca(i), Vc, Vd) ;                                                 % Volume for current crank-angle [m^3]
    m(i)=m(i-1);                                                                % Mass is constant, valves are closed
    dV(i)=V(i)-V(i-1);                                                          % Volume change
    if Ca(i) == 345
        pr = p(i-1);
        Tr = T(i-1);
        Vr = V(i-1);
    end
    
    [Qcom(i), dXb(i)] = HeatReleased(Ca(i), Qlhv_mix, mFuel, dCa); 
    x(i) = x(i-1) + dXb(i);
    for ii = 1:length(Elements)
       Cv(ii) = CvNasa(T(i-1), Elements(ii));
       Cp(ii) = CpNasa(T(i-1), Elements(ii));
    end
    Cv_In   = Cv([1:6])*YmixIn';
    Cv_Out  = Cv([1:6])*YMixOut';
    Cv_mix(i)   = (1-x(i))*Cv_In + x(i)*Cv_Out;
    
    Cp_In   = Cp([1:6])*YmixIn';
    Cp_Out  = Cp([1:6])*YMixOut';
    Cp_mix(i)   = (1-x(i))*Cp_In + x(i)*Cp_Out;
    
    Rg(i)   = Cp_mix(i) - Cv_mix(i);
    k(i) = Cp_mix(i)/Cv_mix(i);
    
    [Qloss(i)] = HeatLoss(Ca(i), T(i-1), p(i-1), pr, Tr, V(i-1), k(i), p1, V1, Vr)*(0.04/(720))*dCa;
    [~, hc(i)] = HeatLoss(Ca(i), T(i-1), p(i-1), pr, Tr, V(i-1), k(i), p1, V1, Vr);
    Q(i) = Qloss(i) + Qcom(i);
    dT(i)=(Q(i)-p(i-1)*dV(i))/Cv_mix(i)/m(i);                                          % 1st Law dU=dQ-pdV (closed system)
    % adiabatic closed system with constant
    % gas composition and constant Cv
    T(i)=T(i-1)+dT(i);
    p(i)=m(i)*Rg(i)*T(i)/V(i); % Gaslaw
end;

%% Efficiency

efficiency = abs(trapz(p, V)/(Qlhv_mix*mFuel));

%% Functions
function V = Vcyl(Ca, Vc, Vd)
    % V         - Volume at give crank angle            - [m^3]
    % Ca        - Crank angle                           - [degree]
    % Vc        - Clearance volume                      - [m^3]
    % Vd        - Displaced volume                      - [m^3]
    
    V=-Vd/2*cos(Ca*(2*pi/360))+Vc+Vd/2;
end

function [Qcomb, dXb] = HeatReleased(Ca, qlhv, mfuel, dCa)
    % Qcomb     - Energy released at Crank angle        - [J]
    % dXb       - To determine the combustion progress  - [-]
    % Ca        - Crank angle                           - [degree]
    % qlhv      - Lower heating value of mixture        - [J/kg]
    % mfuel     - Mass of fuel                          - [kg]
    % dCa       - step angle                            - [degree]
    % n, a      - Constants                             - [-]
    % Theta_d   - Duration of combustion                - [degree]
    % Theta_s   - Start of combustion                   - [degree]
    
    n =3;
    a = 5;
    Theta_d = 35;
    Theta_s = (360-10);

    if Ca >= Theta_s
        Xb = 1 - exp(-a*((Ca-Theta_s)/Theta_d)^n);
        dXb = n*a*(1-Xb)/Theta_d*((Ca-Theta_s)/Theta_d)^(n-1)*dCa;
        dQcomb_dTheta = qlhv*mfuel*dXb;
        Qcomb = dQcomb_dTheta*dCa;
    elseif Ca < Theta_s  
        Qcomb =0;
        Xb = 0;
        dXb = 0;

end
end

function [Qloss, hc] = HeatLoss(Ca, T, p, pr, Tr, V, k, p1, V1, Vr)
%% variables
% Ca    - Crank angle [deg]
% Sp    - Average Piston Speed [m/s]
% Tr    - temp. start of combustion [K]
% T     - Instantaneous temp. [K]
% Tw    - Temperature of the wall [K]
% pr    - pressure start of combustion [Pa]
% p1    - Pressure start of adiabatic compression [Pa]
% V1    - Volume start of adiabatic compression [m^3]
% Vr    - Volume start of combustion [m^3]
% pm    - motored cylinder pressure at the same crank angle as 'p' [Pa]
% p     - Instantaneous pressure [pa] 
% S     - Lenght of the stroke [m]
% B     - Bore diameter [m]
% RpS   - Rotations per second [deg/s]
% k     - Heat capacity ratio (Cp/Cv) [-]
% C1,C2 - Constants for in otto cycle
% omega - The flame speed [m/s]
% hc    - Heat convection coefficient [w/m^2/K]
% h     - Height of the cylinder exposed [m}
% A     - Area of the cylinder exposed [m^2]
% Qloss - Power lost to the wall [W, J/s]
% pm    - Pressure of motored engine [Pa]
if pr == 0
    Qloss =0;
    hc = 0;
else
    if Ca >= 345 && Ca <= 540
        C1 = 2.28;
        C2 = 3.24e-03;
    elseif Ca > 180 && Ca <345
        C1 = 2.28; 
        C2 =0;
    else
        C1 = 0;
        C2 = 0;
    end

    S   = 0.055;        %[m]
    RpS = 50;
    Sp  = 2*S*RpS;      %[m/s]
    B   = 0.067;        %[m]
    Tw  = 273.15+95;    %%% CHECK %%%
    %% Computations

    Vt  = 196*10^(-6);
    Vc  = Vt/8.5;
    Vd  = Vt-Vc;

    pm = p1*(V1/V)^k;
    omega   = C1*Sp + C2*(Vd*Tr)/(pr*Vr)*(p - pm);
    hc      = 3.26*B^(-0.2)*(p/1000)^(0.8)*T^(-0.55)*omega^(0.8);           
    h       = (V - Vc)/(Vd)*S;
    A       = pi/2*B^2 + 2*pi*B*h;
    if T < Tw
        Qloss = 0;
    else 
        Qloss   = -hc*A*(T-Tw);
    end
    end
end
