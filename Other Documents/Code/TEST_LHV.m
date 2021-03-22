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

E           = 10;                                                           % E number of fuel [%]
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
%Mi(1)      + 3Mi(6) 3.76(Mi(5))       3(Mi(3))     2(Mi(4))  3*3.76Mi(5)
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

%%%=============================================
%
% Qlhv = m(exhaust)/m(fuel) * Cp(exhaust) * (T(exhaust) - T(ambient))
%
%%%=============================================
Texh = 478.15                                                               % Exhaust temperatuer from literature

for i = 1:NElements
    hf(i)   = HNasa(Tref, Elements(i));                                     % Enthalpies of formation    
    Cp(i)   = CpNasa(Texh, Elements(i));                                    % Specific heat capacity at Texh
    hexh(i) = HNasa(Texh, Elements(i));                                     % Enthalpies at Texh
end

hf_GAS_in   = hf*YGAS_in';                                                  % Enthalpy of formation of initial compisition Gasoline 
hf_GAS_out  = hf*YGAS_out';                                                 % Enthalpy of formation of exhaust compisition Gasoline            
hf_ETH_in   = hf*YETH_in';                                                  % Enthalpy of formation of initial compisition Ethanol
hf_ETH_out  = hf*YETH_out';                                                 % Enthalpy of formation of exhaust compisition Gasoline        

hc_GAS = (-hf_GAS_out + hf_GAS_in)/YGAS_in(2)                               % Enthalpy of combustion of Gasoline
hc_ETH = (-hf_ETH_out + hf_ETH_in)/YETH_in(1)                               % Enthalpy of combustion of Ethanol

hc_MIX = [hc_ETH, hc_GAS, 0, 0, 0, 0];                                      
Qlhv_mix = Yfuel*hc_MIX'                                                    % Enthalpy of combustion of mixture - Lower Heating Value of mixture

Cp_GAS = Cp*YGAS_out';
Cp_ETH = Cp*YETH_out'; 
Cp_mix = [Cp_ETH, Cp_GAS, 0,0,0,0]*Yfuel';

loss = mMixIn*Cp_mix*(Texh-Tref); 
Qcomb = Qlhv_mix*mFuel -loss

Qlhv_GAS = mMixIn/mGAS*Cp_GAS*(Texh-Tref) ;
Qlhv_ETH = mMixIn/mETH*Cp_ETH*(Texh-Tref) ;


%% Combustion Cycle
Ca(1) = 180;                                                                % Initial crank angle [degree]
NCa=360;                                                                    % Number of crank-angles 
dCa=0.5;                                                                    % Stepsize [degree]
NSteps=NCa/dCa;
m(1) = mMixIn;

for i=2:NSteps,
Ca(i)=Ca(i-1)+dCa;
V(i)= Vcyl(Ca(i), Vc, Vd) ;                                                 % Volume for current crank-angle [m^3]
m(i)=m(i-1);                                                                % Mass is constant, valves are closed
dV(i)=V(i)-V(i-1);                                                          % Volume change
dQcom = YourModel(Ca(i));                                                   % Heat Release by combustion
dT=(dQcom-p(i-1)*dV(i))/Cv/m(i-1);                                          % 1st Law dU=dQ-pdV (closed system)
% adiabatic closed system with constant
% gas composition and constant Cv
T(i)=T(i-1)+dT;
p(i)=m(i)*Rg*T(i)/V(i); % Gaslaw
end;

%% Functions
function V = Vcyl(Ca, Vc, Vd)
    % V         - Volume at give crank angle            - [m^3]
    % Ca        - Crank angle                           - [degree]
    % Vc        - Clearance volume                      - [m^3]
    % Vd        - Displaced volume                      - [m^3]
    phi = 0;
    V=-Vd/2*cos(Ca*(2*pi/360))+Vc+Vd/2;
end


