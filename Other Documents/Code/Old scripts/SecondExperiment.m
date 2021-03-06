clear all;close all;clc

%% Importing Nasa tables

addpath('General/Nasa'); % Add directory of Nasa routines to Matlab-path
global Runiv 
Runiv = 8.314; %universal gas constant
DBdir = '4GB10_matlab_and_trainingdata\matlab\General\Nasa';
DBname = 'NasaThermalDatabase.mat';
load(fullfile(DBdir,DBname));
Sp
El

%% Importing data
%for now, only example data can be used.

addpath('General');
DataDir='../Data/Training Set';ColumnOrder={'time','Sensor','Encoder'};
% DataDir='../Data/Gasoline';cOrder={'time'','Encoder','Sensor};
col = lines(3);
% Loading all measurments in DataDir
figure(1)
Files=dir(fullfile(DataDir,'*.txt'));nFiles=length(Files);                  % dir gives a directory listing, only *.txt files in this case
for i=1:nFiles
    fname       = Files(i).name;                                            % Take a name from the list
    curfilename = fullfile(DataDir,fname);                                  % Create the full name
    Data        = ImportData4GB10(curfilename,ColumnOrder);                             % Read the data. Type help ImportData4GB10
    % store it for each case. Yes a struct that contains a struct and other
    % things. Why not?
    Case(i).Data     = Data;
    Case(i).filename = fname;
    Case(i).DataDir  = DataDir;
    % preamble, put data in easy to use arrays
    t      = Data.t;
    p      = Data.pulse;
    V      = Data.Volt;
    RevEnd = Data.RevEnds;
end

%% Specifying dimensions

%INSERT MEASUREMENTS HERE

%Honda GX200 engine data sheet from online
bore                    = 0.068; %[m]; diameter of the piston cylinder
crankRadius             = 0.027; %[m]; radius of crankshaft
deadVolume              = 0; %[m^3]; volume that is left when the volume in front of the piston is minimal
connectingRodLength     = 0; %[m]; length of the rod that connects the crankshaft to the piston
RPM                     = 3000; %[revolutions per minute]; how many times the crankshaft makes a 360 degree turn in one minute
angleCrank              = 0; %[rad]; angle the crankshaft makes, MAY BE A FUNCTION OF TIME
%YOU'RE DONE. STOP INSERTING AFTER HERE

%See slide 19 of the first presentation
stroke                  = 2 * crankRadius; %[m]; length that the piston covers
displacementVolume      = 0.000196; %(pi / 4) * bore * bore * stroke; %[m^3]; volume displaced by piston
compressionRatio        = 8.5; %(displacementVolume + deadVolume) / deadVolume; %[ratio]; how many times is the maximal volume larger than the smallest
pistonVelocity          = 2 * stroke * RPM / 60; %[m/s]; how fast the piston moves up and down

%See slide 20 of the first presentation
distancePistonCrank 	= crankRadius * cos(angleCrank) + sqrt( connectingRodLength * connectingRodLength - crankRadius * crankRadius * sin(angleCrank) * sin(angleCrank) ); %[m]; distance between the crankshaft and the piston
distancePistonCrankMax  = connectingRodLength + crankRadius; %[m]; maximum distance between the crankshaft and the piston
distancePistonCrankMin  = connectingRodLength + crankRadius; %[m]; minimum distance between the crankshaft and the piston
distancePistonTop       = distancePistonCrankMax - distancePistonCrank; %[m]; distance from the piston to the maximum extended position
volume                  = pi * (bore / 2) ^ 2 * distancePistonTop + deadVolume; %[m^3]; volume compressed by the piston

%% Thermodynamic model:

%central equation:
% 0 = Qdot - Wdot + QLHV * mdotfuel - cpexh * ( Texh - Tamb ) * mdotexh
Qdot                = 0; %[J]; Generated heat, ONLY HOLDS WHEN ADIABATIC

%% Part I: Chemical bond energy

%CxHyOz + (x + y / 4 - z / 2) O2 + other stuff in the air relative to the oxygen --> xCO2 + y / 2 H2O + other stuff in the air
%ALL oxygen in the air is burnt up (stoichiometric)

%ASSSUMPTION: 100% octane
x                   = 8; %[atoms] number of CARBON atoms in a fuel molecule
y                   = 18; %[atoms] number of HYDROGEN atoms in a fuel molecule
z                   = 0; %[atoms] number of OXYGEN atoms in a fuel molecule
molarmassFuel       = 114.232; %[g/mol] G/MOL, NOT KG/MOL!!! how heavy one mole of fuel is
airDensity          = 101325; %[kg/m^3] density of the air at the specified temperature and elevation without fuel


fuelConsumption     = 2.4333 * 10^(-4); %[kg/s] how many kilogrammes of fuel are needed to let the engine operate for one second
molarRatioOxygenFuel= x + y / 4 - z / 2; %[atoms] how many oxygen atoms are needed to completely combust one fuel molecule, may be a fraction.
ratioNitrogenOxygen = 78.084 / 20.946; %[ratio] mass ratio of nitrogen to oxygen in the air 
ratioCO2Oxygen      = 0.0407 / 20.946; %[ratio] mass ratio of nitrogen to oxygen in the air 
%I neglect the other noble gasses present in the air
AFstoi              = ( molarRatioOxygenFuel * 31.9988 + ratioNitrogenOxygen * molarRatioOxygenFuel * 28.0134 + ratioCO2Oxygen * molarRatioOxygenFuel * 44.01) / (molarmassFuel * 1); 
%[ratio]; Air/fuel stoichiometric mixture determined by computing the moles
%of each component in the air needed to burn the fuel stoiciometrically
%multiplied by the molar mass divided by the mass of one mole of fuel.
mair                = displacementVolume * airDensity; %[kg/cycle]
mdotfuel            = mair / AFstoi; %[kg/cycle]; inflowing mass per second of the fuel

%STILL NEEDS TO BE DONE, SEE SLIDE 43 OF FIRST PRESENTATION
%I do not know what each variable means, or how to get that from the NASA
%tables.

QLHV                = 0; %[J/kg/cycle]; chemical energy obtainable form fuel type


ChemicalBondEnergy  = QLHV * mdotfuel; %[J/cycle]; Energy stored in the fuel
%alternative formula from slide 33:
%ChemicalBondEnergy  = mfuel * cVfuel * (T3 - T2);
%Source molar masses:
%https://www.engineeringtoolbox.com/molecular-weight-gas-vapor-d_1156.html

%% Part II: Work
%Power               = 0; %[J]; Power provided by the engine 

%Wdot                = Power; %[J]; Work done by engine

%% Part III: Sensible heat loss

%USE NASA TABLE FOR cpexh
cpexh               = 0; %[J/K]; constant pressure heat capacity of the combustion products
%At 3000 RPM, the engine provides 4.4 HP which correspons to an exhaust
%temperature of roughly 205 C
Texh                = 478.15; %[K]; temperature at the exit of the combustion chambre
Tamb                = 293.15; %[K]; temperature of the rest of the environment (NOT NECESSARILLY FUEL INTAKE TEMPERATURE)

mdotexh             = mair + mdotfuel; %[kg/cycle]; mass of the air the leaves the combustion chambre (conservation of mass)
SensibleHeatLoss    = cpexh * ( Texh - Tamb ) * mdotexh; %[J/cycle]; heat lost during each cycle (due to passive heat dissipation and such)

%% Part IV: Results

%rewriting the first equation
Wdot                = 3300; %ChemicalBondEnergy - SensibleHeatLoss; %[J/cycle]; work done by the engine each cycle

%bsfc                = mdotfuel / Power * 1000 / 3600000 ; %[g/kWhr]; brake specific fuel consumption
%SPi                = mdot(OF YOUR SPECIFIC FUEL TYPE) / Power * 1000 / 3600000 %[g/kWhr];
%specific emissions
%efficiencyCarnot    = Wdot * ChemicalBondEnergy; %[ratio]; theoretically highest obtainable efficiency


%% Theoretical pV-diagram

AirPressure = 100960; %[N/m^2] atmospheric pressure

%table
%           p           V           T
%   1       p1          V1          Tamb
%   2       p2          V2          T2
%   3       p3          V2          T3
%   4       pamb        V1          Texh

%THIS MUST BE CHANGED
pamb = AirPressure;

molarmassAir = 28.9647; %[g/mol] https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html
mmRatioWater = y / 2 * 18.02 / molarmassFuel;%[ratio] the mass of the formed water with respect to the fuel
mmRatioCO2 = x * 44.01 / molarmassFuel; %[ratio] the mass of the formed carbon dioxide with respect to the fuel

cpT1 = AFstoi * CpNasa(Tamb, Sp(55) ) + (1 - AFstoi) * (0.78 * CpNasa(Tamb, Sp(48)) + 0.21 * CpNasa(Tamb, Sp(4)) + 0.01 * CpNasa(Tamb, Sp(49))); %55 corresponds to gasoline, 48 = N2, 4 = O2, 49 = Ar
cvT1 = AFstoi * CvNasa(Tamb, Sp(55) ) + (1 - AFstoi) * (0.78 * CvNasa(Tamb, Sp(48)) + 0.21 * CvNasa(Tamb, Sp(4)) + 0.01 * CvNasa(Tamb, Sp(49)));
cpTexh = AFstoi * mmRatioWater * CpNasa(Texh, Sp(6) ) + AFstoi * mmRatioCO2 * CpNasa(Texh, Sp(16) ) + (1 - AFstoi) * (0.78 * CpNasa(Texh, Sp(48)) + 0.01 * CpNasa(Texh, Sp(49))); %55 corresponds to gasoline
cvTexh = AFstoi * mmRatioWater * CvNasa(Texh, Sp(6) ) + AFstoi * mmRatioCO2 * CvNasa(Texh, Sp(16) ) + (1 - AFstoi) * (0.78 * CvNasa(Texh, Sp(48)) + 0.01 * CvNasa(Texh, Sp(49)));

%Initial


%Step41
%Isochoric release
%m * cv * (T4 - T1) = Qc
%p / T = constant
%p1 = p4 * T1 / T4
p1 = pamb * Tamb / Texh;
Qc = cvT1 * ((mdotfuel + mair) / RPM / 0.0167) * (Texh - Tamb);
Qh = Qc - Wdot;
%ideal gas law
V1 = Tamb * 8.314 * ((mdotfuel * 1000 / molarmassFuel + mair * 1000 / molarmassAir) / RPM / 0.0167) / p1;
%it took me two hours to realise this...
V2 = V1 - displacementVolume;
p2 = p1 * V1^(cpT1/cvT1) / V2^(cpT1/cvT1);
T2 = p2 * V2 / 8.314 / ((mdotfuel * 1000 / molarmassFuel + mair * 1000 / molarmassAir) / RPM / 0.0167);

%Step23
%Isochoric combustion
%m * cv * (T2 - T3) = QLHV * mfuel = ChemicalBondEnergy
%(T2 - T3) = ChemicalBondEnergy / m / cv
cvT2 = AFstoi * CvNasa(T2, Sp(55) ) + (1 - AFstoi) * (0.78 * CvNasa(T2, Sp(48)) + 0.21 * CvNasa(T2, Sp(4)) + 0.01 * CvNasa(T2, Sp(49)));
T3 = T2 - ChemicalBondEnergy / ((mdotfuel + mair) / RPM / 0.0167) / cvT2;
%p / T = constant
%p2 / T2 = p3 / T3
p3 = p2 * T3 / T2;


%Step41
%Isochoric release
%m * cv * (T4 - T1) = Qout
%p / T = constant
%p1 = p4 * T1 / T4
p1 = pamb * Tamb / Texh;

%% Efficiency
% Calculate the ideal thermal efficiency (Ott efficiency)
gamma = 1 %Since gamma (Cp/Cv) can be calculated using this matlab model I gave it a value of 1. SHOULD BE CHANGED
eta_O = 1-(1/compressionRatio)^(mean(gamma)-1); 


% Calculate efficiency by dividing the work W by the input heat Q
% Work is equal to the area under the curve (trapz command)

p=1 %p isn't calculated yet.
eta = trapz(displacementVolume,p)/(QLHV*mdotfuel); 

%% plots

T12 = linspace(273, 323, 10); %to be replaced with the temperature in step12
T34 = linspace(323, 373, 10); %to be replaced with the temperature in step34
Vv1 = linspace(150, 200, 10); 
[T23,V1m] = ndgrid(T12,Vv1); %to be replaced with the temperature in step23
[T41,V1m] = ndgrid(T34,Vv1); %to be replaced with the temperature in step41
R = 0.2871;
thank_capacity = 0.0031 %m??;
desity_gasoline = 800 %kg/m??;
total_mass_fuel = thank_capacity * desity_gasoline; %kg
m = molarmassFuel  
% p1 = 7.8378;
% V1 = 50;
p = @(T,V) m*R*T./V;
figure
surf(T23, V1m, p(T23,V1m))
hold on
surf(T41, V1m, p(T41,V1m))
hold off
grid on
xlabel('T')
ylabel('V')
zlabel('p')
view(90,0)
figure
plot(Vv1, p(T12,Vv1), 'LineWidth',2)
hold on
plot(Vv1, p(T34,Vv1), 'LineWidth',2)
plot(Vv1(1)*[1 1], p([T12(1) T12(end)],[1 1]*Vv1(1)), '-g', 'LineWidth',2)
plot(Vv1(end)*[1 1], p([T34(1) T34(end)],[1 1]*Vv1(end)), '-g', 'LineWidth',2)
hold off
grid
xlabel('V')
ylabel('p')
axis([125  225    2.0  3.3])
legend('T12(273 - 323)', 'T34(323 - 373)')
