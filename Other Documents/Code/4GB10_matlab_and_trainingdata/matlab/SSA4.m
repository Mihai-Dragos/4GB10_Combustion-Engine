clear all;close all;clc

%%
%Importing Nasa tables

addpath('C:\Users\20192303\Documents\Year 2\Quartile 3\Combustion engine\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\matlab\General\Nasa'); % Add directory of Nasa routines to Matlab-path
global Runiv 
Runiv = 8.314; %universal gas constant
DBdir = 'General\Nasa';
DBname = 'NasaThermalDatabase';
load(fullfile(DBdir,DBname));
Sp
El

%%
%Importing data
%for now, only example data can be used.

addpath('General');
DataDir='C:\Users\20192303\Documents\Year 2\Quartile 3\Combustion engine\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set';ColumnOrder={'time','Sensor','Encoder'};
% DataDir='../Data/Gasoline';cOrder={'time'','Encoder','Sensor};
col = lines(3);
% Loading all measurments in DataDir
figure(1)
Files=dir(fullfile(DataDir,'C:\Users\20192303\Documents\Year 2\Quartile 3\Combustion engine\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\*.txt'));nFiles=length(Files);                  % dir gives a directory listing, only *.txt files in this case
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

%%
%Specifying dimensions

%INSERT MEASUREMENTS HERE
bore                    = 0; %[m]; diameter of the piston cylinder
crankRadius             = 0; %[m]; radius of crankshaft
deadVolume              = 0; %[m^3]; volume that is left when the volume in front of the piston is minimal
connectingRodLength     = 0; %[m]; length of the rod that connects the crankshaft to the piston
RPM                     = 0; %[revolutions per minute]; how many times the crankshaft makes a 360 degree turn in one minute
angleCrank              = 0; %[rad]; angle the crankshaft makes, MAY BE A FUNCTION OF TIME
%YOU'RE DONE. STOP INSERTING AFTER HERE

%See slide 19 of the first presentation
stroke                  = 2 * crankRadius; %[m]; length that the piston covers
displacementVolume      = (pi / 4) * bore * bore * stroke; %[m^3]; volume displaced by piston
compressionRatio        = (displacementVolume + deadVolume) / deadVolume; %[ratio]; how many times is the maximal volume larger than the smallest
pistonVelocity          = 2 * stroke * RPM / 60; %[m/s]; how fast the piston moves up and down

%See slide 20 of the first presentation
distancePistonCrank 	= crankRadius * cos(angleCrank) + sqrt( connectingRodLength * connectingRodLength - crankRadius * crankRadius * sin(angleCrank) * sin(angleCrank) ); %[m]; distance between the crankshaft and the piston
distancePistonCrankMax  = connectingRodLength + crankRadius; %[m]; maximum distance between the crankshaft and the piston
distancePistonCrankMin  = connectingRodLength + crankRadius; %[m]; minimum distance between the crankshaft and the piston
distancePistonTop       = distancePistonCrankMax - distancePistonCrank; %[m]; distance from the piston to the maximum extended position
volume                  = pi * (bore / 2) ^ 2 * distancePistonTop + deadVolume; %[m^3]; volume compressed by the piston
%%


S = 0.055;
r= 0.035;
l= 0.088;%moet nog gemeten worden%
rc=8.5;
%%
for i=3:length(Case(1).Data.RevEnds)-2
    
    s_begin = RevEnd_1(i);%startpunt revolution%
    s_end = RevEnd_1(i);
   revs=((3000/60)*(max(Case(1).Data.t)))
   arevs=revs*2
   graden=arevs*360
   cas_graden=(full_time-(s_begin/4777))*graden
   cas_volume = (1+0.5*S - 0.5*S*cosd(cas_graden)-sqrt(l^2-(0.25*S^2*sind(cas_graden).^2))+(S/(rc-1)))*0.25*pi*0.067;
   
   figure(9)
   hold on
   loglog(cas_volume(592:791), P_full(162:361))
   grid on
   i=i+2
end




%%
%Thermodynamic model:

%central equation:
% 0 = Qdot - Wdot + QLHV * mdotfuel - cpexh * ( Texh - Tamb ) * mdotexh
Qdot                = 0; %[J]; Generated heat, ONLY HOLDS WHEN ADIABATIC

%%
%Part I: Chemical bond energy

%CxHyOz + (x + y / 4 - z / 2) O2 + other stuff in the air relative to the oxygen --> xCO2 + y / 2 H2O + other stuff in the air
%ALL oxygen in the air is burnt up (stoichiometric)
x                   = 0; %[atoms] number of CARBON atoms in a fuel molecule
y                   = 0; %[atoms] number of HYDROGEN atoms in a fuel molecule
z                   = 0; %[atoms] number of OXYGEN atoms in a fuel molecule
molarmassFuel       = 0; %[g/mol] G/MOL, NOT KG/MOL!!! how heavy one mole of fuel is
airDensity          = 101325; %[kg/m^3] density of the air at the specified temperature and elevation without fuel


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

%%
%Part II: Work
%Power               = 0; %[J]; Power provided by the engine 

%Wdot                = Power; %[J]; Work done by engine

%%
%Part III: Sensible heat loss

%USE NASA TABLE FOR cpexh
cpexh               = 0; %[J/K]; constant pressure heat capacity of the combustion products
Texh                = 0; %[K]; temperature at the exit of the combustion chambre
Tamb                = 0; %[K]; temperature of the rest of the environment (NOT NECESSARILLY FUEL INTAKE TEMPERATURE)

mdotexh             = mair + mdotFuel; %[kg/cycle]; mass of the air the leaves the combustion chambre (conservation of mass)
SensibleHeatLoss    = cpexh * ( Texh - Tamb ) * mdotexh; %[J/cycle]; heat lost during each cycle (due to passive heat dissipation and such)

%%
%Part IV: Results

%rewriting the first equation
Wdot                = ChemicalBondEnergy - Sensible; %[J/cycle]; work done by the engine each cycle

bsfc                = mdotfuel / Power * 1000 / 3600000 ; %[g/kWhr]; brake specific fuel consumption
%SPi                = mdot(OF YOUR SPECIFIC FUEL TYPE) / Power * 1000 / 3600000 %[g/kWhr];
%specific emissions
efficiencyCarnot    = Wdot * ChemicalBondEnergy; %[ratio]; theoretically highest obtainable efficiency