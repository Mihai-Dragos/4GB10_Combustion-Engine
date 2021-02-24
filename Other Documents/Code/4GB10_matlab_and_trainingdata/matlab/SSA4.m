clear all;close all;clc
%% add general to matlab path 
addpath('C:\Users\20192303\Desktop\matlab\General\');
%%
DataDir='C:\Users\20192303\Desktop\matlab\gasoline';ColumnOrder={'time','Sensor','Encoder'};
% DataDir='../Data/Gasoline';cOrder={'time'','Encoder','Sensor};
col = lines(3);
%% Loading all measurments in DataDir
figure(1)
Files=dir(fullfile(DataDir,'*.txt'));nFiles=length(Files);                  % dir gives a directory listing, only *.txt files in this case
for i=1:nFiles
    fname       = Files(i).name;                                            % Take a name from the list
    curfilename = fullfile(DataDir,fname);                                  % Create the full name
    Data        = ImportData4GB10(curfilename,ColumnOrder);                 % Read the data. Type help ImportData4GB10
    % Store it for each case. Yes a struct that contains a struct and other
    % things. Why not?
    Case(i).Data     = Data;
    Case(i).filename = fname;
    Case(i).DataDir  = DataDir;
    % preamble, put data in easy to use arrays
    t      = Data.t;
    Pressure      = Data.pulse;
    V      = Data.Volt;
    RevEnd = Data.RevEnds;
    % Actual Plotting use of subplot
end
%% All measurements are loaded

%figure(1)                                                                 %This figure plots the pulse sensor.
%plot(t_2, V_2, 'r')
%xlim([0.01800, 0.09800])
%line(t(RevEnd_1-35),p_2(RevEnd_1-30),'LineStyle','none','Marker','s','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','y');
%xlabel('t [s]');
%ylabel('pulse [V]');
%title("Pulse sensor, 0% ethanol, full load, ")
  
%figure(2)                                                                 %This figure plots the pressure sensor
%plot(t_2, V_2, 'b')
%xlim([0, 0.3])
%line(t(RevEnd_1),p_2(RevEnd_1),'LineStyle','none','Marker','s','MarkerSize',5,'MarkerEdgeColor','r','MarkerFaceColor','y');
%xlabel('t [s]');ylabel('p Sensor Signal [V]');
%title("Pressure sensor, 0% ethanol, no load, ")

%%


 
%%
%This are some measurements from the engine.
S=0.0190;
rc=9.5;
l=0.088;
B=0.065;
%%
full_table = dlmread("C:\Users\20192303\Desktop\matlab\gasoline\1.txt");
full_time = full_table(:,1);
full_voltage = full_table(:,2);
full_pulse = full_table(:,3);


V_s = 5; %V, constant

V_full = full_voltage * V_s;
%2nd voltage 
%3rd pulse

%1 complete cycle is 0.04s = 1/25

P_full = ((V_full/V_s-0.115)/0.0385); 





%% Volume and Pv's
%This for loop determines the volume. 
for i=3:length(Case(1).Data.RevEnds)-2
  
    s_begin(i) = Case(1).Data.RevEnds(i-2);
    s_end(i) = Case(1).Data.RevEnds(i);
    revs=((3000/60)*(max(Case(1).Data.t)))
    arevs=revs*2;
    graden=arevs*360;
    Case(1).Data.cas_graden=(Case(1).Data.t-(s_begin(3)/4777))*graden
    Case(1).Data.volume=  (l+0.5*S -0.5*S*cosd(Case(1).Data.cas_graden)- sqrt(l^2-(0.25*S^2*sind(Case(1).Data.cas_graden).^2)) +(S/(rc-1)))*0.25*pi*B^2; 
   
%Plotting the Pv diagrams (intervals are determined by hand)    
    figure(7)
    hold on
    plot(Case(1).Data.volume(350:750), P_full(350:750))
    i=i+2
    
    grid on
title("PV-diagram Half load 10 percent")
ylim([0, 50])
xlabel("Volume [m^3]")
ylabel("Pressure [Pa]")
 
end
%% Loglog
%The Pv diagram in loglog scale
figure()
loglog(Case(1).Data.volume(350:750), P_full(350:750))
ylim([0.8,23])
grid on
title("PV-diagram Half load 10 percent, loglog scale")
xlabel("Volume [m^3]")
ylabel("Pressure [Pa]")

%figure(8)
%loglog(Case(a).Data.volume(1800:2000), pres(1800:2000))

%% Work 
%The work of the Pv diagram
Work =trapz(Case(1).Data.volume(350:750), P_full(350:750))

Work_real = Work* 10^5

