%% Setup
clear all; close all; clc; 

%% Properties
V = 0.000195; %Intake volume [m^3]
R = 8.314; %Gas Constant [J/K mol]
BaPas = 100000; %Bar to Pascal

%% Variables
%%Due to the errors while calculating the average plots, the value for
%%P0_E0 half load is set by default to 0.40 (even if the value is not
%%right)
%P0_E15 full load is set by default to 0.55 (even if the value is not
%%right)
%Intake pressure [Bar] [no half full]
%For corectness everything was multiplied by 10^4
P0_E0 = [2000 4000 5600];

P0_E5 = [2400 4000 5500];

P0_E10 = [2600 4500 5600];

P0_E15 = [2200 4200 5500];

%n, number of moles
n_E0 = 8.5214;
n_E5 = 8.9517;
n_E10 = 9.3819;
n_E15 = 9.8121;

%% Intake temperature [K]
T_E0 = (P0_E0*BaPas*V)/(n_E0*R)

T_E05 = (P0_E5*BaPas*V)/(n_E5*R)

T_E10 = (P0_E10*BaPas*V)/(n_E10*R)

T_E15 = (P0_E15*BaPas*V)/(n_E15*R)