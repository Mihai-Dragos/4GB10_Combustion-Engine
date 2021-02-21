%%
%%For full_load dataset

%Please change the path below to the path on your computer where the
%dataset .txt files are 
full_table = dlmread("D:\DESKTOP DE MUTAT INAPOI\TUe University\Uni courses\3.3 Third best Year\Q3\4GB10 Combustion Engine\Github repo\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\full_load.txt");
full_time = full_table(:,1); %time
full_voltage = full_table(:,2); %voltage
full_pulse = full_table(:,3); %pulse



V_s = 5; %V, constant

V_full = full_voltage * V_s;
%2nd voltage 
%3rd pulse

%1 complete cycle is 0.04s = 1/25

P_full = ((V_full/V_s-0.115)/0.0385); 


%%
%%For half_load dataset


%Please change the path below to the path on your computer where the
%dataset .txt files are
half_table = dlmread("D:\DESKTOP DE MUTAT INAPOI\TUe University\Uni courses\3.3 Third best Year\Q3\4GB10 Combustion Engine\Github repo\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\half_load.txt");
half_time = half_table(:,1); %time
half_voltage = half_table(:,2); %voltage
half_pulse = half_table(:,3); %pulse



V_s = 5; %V, constant

V_half = half_voltage * V_s;


%1 complete cycle is 0.04s = 1/25

P_half = ((V_half/V_s-0.115)/0.0385); 

%%
%%For no_load dataset


%Please change the path below to the path on your computer where the
%dataset .txt files are
no_table = dlmread("D:\DESKTOP DE MUTAT INAPOI\TUe University\Uni courses\3.3 Third best Year\Q3\4GB10 Combustion Engine\Github repo\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\no_load.txt");
no_time = no_table(:,1); %time
no_voltage = no_table(:,2); %voltage
no_pulse = no_table(:,3); %pulse



V_s = 5; %V, constant


V_no = no_voltage * V_s;


%1 complete cycle is 0.04s = 1/25

P_no = ((V_no/V_s-0.115)/0.0385); 

