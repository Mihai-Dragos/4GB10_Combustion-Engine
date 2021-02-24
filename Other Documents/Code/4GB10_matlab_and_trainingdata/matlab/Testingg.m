%%
%%For full_load dataset

full_table = dlmread("C:\Users\20192303\Documents\Year 2\Quartile 3\Combustion engine\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\full_load.txt");
full_time = full_table(:,1);
full_voltage = full_table(:,2);
full_pulse = full_table(:,3);
full_volume = full_table(:,4);


V_s = 5; %V, constant

V_full = full_voltage * V_s;
%2nd voltage 
%3rd pulse

%1 complete cycle is 0.04s = 1/25

P_full = ((V_full/V_s-0.115)/0.0385); 

figure()
hold on
plot(full.t,P_full)
xlim([0 0.04])
hold off

%%
%%For half_load dataset

half_table = dlmread("C:\Users\20192303\Documents\Year 2\Quartile 3\Combustion engine\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\half_load.txt");
half_time = half_table(:,1);
half_voltage = half_table(:,2);
half_pulse = half_table(:,3);



V_s = 5; %V, constant

V_half = half_voltage * V_s;
%2nd voltage 
%3rd pulse

%1 complete cycle is 0.04s = 1/25

P_half = ((V_half/V_s-0.115)/0.0385); 

figure()
hold on
plot(full.t,P_half)
xlim([0 0.04])
hold off




%%
%%For no_load dataset

no_table = dlmread("C:\Users\20192303\Documents\Year 2\Quartile 3\Combustion engine\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\no_load.txt");
no_time = no_table(:,1);
no_voltage = no_table(:,2);
no_pulse = no_table(:,3);



V_s = 5; %V, constant


V_no = no_voltage * V_s;%2nd voltage 
%3rd pulse

%1 complete cycle is 0.04s = 1/25

P_no = ((V_no/V_s-0.115)/0.0385); 
figure()
hold on
plot(full.t,P_no)
xlim([0 0.04])
hold off


%% Pulse

figure()
hold on
plot(full_time, full_pulse)
hold off

figure()
hold on
plot(full_time, full_pulse)
xlim([0 0.04])
hold off
%%
figure()
hold on
plot(full_time, full_pulse)
hold off

figure()
hold on
plot(full_time, full_pulse)
xlim([0 0.04])
hold off

%%
figure()
hold on
plot(no_time, no_pulse)
hold off

figure()
hold on
plot(no_time, no_pulse)
xlim([0 0.04])
hold off

%%
figure()
hold on
plot(no_time, P_no)
hold off

figure()
hold on
plot(no_time, P_no)
xlim([0 0.04])
hold off

%%
figure()
hold on
plot(half_time, P_half)
hold off

figure()
hold on
plot(half_time, P_half)
xlim([0 0.04])
hold off


%%
figure()
hold on
plot(full_pulse,P_full)
xlim([0 0.04])
hold off


%%
figure()
hold on
plot(full.Vol,P_full)

hold off


