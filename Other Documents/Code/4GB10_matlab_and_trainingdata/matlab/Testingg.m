%%
%%For full_load dataset


full_table = dlmread("C:\Users\alexa\OneDrive\Desktop\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\full_load.txt");
full_time = full_table(:,1);
full_voltage = full_table(:,2);
full_pulse = full_table(:,3);
%full_volume = full_table(:,4);



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

half_table = dlmread("C:\Users\alexa\OneDrive\Desktop\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\half_load.txt");

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


no_table = dlmread("C:\Users\alexa\OneDrive\Desktop\4GB10_Combustion-Engine\Other Documents\Code\4GB10_matlab_and_trainingdata\Data\Training Set\no_load.txt");


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



%% Formula for angle and volume
tetha_0 = 180;
v_tetha= 18000;
k=1;
for i=0:full_time
    tetha = v_tetha*full_time(i) + tetha_0;
if tetha <= 360
    tetha =tetha - 360*k;
    k=k+1;
end
end

% d_tetha = l+r-r*cos(tetha)-sqrt(l^2-r^2*sin*(tetha)^2);

%V_tetha = pi*(B/2)^2*d_tetha+V_c;

%% Pulse full load

% figure()
% hold on
% plot(full_time, full_pulse)
% hold off

%% Pulse

figure()
hold on
plot(full_time, full_pulse)
hold off

figure()
hold on
plot(full_time, full_pulse)

title("Full load pulse sensor");
xlim([0 0.04])
hold off
%% Pulse graphs half load
% figure()
% hold on
% plot(half_time, half_pulse)
% hold off

figure()
hold on
plot(half_time, half_pulse)
title("Half load pulse sensor");
xlim([0 0.04])
hold off

%% Pulse graphs no load
% figure()
% hold on
% plot(no_time, no_pulse)
% hold off

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

title("No load pulse sensor");
xlim([0 0.04])
hold off

%% Pressure graphs no load
% figure()
% hold on
% plot(no_time, P_no)
% hold off

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

title("No load pressure sensor");
hold off

%% Pressure graphs half load
% figure()
% hold on
% plot(half_time, P_half)
% hold off

figure()
hold on
x=1;
plot(half_time, P_half)
plot(x)
xlim([0 0.04])
title("Half load pressure sensor");
hold off

%% Pressure graphs full load

% figure()
% hold on
% plot(full_time, P_full)
% hold off

figure()
hold on
plot(full_time, P_full)
xlim([0 0.04])
title("Full load pressure sensor");

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

% figure()
% hold on
% plot(full_pulse,P_full)
% xlim([0 0.04])
% hold off
% 
% 
%%
% figure()
% hold on
% plot(full.Vol,P_full)
%  hold off

figure()
hold on
plot(full.Vol,P_full)

hold off



