%Do not forget to import the "Group07" excel file
%% Task 1
%%Task A

dissolution_table = Group07(:,"Dissolution");
belnd_time_table = Group07(:,"BlendTime");
coating_viscosity_table = Group07(:,"CoatingViscosity");
tempinlet_table = Group07(:,"TempInlet");
belnd_speed_table = Group07(:,"BlendSpeed");

dissolution = table2array(dissolution_table);
BlendTime= table2array(belnd_time_table);
BlendSpeed= table2array(belnd_speed_table);
CoatingViscosity= table2array(coating_viscosity_table);
TempInlet= table2array(tempinlet_table);


figure;
plot(dissolution,'o')
xlabel('Index');
ylabel('Result');
grid on;
hold on;

figure();
subplot(2,2,1);
histogram(dissolution)
title("Histogram of the dissolution of tablets");
xlabel('Dissolution');
ylabel('Frequency');
grid on;
%hold on;

subplot(2,2,2);
ksdensity(dissolution)
title("Kernel density plot of the dissolution of tablets");
xlabel('Dissolution');
ylabel('Frequency');
grid on;
%hold on;

subplot(2,2,3);
boxplot(dissolution) 
title("Box plot dissolution");
xlabel('Result'); 
grid on;

summary = ["Data points:" length(dissolution);"Minimum:" min(dissolution);"Maximum:" max(dissolution);"Range:" range(dissolution);"Mean:" mean(dissolution);"Median:" median(dissolution);"Standard deviation:" std(dissolution);"Variance:" var(dissolution);"Inter quartile range:" iqr(dissolution)];



%%Task C
log_ratio_dissolution_table = log( (dissolution./100) ./ (ones(87) - dissolution./100) );
log_ratio_dissolution = log_ratio_dissolution_table(:,1);

log_ratio_std = std(log_ratio_dissolution);
log_ratio_n = length(log_ratio_dissolution);
log_ratio_avg = mean(log_ratio_dissolution);

%% Task 2

%%Task D

EDA_dissolution=[min(dissolution),max(dissolution),mean(dissolution),median(dissolution),std(dissolution),var(dissolution),iqr(dissolution),length(dissolution)]
EDA_BlendTime=[min(BlendTime),max(BlendTime),mean(BlendTime),median(BlendTime),std(BlendTime),var(BlendTime),iqr(BlendTime),length(BlendTime)]
EDA_BlendSpeed=[min(BlendSpeed),max(BlendSpeed),mean(BlendSpeed),median(BlendSpeed),std(BlendSpeed),var(BlendSpeed),iqr(BlendSpeed),length(BlendSpeed)]
EDA_CoatingViscosity=[min(CoatingViscosity),max(CoatingViscosity),mean(CoatingViscosity),median(CoatingViscosity),std(CoatingViscosity),var(CoatingViscosity),iqr(CoatingViscosity),length(CoatingViscosity)]
EDA_BlendTime=[min(TempInlet),max(TempInlet),mean(TempInlet),median(TempInlet),std(TempInlet),var(TempInlet),iqr(TempInlet),length(TempInlet)]

figure(1);
boxplot(dissolution) 
title("Box plot dissolution");
xlabel('Result'); 
grid on;

figure(2);
subplot(2,2,1);
boxplot(BlendTime);
xlabel('Result');
title('Blend time');
grid on;

subplot(2,2,2);
boxplot(BlendSpeed);
xlabel('Result');
title('Blend speed');
grid on;

subplot(2,2,3);
boxplot(CoatingViscosity);
xlabel('Result');
title('Coating viscosity');
grid on;

subplot(2,2,4);
boxplot(TempInlet);
xlabel('Result');
title('Temperature inlet');
grid on;
hold on;

figure(3);
subplot(2,2,1);
scatter(BlendTime, dissolution); xlabel('BlendTime');ylabel('dissolution');
subplot(2,2,2);
scatter(BlendSpeed, dissolution); xlabel('BlendSpeed');ylabel('dissolution');
subplot(2,2,3);
scatter(CoatingViscosity, dissolution); xlabel('CoatingViscosity');ylabel('dissolution');
subplot(2,2,4);
scatter(TempInlet, dissolution); xlabel('TempInlet');ylabel('dissolution');

corr([BlendTime,dissolution])
corr([BlendSpeed,dissolution])
corr([CoatingViscosity,dissolution])
corr([TempInlet,dissolution])


%% Task E

RegModel0502=fitlm(CoatingViscosity,dissolution)

figure();
plot(RegModel0502)
grid on

coefCI(RegModel0502)
%% Residuals
figure();
plotResiduals(RegModel0502,'fitted')
plotResiduals(RegModel0502,'caseorder')

%% Polynomial model

Example0503=table(CoatingViscosity,dissolution);
RegModel0503=fitlm(Example0503,'dissolution~CoatingViscosity+CoatingViscosity^2')

figure();
plot(RegModel0503)
grid on

coefCI(RegModel0503)

figure();
plotResiduals(RegModel0503,'fitted')
plotResiduals(RegModel0503,'caseorder')

%% This is just me trying do EDA, you can ingnore it
% figure(1);
% 
% subplot(2,2,1);
% plot(BlendTime,'o');
% hold on;
% plot(dissolution,'o');
% legend("BlendTime", "dissolution");
% xlabel('Index');
% ylabel('Result values');
% grid on;
% 
% subplot(2,2,2);
% plot(BlendSpeed,'o');
% hold on;
% plot(dissolution,'o');
% legend("Blend speed", "dissolution");
% xlabel('Index');
% ylabel('Result values');
% grid on;
% 
% subplot(2,2,3);
% plot(CoatingViscosity,'o');
% hold on;
% plot(dissolution,'o');
% legend("Coating viscosity", "dissolution");
% xlabel('Index');
% ylabel('Result values');
% grid on;
% 
% subplot(2,2,4);
% plot(TempInlet,'o');
% hold on;
% plot(dissolution,'o');
% legend("Temperature inlet", "dissolution");
% xlabel('Index');
% ylabel('Result values');
% grid on;
% 
% %histograms
% 
% figure(2);
% subplot(2,2,1);
% histogram(BlendTime);
% hold on;
% histogram(dissolution);
% legend("BlendTime", "dissolution");
% xlabel('Frequency');
% ylabel('Blend time');
% grid on;
% 
% subplot(2,2,2);
% histogram(BlendSpeed);
% hold on;
% histogram(dissolution);
% legend("Blend speed", "dissolution");
% xlabel('Frequency');
% ylabel('Blend speed');
% grid on;
% 
% subplot(2,2,3);
% histogram(CoatingViscosity);
% hold on;
% histogram(dissolution);
% legend("viscosity", "dissolution");
% xlabel('Frequency');
% ylabel('Coating viscosity');
% grid on;
% 
% subplot(2,2,4);
% histogram(TempInlet);
% hold on;
% histogram(dissolution);
% legend("temp", "dissolution");
% xlabel('Frequency');
% ylabel('Temperature inlet');
% grid on;
% 
% 
% %kdens
% 
% figure(3);
% subplot(2,2,1);
% ksdensity(BlendTime);
% hold on;
% ksdensity(dissolution);
% legend("BlendTime", "dissolution");
% xlabel('Density');
% ylabel('Blend time');
% grid on;
% 
% subplot(2,2,2);
% ksdensity(BlendSpeed);
% hold on;
% ksdensity(dissolution);
% legend("Blend speed", "dissolution");
% xlabel('Density');
% ylabel('Blend speed');
% grid on;
% 
% subplot(2,2,3);
% ksdensity(CoatingViscosity);
% hold on;
% ksdensity(dissolution);
% legend("viscosity", "dissolution");
% xlabel('Density');
% ylabel('Coating viscosity');
% grid on;
% 
% subplot(2,2,4);
% ksdensity(TempInlet);
% hold on;
% ksdensity(dissolution);
% legend("temp", "dissolution");
% xlabel('Density');
% ylabel('Temperature inlet');
% grid on;
% 
% %boxplot
% 
% figure(4);
% subplot(2,2,1);
% boxplot(BlendTime);
% % hold on;
% % boxplot(dissolution);
% % legend("BlendTime", "dissolution");
% xlabel('Result');
% title('Blend time');
% grid on;
% 
% subplot(2,2,2);
% boxplot(BlendSpeed);
% % hold on;
% % boxplot(dissolution);
% % legend("Blendspeed", "dissolution");
% xlabel('Result');
% title('Blend speed');
% grid on;
% 
% subplot(2,2,3);
% boxplot(CoatingViscosity);
% % hold on;
% % boxplot(dissolution);
% % legend("viscosity", "dissolution");
% xlabel('Result');
% title('Coating viscosity');
% grid on;
% 
% subplot(2,2,4);
% boxplot(TempInlet);
% % hold on;
% % boxplot(dissolution);
% % legend("temperature", "dissolution");
% xlabel('Result');
% title('Temperature inlet');
% grid on;




