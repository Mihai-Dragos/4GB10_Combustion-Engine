%Do not forget to import the "Group07" excel file

%%Task A
dissolution_table = Group07(:,"Dissolution");
dissolution = table2array(dissolution_table);

figure;
plot(dissolution,'o');
xlabel('Index');
ylabel('Result');
grid on;
hold on;

histogram(dissolution);
title("Histogram of the dissolution of tablets");
xlabel('Dissolution');
ylabel('Frequency');
grid on;
hold on;

ksdensity(dissolution);
title("Kernel density plot of the dissolution of tablets");
xlabel('Dissolution');
ylabel('Frequency');
grid on;
hold on;

boxplot(dissolution); 
xlabel('Result'); 
grid on;

summary = ["Data points:" length(dissolution);"Minimum:" min(dissolution);"Maximum:" max(dissolution);"Range:" range(dissolution);"Mean:" mean(dissolution);"Median:" median(dissolution);"Standard deviation:" std(dissolution);"Variance:" var(dissolution);"Inter quartile range:" iqr(dissolution)];

%%Task C
log_ratio_dissolution_table = log( (dissolution./100) ./ (ones(87) - dissolution./100) );
log_ratio_dissolution = log_ratio_dissolution_table(:,1);

log_ratio_std = std(log_ratio_dissolution);
log_ratio_n = length(log_ratio_dissolution);
log_ratio_avg = mean(log_ratio_dissolution);