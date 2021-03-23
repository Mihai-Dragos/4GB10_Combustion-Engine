h1 = openfig('test.fig','reuse'); % open figure
ax1 = gca; % get handle to axes of figure
h2 = openfig('E0_fl_mean.fig','reuse');
ax2 = gca;
h3 = openfig('test.fig','reuse'); % open figure
ax3 = gca; % get handle to axes of figure
% test1.fig and test2.fig are the names of the figure files which you would % like to copy into multiple subplots
h5 = figure; %create new figure
s1 = subplot(3,1,1); %create and get handle to the subplot axes
s2 = subplot(3,1,2);
s3 = subplot(3,1,3);

fig1 = get(ax1,'children'); %get handle to all the children in the figure
fig2 = get(ax2,'children');
fig3 = get(ax3,'children');
copyobj(fig1,s1); %copy children to new parent axes i.e. the subplot axes
copyobj(fig2,s2);
copyobj(fig2,s3);