a = 5; %Constant
n = 3; %Constant
theta_s = 0;  % [rad] Start of ignition
theta_d = 15;        % [rad] Ignition duration
f = 5;                      %Constant
dtheta = pi/(f*180);            % [rad] Step in crank angle

i2 = (theta_s)/dtheta + 1; 
i1 = (theta_s + theta_d)/dtheta + 1;

theta_b = i1 - i2;
% Wiebe
x_b = 1-exp(-a*((theta_b-theta_s)/theta_d).^n); %Wiebe function
x_b

%%
%Woschni
B = 0.067;
C1 = 6.18;
C2 = 0;
Tmax = 1031;
pmax = 2000;

 

%%
RpS = 3000/60;
S = 0.056;
Sp= 2*S*RpS;
w = C1*Sp; %+ C2*((Vd*Tr)/(pr*Vr))*(p-pm);%%

hc = 3.26*(B^-0.2)*(pmax^0.8)*(Tmax^-0.55)*(w^0.8) %W/m^2*K

%%
% heat loss through walls:
r = 0.5*S;
l = 0.88  ;                                                                %Length Connecting rod [m]
Cai = 360
%%
    D(i) = (l+r-(r*cosd(Cai)+sqrt(l^2-(r^2*(sind(Cai))^2))));              %length of variable part [m]
    %%
    A(i) = 2*pi*Rin*D(i);                                                  %[m^2]
    Q_dot_cylinder(i) = -1*((k(i)/dx*A(i)*(T(i-1)-T_in))) ;                %[J/s] 
    Q_dot_top(i) = -1*((k(i)/dx*A_top_cylinder*(T(i-1)-T_in)));            %[J/s]
    Q_dot_bot(i)=-1*(k(i)/thickness_bot*A_bot_cylinder*(T(i-1)-T_in));     %[J/s]
    
 %% Heat loss formula from handbook
 alpha = 15; %[Heat coefficient Aluminum]
 A = ? %[Surface of inside combustion engine]
 dTloss =  Tmax - Twall; %[Difference between gas temperature and wall temperature]
 Qloss = -alpha * A * (dTloss)