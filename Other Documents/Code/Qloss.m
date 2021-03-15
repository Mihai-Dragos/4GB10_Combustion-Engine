%%
%Woschni
B = 0.067;
C1 = 6.18;
C2 = 0;
Tmax = 1030;
pmax = 2000;

 

%%
RpS = 3000/60;
S = 0.056;
Sp= 2*S*RpS;
w = C1*Sp; %+ C2*((Vd*Tr)/(pr*Vr))*(p-pm);%%

hc = 3.26*(B^-0.2)*(pmax^0.8)*(Tmax^-0.55)*(w^0.8) %W/m^2*K

    
 %% Heat loss formula from handbook
 Twall = 373.15;
 A = 0.0186 %[Surface of inside combustion engine]
 dTloss =  Tmax - Twall; %[Difference between gas temperature and wall temperature]
 Qloss2 = -hc * A * (dTloss)
 
 
 %% Q useful
% Q_useful = Qlhv - Qloss2