clear all, clc
%% Parameters
mm=0.5; %[kg] Mass
beta=0.3; %[N/m^2] Dynamic Friction coefficient
kk=1; % [N/m] Spring stiffness constant

%%System dynamics
% dot{x_1}=x_2
% dot{x_2}mm=-kk x_1- beta x_2 + u
% y=x_1

A=[0,1; -kk/mm, -beta/mm];
B=[0;1];
C=[1,0];
D=0;
%% Equivalent transfer function
s=tf('s');
GG=1/(s^2+beta/mm*s+kk/mm);

figure(1)
bode(GG)
figure(2)
step(GG)
figure(3)
margin(GG)
%% Requirements/specifications
% -> w=W 1(t)
W=2;% [m] Reference position
% -> Allowed overshoot S%<=3%=0.03
S_100= 0.03;
xi= sqrt(log(S_100)^2/(pi^2+log(S_100)^2));
Mf=xi*100; %in good approximation
display(Mf)
% -> Settling time
%Ta,1<=2 [s]
T_a1Max=2;
%
%T_a1=4.6/(xi*omega_c)<2
%omega_c>4.6/(2*xi)
omega_cMin=4.6/(T_a1Max*xi);
% Measure Noise requirements 
omega_n=1500;
A_n=-26; %20 

%% Graphical constraints of the step response
figure(2); 
T_simulation=20;
step(GG,T_simulation,'b')
hold on
% add overshoot constraint
patch([0,T_simulation,T_simulation,0],[W*(1+S_100),W*(1+S_100),W+1,W+1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
hold on; ylim([0,W+1]);
% add Settling time constraint
patch([T_a1Max,T_simulation,T_simulation,T_a1Max],[W*(1-0.01),W*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a1Max,T_simulation,T_simulation,T_a1Max],[W*(1+0.01),W*(1+0.01),W+1,W+1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

hold off
Legend=["G(s) Step Response";"Overshoot Constraint"; "Settling time Cons"];
legend(Legend);
%% Graphical constraints on Bode's plot
omega_plot_min=10^(-1);
omega_plot_max=10^4;

omega_cMax=omega_n;

[Mag,phase,omega]=bode(GG,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "G(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Nominal plant TF Bode diagram")

% % %% Adjusting the gain we cannot acheive the right performance
% % G_temp=GG/0.5*W;
% % 
% % [Mag,phase,omega]=bode(G_temp,{omega_plot_min,omega_plot_max},'k');
% % figure()
% % patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
% % patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
% % Legend_dB=["A_n";
% %     "\omega_{cMin}"
% %     "G(s)" ];
% % legend(Legend_dB);
% % hold on;
% % margin(Mag,phase,omega);grid on;
% % patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on
% % Legend_arg=["G(s)";"M_f"] ;
% % legend(Legend_arg)
% % hold off
% % step(G_temp)

%% Inserting the internal model
GG_e=GG/s;
[Mag,phase,omega]=bode(GG_e,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "G_e(j\omega)" ];

legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G_e(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Extended plant TF Bode diagram")
%% Prima parte della rete dinamica
%recupero di fase maggiore di 90°, devo avere due zeri o due reti anticipatrici

%we can add for free a zero
tau_z=1/(omega_cMin*0.5); %= 1/omega_z
R_1=1+tau_z*s; %GG_e*R_1=G(s)* R_1/s; thus respecting the causality property

[Mag,phase,omega]=bode(GG_e*R_1,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "G_e R_1(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G_e R_1(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Plant with PI TF Bode diagram")

%% Rete anticipatrice

%R_ant=(1+tau*s)/(1+atau*s), a<1
%tau= (M_star - cos(phi_star))/(omega_cStar * sin(phi_star))
omega_cStar=(omega_cMin+50);
phi_star= (Mf+10) - 180 -(-180);display(phi_star);
M_star=20;% so to fulfill cos(phi_star)>1/M_star
display(1/M_star)
display (cosd(phi_star))

tau= (M_star - cosd(phi_star))/(omega_cStar * sind(phi_star));
atau=(cosd(phi_star) - 1/M_star)/(omega_cStar * sind(phi_star));

R_ant=(1+tau*s)/(1+atau*s);
figure()
bode(R_ant)
title("Rete anticipatrice")

%% 
[Mag,phase,omega]=bode(GG_e*R_1*R_ant,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "G_e R_1 R_{ant}(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G_e R_1 R_{ant}(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Plant with PI+ R_{ant} TF Bode diagram")
zoom on

figure()
L=GG*R_1/s*R_ant;
step(L/(1+L))

%% Define the gain to adjust the cross frequency
gain=1/abs(evalfr(GG_e*R_1*R_ant,1i*omega_cStar));
L=GG_e*R_1*R_ant*gain;


[Mag,phase,omega]=bode(L,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "L(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["L(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Open loop TF Bode diagram")
zoom on

figure()
step(L/(1+L),5)
grid on
zoom on % ylim([0.94 1.07]);xlim([0 2])
%% USE CONTROL SYSTEM DESIGNER
%addictional gain
gain_CSD=2.01;
L_CSD=L*gain_CSD;
figure()
[y_step,t_step]=step(W*L_CSD/(1+L_CSD),5);
plot(t_step,y_step)
patch([0,5,5,0],[W*(1+S_100),W*(1+S_100),W+1,W+1],'r','FaceAlpha',0.3,'EdgeAlpha',0.5);
hold on; ylim([0,W+1]);
% add Settling time constraint
patch([T_a1Max,5,5,T_a1Max],[W*(1-0.01),W*(1-0.01),0,0],'g','FaceAlpha',0.1,'EdgeAlpha',0.5);
patch([T_a1Max,5,5,T_a1Max],[W*(1+0.01),W*(1+0.01),W+1,W+1],'g','FaceAlpha',0.1,'EdgeAlpha',0.1);

hold off
Legend=["F(s) Step Response";"Overshoot Constraint"; "Settling time Constraint"];
legend(Legend,'location','SouthEast');

zoom on

%% 

[Mag,phase,omega]=bode(L_CSD,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{cMin}"
    "L_{CSD}(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["L_{CSD}(j\omega)";"M_f"] ;
legend(Legend_arg)
hold off
title("Open loop L_{CSD} TF Bode diagram")
zoom on


%% Addictional noise requirements
omega_n=1300;
%To fulfill noise attenuation requirements work on Control system Designer
%by adding a high frequency pole (trail and error synthesis, 'sintesi per tentativi') 
%or designing a "rete ritardatrice"

%% RECAP
%progetto in tre fasce
%-> Bassa frequenza : static gain (abbattimento disturbo sull'uscita/ inseguimento riferimento)
%-> Media frequenza : loop shaping (specifiche di performance, sovraelongazione e tempo di assestamento)
%-> Alta frequenza : poli in alta frequenza (riduzione rumore di misura da sensore)

%% Simulink addictional part
open("MassaMollaSmorzatore")
% Define the controller coefficients to run the simulink file
R=gain*gain_CSD*R_ant*R_1/s;
[n_r,d_r]=tfdata(R);
num_r=n_r{1};
den_r=d_r{1};
%initial condition
x0=[0;0];