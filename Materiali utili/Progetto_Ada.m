clear all, clc
%% Parameters

mm=0.15; %[kg] Mass
ee=0.47; 
gg=9.8065; 
beta=0.3; %[N/m^2] Dynamic Friction coefficient
Ie=0.75;
%Ie variations +0.1 -0.1
Ie_min=0.65;
Ie_max=0.85;


%%System dynamics
% dot{x_1}=x_2
% dot{x_2}mm=-kk x_1- beta x_2 + u
% y=x_1

A=[0, 1; ((-gg*mm*ee)*cos(pi/2))/((mm*ee*ee+Ie)), -beta/(mm*ee*ee +Ie)];
B=[0;1/(mm*ee*ee + Ie)];
C=[0, 1];
D=0;
%% Equivalent transfer function
omega_plot_min=10^(-3);
omega_plot_max=10^4;

mi = 1;
s=tf('s');
GG=mi/(s*(mm*ee*ee+Ie)+beta);
GG_Ie_min=mi/(s*(mm*ee*ee+Ie_min)+beta);
GG_Ie_max=mi/(s*(mm*ee*ee+Ie_max)+beta);

figure(1)
h_G = bodeplot(GG,{omega_plot_min,omega_plot_max});
grid on, zoom on;
hold on; title("G")
hold off;
h_G2 = bodeplot(GG,{omega_plot_min,omega_plot_max});
grid on, zoom on;
hold on; title("G con Ie min e max")
h_G_min = bodeplot(GG_Ie_min,{omega_plot_min,omega_plot_max});
h_G_max = bodeplot(GG_Ie_max,{omega_plot_min,omega_plot_max});
Legend=["G(s) Ie=0.75";"G(s) Ie=0.85"; "G(s) Ie=0.65"];
legend(Legend);
hold off;
figure(2)
step(GG)
title("GG Step Response")
figure(3)
margin(GG)
Legend=["G(s) Ie=0.75"];
legend(Legend);
%% Requirements/specifications
% -> w=W 1(t)
W=4;% [m] Reference position
% -> Allowed overshoot S%<=1%=0.01
S_100= 0.01;
xi= sqrt(log(S_100)^2/(pi^2+log(S_100)^2));
Mf=xi*100; %in good approximation
display(Mf)
% -> Settling time
%Ta,1<=6 [s]
T_a1Max=6;
%
%T_a1=4.6/(xi*omega_c)<6
%omega_c>4.6/(6*xi)
omega_cMin=4.6/(T_a1Max*xi);
% Measure Noise requirements 
omega_n=120; 
%A_n=-26; %20  %%TO CHECK
A_n=-30;

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
figure(95)
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{c,min}= 0.93"
    "G(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G(j\omega)";"M_{f,min}= 82"] ;
legend(Legend_arg)
hold off
title("G Bode diagram")

% %% Adjusting the gain we cannot acheive the right performance
% G_temp=GG/0.5*W;
% 
% [Mag,phase,omega]=bode(G_temp,{omega_plot_min,omega_plot_max},'k');
% figure()
% patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
% patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
% Legend_dB=["A_n";
%     "\omega_{c,min}= 0.93"
%     "G(s)" ];
% legend(Legend_dB);
% hold on;
% margin(Mag,phase,omega);grid on;
% patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on
% Legend_arg=["G(s)";"M_{f,min}= 82"] ;
% legend(Legend_arg)
% hold off
% step(G_temp)
% title("Step G")

%% Inserting the internal model
Rs_s=1/s;
GG_e=Rs_s*GG;
GGe_Ie_min=Rs_s*GG_Ie_min;
GGe_Ie_max=Rs_s*GG_Ie_max;

[Mag,phase,omega]=bode(GG_e,{omega_plot_min,omega_plot_max},'k');
figure(138)
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'r','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{c,min}"
    "G_e(j\omega)" ];

legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G_e(j\omega)";"M_{f,min}"] ;
legend(Legend_arg)
hold off
title("Extended G Bode diagram")
%% Rete dinamica
%% Rete anticipatrice

%R_ant=(1+tau*s)/(1+atau*s), a<1
%tau= (M_star - cos(phi_star))/(omega_cStar * sin(phi_star))
omega_cStar=2.83; %% Insert here
phi_Ge = -172; %Insert here
Ge_dB = -16; %Insert here

%phi_star= (Mf_spec+10) - 180 -(-180);display(phi_star);
phi_star= (Mf + 2) - 180 - phi_Ge; 
display(phi_star);


espon = (-Ge_dB)/20;
M_star = 10.^espon;

display(M_star)
M_rec = 1/M_star;%cos(phi_star)>1/M_star
display(M_rec)
cos_phi = cosd(phi_star);
display(cos_phi)

tau= (M_star - cosd(phi_star))/(omega_cStar * sind(phi_star));
atau=(cosd(phi_star) - 1/M_star)/(omega_cStar * sind(phi_star));
display(tau)
display(atau)

R_ant=(1+tau*s)/(1+atau*s);
figure()
bode(R_ant)
title("Rete anticipatrice")

%% G_e R_ant
[Mag,phase,omega]=bode(GG_e*R_ant,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{c,min}= 0.93"
    "G_e R_{ant}(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["G_e R_{ant}(j\omega)";"M_{f,Min}= 82"] ;
legend(Legend_arg)
hold off
title("Extended G + R_{ant} Bode diagram")
zoom on

figure()
LL=GG_e*R_ant;
step(LL/(1+LL))
title("Step F")


%% L open loop

[Mag,phase,omega]=bode(LL,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{c,min}= 0.93"
    "L(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["L(j\omega)";"M_{f,min}= 82"] ;
legend(Legend_arg)
hold off
title("Open loop L Bode diagram")
zoom on

%% L graph with Ie variations

[Mag,phase,omega]=bode(LL,{omega_plot_min,omega_plot_max},'k');
figure()
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[A_n,A_n,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
patch([omega_plot_min,omega_cMin,omega_cMin,omega_plot_min],[0,0,-150,-150],'o','FaceAlpha',0.3,'EdgeAlpha',0); grid on
Legend_dB=["A_n";
    "\omega_{c,min}= 0.93"
    "L(j\omega)" ];
legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on; hold on
Ie1_LL=GGe_Ie_min*R_ant;
Ie1h_LL = bodeplot(Ie1_LL); hold on
 
Ie2_LL=GGe_Ie_max*R_ant;
Ie2h_LL = bodeplot(Ie2_LL); hold on
patch([omega_cMin,omega_cMax,omega_cMax,omega_cMin],[-180+Mf,-180+Mf,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); grid on

Legend_arg=["L(j\omega)";"L(j\omega_{Ie,min=0.65}";"L(j\omega_{Ie,max=0.85}";"M_{f,min}= 82"] ;
legend(Legend_arg)
hold off
title("Open loop L with Ie variations Bode diagram")
zoom on


%% F 
FF = LL/(1+LL);
h_FF = bodeplot(FF,{omega_plot_min,omega_plot_max})
title("Rete F");grid on; zoom on;

figure()
step(FF,10)
grid on
zoom on
title("F Step")


%% RECAP
%progetto in tre fasce
%-> Bassa frequenza : static gain (abbattimento disturbo sull'uscita/ inseguimento riferimento)
%-> Media frequenza : loop shaping (specifiche di performance, sovraelongazione e tempo di assestamento)
%-> Alta frequenza : poli in alta frequenza (riduzione rumore di misura da sensore)

%% Simulink addictional part
open("Progetto4cSimulink")
% Define the controller coefficients to run the simulink file
R=R_ant/s;
[n_r,d_r]=tfdata(R);
num_r=n_r{1};
den_r=d_r{1};
%initial condition
x0=[0;0];