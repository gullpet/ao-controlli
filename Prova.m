clear all, clc
%% Caratteristiche

Mi = 0.1; %[kg] Mass
Ei = 0.5; 
bb = 0.3; %[N/m^2] Dynamic Friction coefficient
Ie = 0.75;
Ud = 0.01;
We = 10000/(2*pi);

Ta_1 = 6;
Bn = -30;
omega_n = 100;

%Ie variations +0.1 -0.1
Ei_min=Ei-0.1;
Ei_max=Ei+0.1;


%%System dynamics
A=[0, 1; 0, (bb/(Mi*Ei*Ei+Ie))-(( (Ud*Mi*Ei*Ei)/(Mi*Ei*Ei+Ie) )*2*We) ];
B=[0, 1/(Mi*Ei*Ei+Ie)];
C=[0, 1];
D=0;
%% Equivalent transfer function
omega_range_min=10^(-2);
omega_range_max=10^4;

dist_mas_iner = Mi*Ei*Ei+Ie;
dist_mas_iner_min = Mi*Ei_min*Ei_min+Ie;
dist_mas_iner_max = Mi*Ei_max*Ei_max+Ie;

s=tf('s');
GG        = 1/( (dist_mas_iner) * ( s + ( (bb/(dist_mas_iner)) + (( (Ud*Mi*Ei*Ei)/(dist_mas_iner) )*2*We) ) ) );
GG_Ei_min = 1/( (dist_mas_iner_min) * ( s + ( (bb/(dist_mas_iner_min)) + (( (Ud*Mi*Ei_min*Ei_min)/(dist_mas_iner_min) )*2*We) ) ) );
GG_Ei_max = 1/( (dist_mas_iner_max) * ( s + ( (bb/(dist_mas_iner_max)) + (( (Ud*Mi*Ei_max*Ei_max)/(dist_mas_iner_max) )*2*We) ) ) );
%{
figure(12)
h_G = bodeplot(GG,{omega_range_min,omega_range_max});
grid on, zoom on;
hold on; title("G")
hold off;

figure(21)
h_G2 = bodeplot(GG,{omega_range_min,omega_range_max});
grid on, zoom on;
hold on; title("G con Ie min e max")
h_G_min = bodeplot(GG_Ei_min,{omega_range_min,omega_range_max});
h_G_max = bodeplot(GG_Ei_max,{omega_range_min,omega_range_max});
Legend=["G(s) Ei=0.5";"G(s) Ei=0.4"; "G(s) Ei=0.6"];
legend(Legend);
hold off;
%
figure(3)
hold on;title("G")
hold off;
margin(GG)
%
figure(4)
step(GG)
title("GG Step Response")

%% punto 4

%}

%% Graphical constraints of the G

Xi = 0.8261;
Mf = Xi * 100;
phase_L_jwc_min = -180 + Mf;
omega_c_min = 460 / (Ta_1 * Mf);

figure(83)

[Mag,phase,omega] = bode(GG, {omega_range_min, omega_range_max}, 'k');

patch([omega_range_min, omega_c_min, omega_c_min, omega_range_min], [0, 0, -200, -200], 'r', 'FaceAlpha', 0.3); grid on
patch([omega_n, omega_range_max, omega_range_max, omega_n], [Bn, Bn, 200, 200], 'r', 'FaceAlpha', 0.3); grid on

Legend_dB=["A_n";
    "\omega_{c,min}"
    "G(j\omega)" ];

legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;

patch([omega_c_min, omega_n, omega_n, omega_c_min], [-200, -200, phase_L_jwc_min, phase_L_jwc_min], 'r', 'FaceAlpha', 0.3); grid on

hold off;

%% Graphical constraints of the step response
%{
W=5;
S_100=0.01;
T_simulation=15;

% figure(52); 

step(GG,T_simulation,'b')
hold on
% add overshoot constraint
patch([0,T_simulation,T_simulation,0],[W*(1+S_100),W*(1+S_100),W+1,W+1],'red','FaceAlpha',0.3);
hold on;
ylim([0,W+1]);


% add Settling time constraint

patch([Ta_1,T_simulation,T_simulation,Ta_1],[W*(1+0.01),W*(1+0.01),W+1,W+1], 'green','FaceAlpha',0.2);

patch([Ta_1,T_simulation,T_simulation,Ta_1],[W*(1-0.01),W*(1-0.01),0,0], 'green','FaceAlpha', 0.2);

hold off
Legend=["G(s) Step Response";"Overshoot Constraint"; "Settling time Cons"];
legend(Legend);
%}

%% Graphical constraints of the Ge

GGe        = (1.5/s)*(1/( (dist_mas_iner) * ( s + ( (bb/(dist_mas_iner)) + (( (Ud*Mi*Ei*Ei)/(dist_mas_iner) )*2*We) ) ) ));
% figure(6)
% step(GGe)

figure(126)

[Mag,phase,omega] = bode(GGe, {omega_range_min, omega_range_max}, 'k');

patch([omega_range_min, omega_c_min, omega_c_min, omega_range_min], [0, 0, -200, -200], 'r', 'FaceAlpha', 0.3); grid on
patch([omega_n, omega_range_max, omega_range_max, omega_n], [Bn, Bn, 200, 200], 'r', 'FaceAlpha', 0.3); grid on

Legend_dB=["A_n";
    "\omega_{c,min}"
    "G_e(j\omega)" ];

legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;

patch([omega_c_min, omega_n, omega_n, omega_c_min], [-200, -200, phase_L_jwc_min, phase_L_jwc_min], 'r', 'FaceAlpha', 0.3); grid on

hold off;
