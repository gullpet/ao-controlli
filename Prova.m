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
A=[0 , 1; 0 , -(bb/(Mi*Ei*Ei+Ie))-(( (Ud*Mi*Ei*Ei)/(Mi*Ei*Ei+Ie) )*2*We) ];
B=[0; 1/(Mi*Ei*Ei+Ie)];
C=[0 , 1];
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


figure(1)
h_G2 = bodeplot(GG,{omega_range_min,omega_range_max});
grid on, zoom on;
hold on; title("G con Ie min e max")
h_G_min = bodeplot(GG_Ei_min,{omega_range_min,omega_range_max});
h_G_max = bodeplot(GG_Ei_max,{omega_range_min,omega_range_max});
Legend=["G(s) Ei=0.5";"G(s) Ei=0.4"; "G(s) Ei=0.6"];
legend(Legend);
hold off;

%% Graphical constraints of the G

Xi = 0.8261;
Mf = Xi * 100;
phase_L_jwc_min = -180 + Mf;
omega_c_min = 460 / (Ta_1 * Mf);


figure(2)
[Mag,phase,omega] = bode(GG, {omega_range_min, omega_range_max}, 'k');

patch([omega_range_min, omega_c_min, omega_c_min, omega_range_min], [0, 0, -200, -200], 'r', 'FaceAlpha', 0.2); grid on
patch([omega_n, omega_range_max, omega_range_max, omega_n], [Bn, Bn, 200, 200], 'r', 'FaceAlpha', 0.4); grid on


Legend_dB=["Attenuazione disturbo in uscita";
    "Attenuazione disturbo di musura" ];
legend(Legend_dB);

hold on;
margin(Mag,phase,omega);grid on;

patch([omega_c_min, omega_n, omega_n, omega_c_min], [-200, -200, phase_L_jwc_min, phase_L_jwc_min], 'r', 'FaceAlpha', 0.3); grid on
Legend_arg=["G(j\omega)";"M_{f,min}" ];
legend(Legend_arg);

title("Diagramma di bode del sistema linearizzato");
hold off;

%% Graphical constraints of the step response

W=5;
S_100=0.01;
T_simulation=15;

figure(3); 

step(GG,T_simulation,'b')
hold on
% add overshoot constraint
patch([0,T_simulation,T_simulation,0],[W*(1+S_100),W*(1+S_100),W+1,W+1],'red','FaceAlpha',0.3);
hold on;
ylim([0,W+1]);


% add Settling time constraint

patch([Ta_1,T_simulation,T_simulation,Ta_1],[W*(1+0.01),W*(1+0.01),W+1,W+1], 'green','FaceAlpha',0.2);

patch([Ta_1,T_simulation,T_simulation,Ta_1],[W*(1-0.01),W*(1-0.01),0,0], 'green','FaceAlpha', 0.2);
Legend=["G(s) Step Response";"Overshoot Constraint"; "Settling time Cons"];
legend(Legend);

title("Risposta di G(s) a ingresso a gradino");
hold off;


%% Graphical constraints of the Ge

GGe        = (1/s) * GG;
GGe_Ei_min = (1/s) * GG_Ei_min;
GGe_Ei_max = (1/s) * GG_Ei_max;

figure(4)

[Mag,phase,omega] = bode(GGe, {omega_range_min, omega_range_max}, 'k');

patch([omega_range_min, omega_c_min, omega_c_min, omega_range_min], [0, 0, -200, -200], 'r', 'FaceAlpha', 0.2); grid on
patch([omega_n, omega_range_max, omega_range_max, omega_n], [Bn, Bn, 200, 200], 'r', 'FaceAlpha', 0.4); grid on

Legend_dB=["Att. disturbo uscita";
    "Att. disturbo misura" ];

legend(Legend_dB);
hold on;
margin(Mag,phase,omega);grid on;

%non si vedrà, troppo in basso nel grafico
patch([omega_c_min, omega_n, omega_n, omega_c_min], [-200, -200, phase_L_jwc_min, phase_L_jwc_min], 'r', 'FaceAlpha', 0.3); grid on 

Legend_arg=["G_{e}(j\omega)";"M_{f,min}"] ;
legend(Legend_arg);
title("G(s) estesa");
hold off;

%% Rete Anticipatrice
clc

omega_c_star=omega_c_min+0.3;

casual_n=( (bb/(dist_mas_iner)) + (( (Ud*Mi*Ei*Ei)/(dist_mas_iner) )*2*We) );
GGe_omega_star = (1/(1i*omega_c_star) )*(1/( (dist_mas_iner) * ( (1i*omega_c_star) +  casual_n) ));

ph=rad2deg(angle(GGe_omega_star));
amplitude=mag2db(abs(GGe_omega_star));

phi_star=(Mf+5)-180-ph;
M_star=10.^(-amplitude/20);


disp("maggiore di 1");
display(M_star);
disp("tra 0 e pi/2");
display(phi_star);
disp("sopra > sotto");
display(cos(phi_star));
display(1/M_star);

tau=(M_star-cos(phi_star))/(omega_c_star*sin(phi_star));
alpha_tau=(cos(phi_star)-1/M_star)/(omega_c_star*sin(phi_star));

R_d=(1+tau*s)/(1+alpha_tau*s);
figure(5);
bode(R_d);
title("Rete Anticipatrice");
%%
gain=1.95;

LL = GGe * R_d * gain;
LL_Ei_min = GGe_Ei_min * R_d * gain;
LL_Ei_max = GGe_Ei_max * R_d * gain;

figure(6);

[Mag,phase,omega] = bode(LL, {omega_range_min, omega_range_max});
[Mag_min,phase_min,omega_min] = bode(LL_Ei_min, {omega_range_min, omega_range_max});
[Mag_max,phase_max,omega_max] = bode(LL_Ei_max, {omega_range_min, omega_range_max});


patch([omega_range_min, omega_c_min, omega_c_min, omega_range_min], [0, 0, -200, -200], 'r', 'FaceAlpha', 0.2); grid on
patch([omega_n, omega_range_max, omega_range_max, omega_n], [Bn, Bn, 200, 200], 'r', 'FaceAlpha', 0.4); grid on
hold on;

Legend_dB=["Att. dist. uscita";
    "Att. dist. misura"
     ];
 
legend(Legend_dB);

margin(Mag,phase,omega);grid on;
bode( LL_Ei_min, 'm', LL_Ei_max, 'k');

patch([omega_c_min, omega_n, omega_n, omega_c_min], [-200, -200, phase_L_jwc_min, phase_L_jwc_min], 'r', 'FaceAlpha', 0.3); grid on
Legend_arg=["L(j\omega)";
    "L_{Ei min}";
    "L_{Eimax}";
    "M_{f,min}"] ;
legend(Legend_arg);

title("L con le variazioni su Ei");
hold off;

%% Graphical constraints of the step response
figure(7); 

step(W*LL/(1+LL),T_simulation,'b')
hold on
% add overshoot constraint
patch([0,T_simulation,T_simulation,0],[W*(1+S_100),W*(1+S_100),W+1,W+1],'red','FaceAlpha',0.3);
hold on;
ylim([0,W+1]);


% add Settling time constraint

patch([Ta_1,T_simulation,T_simulation,Ta_1],[W*(1+0.01),W*(1+0.01),W+1,W+1], 'green','FaceAlpha',0.2);

patch([Ta_1,T_simulation,T_simulation,Ta_1],[W*(1-0.01),W*(1-0.01),0,0], 'green','FaceAlpha', 0.2);
Legend=["L(s) Step Response";"Overshoot Constraint"; "Settling time Cons"];
legend(Legend);
title("Risposta di L(s) a ingresso a gradino");
hold off;

%% Graphical constraints of the step response
figure(8); 

step(W*LL/(1+LL),T_simulation,'b')
step(W*LL_Ei_min/(1+LL_Ei_min),T_simulation,'b')
step(W*LL_Ei_max/(1+LL_Ei_max),T_simulation,'b')
hold on
% add overshoot constraint
patch([0,T_simulation,T_simulation,0],[W*(1+S_100),W*(1+S_100),W+1,W+1],'red','FaceAlpha',0.3);
hold on;
ylim([0,W+1]);


% add Settling time constraint

patch([Ta_1,T_simulation,T_simulation,Ta_1],[W*(1+0.01),W*(1+0.01),W+1,W+1], 'green','FaceAlpha',0.2);

patch([Ta_1,T_simulation,T_simulation,Ta_1],[W*(1-0.01),W*(1-0.01),0,0], 'green','FaceAlpha', 0.2);

Legend=["L(s) Step Response";"Overshoot Constraint"; "Settling time Cons"];
legend(Legend);
title("Risposta di L(s) a ingresso a gradino");
hold off;

%% Simulink
clc

%open("Progetto3cSimulink")

R=(R_d*gain)/s;
[n_r,d_r]=tfdata(R);
num_r=n_r{1};
den_r=d_r{1};

teta_e=0; %non modifica il sistema
x0=[0,We];
u_e=1110.7222;
y_e=We;


