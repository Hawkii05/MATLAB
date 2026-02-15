clear; clc;

%% Problem 1
b = 10.3; %m wing span
e = 0.9; %oswald efficiency factor
S = 16.91; %m^2 wing SA
S_t = 3.902; %m^2 tail SA
m = 2900; %kg total mass with fuel tank 
T_sl = 4.5 * 2 * 1000; %kN two turbo jet engines max thrust
TSFC = 1.2; %lb/lb*hr
C_d0 = 0.02;
C_lmax = 1.4; 
C_macw = 0.01;
C_l0w = 0.1; 
C_law = 6.06;
C_lat = 4.7;
l_t = 4.877; %m
i_t = deg2rad(-1); %rad incidence angle
deda = 0.17; %down wash derivative
e_0 = 0;
rho_sl = 1.225; %kg/m^3
AR = b^2/S;
g=9.81; %m/s^2
W = m * g; %N

%altitude
h = 0; %m 

%density
t_sl = 288.16; %kelvin
a = -6.5*10^-3; %k/m
R = 287; %J/kg*k
t=t_sl-a*h;
rho = rho_sl*(t/t_sl).^(-1-g/(a*R));

%Part a
V = 0:350:50; % m/s velocity range for analysis
k = 1/(pi*e*AR);
Q = 0.5*rho.*V.^2;
RC = V.*(T_sl/(W/S)-Q*C_d0/(W/S)-(W/S)*k/Q); 
RC_SC = 100 * 0.508; %m/s

figure(1);
plot(V, RC, 'b-', 'LineWidth', 2)
hold on;
yline(RC_SC, 'r--', 'LineWidth', 1.5, 'Label', 'Service Ceiling Criterion');
grid on;
xlabel('Velocity (m/s)', 'FontSize', 12);
ylabel('Rate of Climb (m/s)', 'FontSize', 12);
title('Rate of Climb vs Velocity at Sea Level', 'FontSize', 14);
legend('Rate of Climb', 'Service Ceiling');


%Part b 
v = 80; %m/s
q = 0.5*rho_sl*v^2;
n_maxT = sqrt((q / (k * (W/S)) * (T_sl/W - q*C_d0/(W/S)))); 
n_maxA = (0.5*rho_sl*v^2*S*C_lmax)/W;

n_max = min(n_maxT, n_maxA);
R_min = v^2/(g*sqrt(n_max^2-1));

display(R_min); 

%% Part c
alt = 3500; %m
V_cruise_3500 = 1;
WorcesterLA = 4108.928; %km
WorcesterDetroit = 924.077; %km
WorcesterChicago = 1323.040; %km
WorcesterDallas = 2431.120; %km 
WorcesterAtlanta = 1450.110; %km
W_0 = 2900; %kg
W_1 = 1900; %kg

Range = V_cruise_3500/TSFC * L/D * ln(W_0/W_1);

%% Problem 2
%determine the values of angle of attack and elevator angle to trim at the
%airspeed for maximum range (which you found in the previous problem), and at altitude 3.5
%km.