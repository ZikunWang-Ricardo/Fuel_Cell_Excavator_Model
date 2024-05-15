%% ***********Constants****************************

%Constants Name Unit

load T_DRCY.mat ;
load Final_Drive_cycle.mat;
load Final_Drive_cycle_2.mat;
%% Thermodynamic constants
F = 96485; %Faraday's constant, C/mol
P_atm = 101.25*10^3; %1atm, Pa
T_atm = 298.15; %atmospheric temperature, K
gamma = 1.4; %ratio of specific heats of air
C_p = 1004; %constant pressure specific heat of air, J/(mol*K)

rho_H2 = 0.0899; 
rho_N2 = 1.165;
R = 8.315; %gas constant, J/(mol*K)
R_a = 286.9; %air gas constant, J/(kg*K)
R_O2 = 259.8;
R_N2 = 296.8;
R_v = 461.5;
R_H2 = 4124.3;
M_O2 = 32*10^-3; %oxygen molar mass, kg/mol
M_H2 = 2*10^-3;
M_N2 = 28*10^-3;
M_v = 18*10^-3;
M_a = 29*10^-3;
g = 9.8; %m/s2

%% Ambient constants
L = 0.0065; %Standard temperature gradient, K/m

T_a = T_atm-h*L;
P_a = P_atm*(1-L*h/T_atm)^(g*M_a/(R*L));
phi_a = 0.1; %relative humidity
rho_a = P_a/(R_a*T_a); %air density, kg/m3

y_O2 = 0.21;

y_O2_ca_in = 0.21;

%% Fuel cell constants
%n_cell = 455; %Number of cells powercell p stack
n_cell = 655; %Number of cells EC250FC
A_fc = 550; %Active area of fuel cell, cm2 powercell p stack&EC250FC
V_ca = n_cell*A_fc*10^(-4)*0.0006*0.3; %Cathode volume, m3
V_an = n_cell*A_fc*10^(-4)*0.0005*0.3; %Anode volume, m3
% V_ca = n_cell*A_fc*10^(-4)*0.0006*0.7; %Cathode volume, m3
 %V_an = n_cell*A_fc*10^(-4)*0.0005*0.7; %Anode volume, m3
T_st = 353.15; %fuel cell temperature, K

P_cp_out = 3.6*P_a;
P_cp_out_ref = 3*P_a;
k_ca_out = 0.2177*10^(-5); %cathode outlet orifice constant, kg/(s*Pa)
P_ca_ref=2.5*P_a;
%% Control constants
K1 = 2.1; %propertional gain Kg/s/kPa
K2 = 0.94; %norminal pressure drop

kp_cp=7.5;%compressor control
ki_cp=0.1;
kd_cp=0.01;

%% membrane constants
t_m = 0.01; %membrane thickness, cm
rho_m_dry = 0.002; %menbrane dry density, kg/cm3
M_m_dry = 1.1; %menbrane dry equivalent weight, kg/mol
i_max = 2.2; %saturate current density, A/cm2
%empircial values:
b11 = 0.005139; 
b12 = 0.00326;
b2 = 350;
c1 = 10;
c3 = 2;

%% intergal initial values
m_H2_initial = P_atm*V_an/(R_H2*T_atm);
m_O2_initial = P_atm*y_O2*V_ca/(R_O2*T_atm);
%m_O2_initial = 0.00001;
m_N2_initial = P_atm*(1-y_O2)*V_ca/(R_N2*T_atm);
%m_N2_initial = P_atm*V_ca/(R_N2*T_atm);
m_H2O_initial = 0.005;

%% humidifier & cooler
phi_des = 1; %desired relative humidification of cathode inlet gas
T_cl = 353.15; %target cooled cathode flow temperature, K

%% anode exhaust manifold 
Cd_rm = 0.0124; %exhaust manifold throttle discharge coefficient
A_T_rm = 0.005; %exhaust manifold throttle area, m2
V_rm = 0.005; %exhaust manifold volume, m3
T_rm = 303.15; %exhaust manifold temperature, K

kp_rm=-8E-7; %exhaust throttle control
ki_rm=-5E-8;
kd_rm=0;

%% hydrogen tank pressure release valve


%% anode purge
Cd_pv = 0.0124; %discharge coefficient of the nozzle
A_noz = 0.0008; %nozzle flow area, m2
VP_purge = 0.3; %purge valve position, 0% is closed
R_Vst = 0.15; %errror voltage over non-loss voltage
V_st_min=400;
%% hydrogen recirculation

%% supply manifold
V_sm=0.02;
K_sm_out= 0.3629*10^-5;
%% Thermal Management
Cp_cl = 4180; %Heat capacity of water, J/(kg*K)
Cp_air = 1005; %heat capacity of air, J/(kg*K)
A_st = 1.1*A_fc*1E-4*n_cell*0.0012; %surface area of stack, m2
rho_cl = 1E3; %density of coolant, kg/m3

Cp_Al = 900; %heat capacity of aluminium, J/(kg*K)
rho_Al = 2.7E3; %density of aluminium, kg/m3
m_st = n_cell*A_fc*t_m*rho_m_dry+(V_ca+V_an)*rho_Al*2.5; %mass of stack, kg
Cp_st = Cp_Al*1.5; %heat capacity of stack, J/(kg*K)
H_conv_st = 3; %convection heat transfer coefficient, from stack to air W/(m2*K) 

%T_cl_in = 333; %inlet coolant flow of stack, K
%dV_cl = 0.0025; %vloume flow rate of coolant, m3/s 
T_st_ref = 353.15; %reference operating temperature, K

A_rad = 80; %surface area of radiator, m2

%radiator fan control
kp_fan = -0.5;
ki_fan = -1E-5;
kd_fan = 0;

%pump control
%kp_pump = -1E-3;
%ki_pump = -1E-4;
%kd_pump = -1E-4;

%% ************Variables***************************
%Variables Name Unit
%W: Mass flow rate, kg/s

%i = [linspace(1,10,10); linspace(0,1.5,10)];
%% Drivetrain
% fc = 20000;
% fnom = 100;
% fsw = 2000;
% J = 0.2;
% L0 = 2e-4;
% Ld = 2e-4;
% Lq = 2e-4;
% N = 6;
% PM = 0.03;
% Pmax = 440e3;
% Rs = 0.013;
% Tmax = 1500;
% Ts = 2e-6;
% Tsc = 5e-5;
Ts = 2e-6;
ratio = 2;
v0 = 0;
wheel_radius = 0.3;
Pmax = 267; % kW
Tmax = 2500; %Nm
Ka = 2;
Kb = 50;
cc=30;
charge_value=30;

%% Driving Cycle



