% Reactor Simulation Script (Updated with dP_sat/dt calculation)

% Clear workspace and command window
clear all;
close all;
clc;

% Add path to XSteam functions if necessary
% addpath('path_to_XSteam_functions');

%% Constants and Parameters

% Reactivity coefficients
alpha_f = -2.16e-5;    % Fuel temperature coefficient (1/°C)
alpha_c = -1.8e-4;     % Moderator temperature coefficient (1/°C)

% Heat transfer areas and coefficients
A_fc = 583;            % Fuel to coolant area (m²)
h_fc = 1135;           % Fuel to coolant coefficient (W/m²·°C)
A_pm = 1123;           % Primary to SG metal area (m²)
h_pm = 20391;          % Primary to SG metal coefficient (W/m²·°C)
A_ms = 1214;           % SG metal to secondary area (m²)
h_ms = 4950;           % SG metal to secondary coefficient (W/m²·°C)

% Volumes
V_core = 1.879;        % Core volume (m³)
V_HL = 9.7;            % Hot leg volume (m³)
V_CL = 26.8;           % Cold leg volume (m³)
V_SG_primary = 3.564;  % SG primary volume (m³)

% Pressures
P_primary = 12.76;     % Primary pressure (MPa)
P_primary_bar = P_primary * 10;  % Convert to bar for XSteam

% Temperatures
T_fi = 170;            % Feedwater inlet temperature (°C)

% Masses and heat capacities
m_m = 7869;            % SG metal mass (kg)
c_m = 0.450 * 1e3;     % SG metal specific heat (J/kg·°C)
m_f = 11252;           % Fuel mass (kg)
c_pf = 0.467 * 1e3;    % Fuel specific heat (J/kg·°C)

% Rated thermal power
P_0 = 160e6;           % Rated power (W)
tau = 0.97;            % Fraction of power to fuel

% Neutron kinetics parameters
beta = 0.007;         % Delayed neutron fraction
Lambda = 2e-5;         % Prompt neutron lifetime (s)
lambda = 0.1;         % Precursor decay constant (1/s)

% Initial Conditions

% Neutron flux and precursor concentration
phi0 = 0.976655;                      % Normalized neutron flux
C0 = (beta / (Lambda * lambda)) * phi0;  % Steady-state precursor concentration

% Temperatures (°C)
T_f0 = 483.249;                 % Fuel temperature
T_c10 = 254.149;                 % Coolant temperature at Node 1
T_c20 = 278.419;                % Coolant temperature at Node 2
T_HL0 = 278.42;                % Hot leg temperature
T_p0 = 252.761;                  % Primary coolant temperature in SG
T_m0 = 246.483;                  % SG metal temperature
T_sat0 = 222.561;               % Saturation temperature in SG secondary

% Convert saturation temperature to pressure using XSteam
P_sat0_bar = XSteam('pSat_T', T_sat0);  % Saturation pressure (bar)
P_sat0 = P_sat0_bar / 10;               % Convert to MPa

% Mass flow rates
m_cp0 = 590.293;                % Primary mass flow rate (kg/s)
m_cs0 = 77.2787;                 % Secondary mass flow rate (kg/s)

% Initial cold leg temperature (assumed equal to T_p0 for initial steady-state)
T_CL0 = 2*T_p0 - T_HL0;

% K_s value (provided)
K_s = 1e6;  % Given value

%% Parameters structure for ODE function
params.alpha_f = alpha_f;
params.alpha_c = alpha_c;
params.beta = beta;
params.Lambda = Lambda;
params.lambda = lambda;
params.P_0 = P_0;
params.tau = tau;
params.h_fc = h_fc;
params.A_fc = A_fc;
params.m_f = m_f;
params.c_pf = c_pf;
params.V_core = V_core;
params.P_primary_bar = P_primary_bar;
params.T_f0 = T_f0;
params.T_c10 = T_c10;
params.T_c20 = T_c20;
params.m_cp0 = m_cp0;
params.h_pm = h_pm;
params.A_pm = A_pm;
params.m_m = m_m;
params.c_m = c_m;
params.h_ms = h_ms;
params.A_ms = A_ms;
params.m_cs0 = m_cs0;
params.T_fi = T_fi;
params.P_sat0 = P_sat0;
params.V_HL = V_HL;
params.V_CL = V_CL;
params.V_SG_primary = V_SG_primary;
params.K_s = K_s;
params.T_sat0 = T_sat0;  % Initial saturation temperature
params.T_HL0 = T_HL0;
params.T_CL0 = T_CL0;

% Pre-calculate K_sm
params.K_sm = h_ms * A_ms;  % W/°C

%% Initial state vector
y0 = [phi0; C0; T_f0; T_c10; T_c20; T_HL0; T_CL0; T_p0; T_m0; P_sat0];

%% Time span
t_span = [0, 500];  % Simulation from 0 to 200 seconds

%% ODE options
options = odeset('RelTol',1e-6,'AbsTol',1e-8);

%% Solve ODEs
[t, y] = ode15s(@(t, y) reactor_ode(t, y, params), t_span, y0, options);

%% Extract variables
phi = y(:,1);
C = y(:,2);
T_f = y(:,3);
T_c1 = y(:,4);
T_c2 = y(:,5);
T_HL = y(:,6);
T_CL = y(:,7);
T_p = y(:,8);
T_m = y(:,9);
P_sat = y(:,10);

%% Calculate reactivity over time
rho_ext = zeros(length(t),1);
rho_ext(t >= 200) = 3.5e-4;
rho = rho_ext + alpha_f * (T_f - T_f0) + alpha_c * ((T_c1 - T_c10) + (T_c2 - T_c20)) / 2;

%% Plot the results
figure;
plot(t, rho, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Core Reactivity \rho(t)', 'FontSize', 12);
title('Core Reactivity over Time', 'FontSize', 14);
grid on;

%% Plot Fuel Temperature (T_f)
figure;
plot(t, T_f, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Fuel Temperature T_f (°C)', 'FontSize', 12);
title('Fuel Temperature over Time', 'FontSize', 14);
grid on;

%% Plot Hot Leg Temperature (T_HL)
figure;
plot(t, T_HL, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Hot Leg Temperature T_{HL} (°C)', 'FontSize', 12);
title('Hot Leg Temperature over Time', 'FontSize', 14);
grid on;

%% Plot Cold Leg Temperature (T_CL)
figure;
plot(t, T_CL, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Cold Leg Temperature T_{CL} (°C)', 'FontSize', 12);
title('Cold Leg Temperature over Time', 'FontSize', 14);
grid on;

%% Plot Saturation Pressure (P_sat)
figure;
plot(t, P_sat, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Saturation Pressure P_{sat} (MPa)', 'FontSize', 12);
title('Saturation Pressure over Time', 'FontSize', 14);
grid on;

%% Plot Thermal Power (P_th)
P_th = params.P_0 * phi;  % Calculate P_th based on neutron flux
figure;
plot(t, P_th, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Thermal Power P_{th} (W)', 'FontSize', 12);
title('Thermal Power over Time', 'FontSize', 14);
grid on;

%% Plot Average Core Temperature (T_avg_core)
T_avg_core = (T_c1 + T_c2) / 2;  % Calculate T_avg_core as average of coolant temperatures at nodes
figure;
plot(t, T_avg_core, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Average Core Temperature T_{avg\_core} (°C)', 'FontSize', 12);
title('Average Core Temperature over Time', 'FontSize', 14);
grid on;

%% Plot Primary Coolant Mass Flow Rate (m_cp)
m_cp = params.m_cp0 * phi.^(1/3);  % Calculate m_cp based on neutron flux
figure;
plot(t, m_cp, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Primary Coolant Mass Flow Rate m_{cp} (kg/s)', 'FontSize', 12);
title('Primary Coolant Mass Flow Rate over Time', 'FontSize', 14);
grid on;

%% Plot Neutron Flux (phi)
figure;
plot(t, phi, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Neutron Flux \phi', 'FontSize', 12);
title('Neutron Flux over Time', 'FontSize', 14);
grid on;

%% Plot Saturation Temperature (T_sat)
T_sat = zeros(length(t), 1);
for i = 1:length(t)
    % Convert P_sat to bar for XSteam function
    P_sat_bar = P_sat(i) * 10;
    T_sat(i) = XSteam('Tsat_p', P_sat_bar);  % Saturation temperature in °C
end
figure;
plot(t, T_sat, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Saturation Temperature T_{sat} (°C)', 'FontSize', 12);
title('Saturation Temperature over Time', 'FontSize', 14);
grid on;

%% Plot (U_v - c_pi * T_fi)
U_v_minus_c_pi_T_fi = zeros(length(t), 1);
for i = 1:length(t)
    % Get saturation temperature T_sat for current P_sat
    P_sat_bar = P_sat(i) * 10;
    T_sat_i = XSteam('Tsat_p', P_sat_bar);

    % Calculate internal energy of saturated vapor at T_sat (U_v)
    U_v = XSteam('uV_T', T_sat_i) * 1e3;  % Convert from kJ/kg to J/kg

    % Calculate internal energy of saturated Liquid at T_sat (U_w)
    U_w = XSteam('uL_T', T_sat_i) * 1e3;  % Convert from kJ/kg to J/kg

    % Specific heat capacity of feedwater at P_sat and T_fi (c_pi)
    c_pi = XSteam('Cp_pT', P_sat_bar, T_fi) * 1e3;  % J/kg·°C

    % Calculate (U_v - c_pi * T_fi)
    U_v_minus_c_pi_T_fi(i) = U_v - (c_pi * T_fi);
end
figure;
plot(t, U_v_minus_c_pi_T_fi, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('(U_v - c_{pi} * T_{fi}) (J/kg)', 'FontSize', 12);
title('Difference (U_v - c_{pi} * T_{fi}) over Time', 'FontSize', 14);
grid on;

%% Plot Primary Coolant Temperature in SG (T_p)
figure;
plot(t, T_p, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Primary Coolant Temperature T_p (°C)', 'FontSize', 12);
title('Primary Coolant Temperature in SG over Time', 'FontSize', 14);
grid on;

%% Plot Coolant Temperature at Node 1 (T_c1)
figure;
plot(t, T_c1, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Coolant Temperature at Node 1, T_{c1} (°C)', 'FontSize', 12);
title('Coolant Temperature at Node 1 over Time', 'FontSize', 14);
grid on;

%% Plot Coolant Temperature at Node 2 (T_c2)
figure;
plot(t, T_c2, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Coolant Temperature at Node 2, T_{c2} (°C)', 'FontSize', 12);
title('Coolant Temperature at Node 2 over Time', 'FontSize', 14);
grid on;

%% Plot Steam Generator Metal Temperature (T_m)
figure;
plot(t, T_m, 'LineWidth', 2);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Steam Generator Metal Temperature T_m (°C)', 'FontSize', 12);
title('Steam Generator Metal Temperature over Time', 'FontSize', 14);
grid on;

%% Function defining the differential equations
function dydt = reactor_ode(t, y, params)
    % Unpack state variables
    phi = y(1);      % Neutron flux
    C = y(2);        % Delayed neutron precursor concentration
    T_f = y(3);      % Fuel temperature
    T_c1 = y(4);     % Coolant temperature at Node 1
    T_c2 = y(5);     % Coolant temperature at Node 2
    T_HL = y(6);     % Hot leg temperature
    T_CL = y(7);     % Cold leg temperature
    T_p = y(8);      % Primary coolant temperature in SG
    T_m = y(9);      % SG metal temperature
    P_sat = y(10);   % Saturation pressure in SG secondary (MPa)

    % Unpack parameters
    alpha_f = params.alpha_f;
    alpha_c = params.alpha_c;
    beta = params.beta;
    Lambda = params.Lambda;
    lambda_p = params.lambda;
    P_0 = params.P_0;
    tau = params.tau;
    h_fc = params.h_fc;
    A_fc = params.A_fc;
    m_f = params.m_f;
    c_pf = params.c_pf;
    V_core = params.V_core;
    P_primary_bar = params.P_primary_bar;
    T_f0 = params.T_f0;
    T_c10 = params.T_c10;
    T_c20 = params.T_c20;
    m_cp0 = params.m_cp0;
    h_pm = params.h_pm;
    A_pm = params.A_pm;
    m_m = params.m_m;
    c_m = params.c_m;
    h_ms = params.h_ms;
    A_ms = params.A_ms;
    m_cs0 = params.m_cs0;
    T_fi = params.T_fi;
    P_sat0 = params.P_sat0;
    V_HL = params.V_HL;
    V_CL = params.V_CL;
    V_SG_primary = params.V_SG_primary;
    K_s = params.K_s;
    T_sat0 = params.T_sat0;
    K_sm = params.K_sm;

    % External reactivity step at t = 200 s
    if t < 200
        rho_ext = 0;
    else
        rho_ext = 3.5e-4;
    end

    % Reactivity feedback
    rho = rho_ext + alpha_f * (T_f - T_f0) + alpha_c * ((T_c1 - T_c10) + (T_c2 - T_c20)) / 2;

    % Neutron kinetics equations
    dphi_dt = ((rho - beta) / Lambda) * phi + lambda_p * C;
    dC_dt = (beta / Lambda) * phi - lambda_p * C;

    % Thermal power
    P_th = P_0 * phi;

    % Mass flow rate of primary coolant
    m_cp = m_cp0 * phi^(1/3);

    % Average core temperature
    T_avg_core = (T_c1 + T_c2) / 2;

    % Specific heat capacity and specific volume (use XSteam)
    c_pc = XSteam('Cp_pT', P_primary_bar, T_avg_core);  % kJ/kg·°C
    c_pc = c_pc * 1e3;  % Convert to J/kg·°C
    v_core = XSteam('v_pT', P_primary_bar, T_avg_core);  % m³/kg

    % Mass of primary coolant in core region
    m_c = V_core / v_core;  % kg

    % Fuel temperature equation
    dT_f_dt = (tau * P_th - h_fc * A_fc * (T_f - T_c1)) / (m_f * c_pf);
    %dT_f_dt = (tau * P_th - h_fc * 2000 * (T_f - T_c1)) / (m_f * 497);


    % Primary coolant temperature at Node 1
    dT_c1_dt = ((1 - tau) * P_th + h_fc * A_fc * (T_f - T_c1)) / (m_c * c_pc) ...
        + (2 * m_cp * (T_CL - T_c1)) / m_c;

    % Primary coolant temperature at Node 2
    dT_c2_dt = ((1 - tau) * P_th + h_fc * A_fc * (T_f - T_c2)) / (m_c * c_pc) ...
        + (2 * m_cp * (T_c1 - T_c2)) / m_c;

    % Time constants
    % Masses in hot leg and cold leg regions
    v_HL = XSteam('v_pT', P_primary_bar, T_HL);
    m_HL = V_HL / v_HL;
    v_CL = XSteam('v_pT', P_primary_bar, T_CL);
    m_CL = V_CL / v_CL;

    tau_HL = m_HL / m_cp;
    tau_CL = m_CL / m_cp;

    % Hot leg temperature
    dT_HL_dt = (T_c2 - T_HL) / tau_HL;

    % Cold leg temperature
    dT_CL_dt = (2 * T_p - T_HL - T_CL) / tau_CL;

    % Primary coolant in SG region
    v_p = XSteam('v_pT', P_primary_bar, T_p);
    m_p = V_SG_primary / v_p;
    c_p = XSteam('Cp_pT', P_primary_bar, T_p);  % kJ/kg·°C
    c_p = c_p * 1e3;  % Convert to J/kg·°C

    % Heat transfer coefficients
    K_HL = m_cp / m_p;
    K_pm = (h_pm * A_pm) / (m_p * c_p);

    % SG metal lump equations
    K_mp = (h_pm * A_pm) / (m_m * c_m);
    K_ms = (h_ms * A_ms) / (m_m * c_m);

    % Update T_sat based on current P_sat
    P_sat_bar = P_sat * 10;  % Convert MPa to bar for XSteam
    T_sat = XSteam('Tsat_p', P_sat_bar);

    % Internal energy of saturated vapor at T_sat
    U_v = XSteam('uV_T', T_sat) * 1e3;  % Convert from kJ/kg to J/kg

    % Internal energy of saturated Liquid at T_sat
    U_w = XSteam('uL_T', T_sat) * 1e3;  % Convert from kJ/kg to J/kg

    % Specific heat capacity of feedwater at P_sat and T_fi
    c_pi = XSteam('Cp_pT', P_sat_bar, T_fi) * 1e3;  % J/kg·°C

    % Primary coolant temperature in SG region
    dT_p_dt = K_HL * (T_HL - T_p) + (K_pm * (T_m - T_p))/2;

    % SG metal lump temperature
    dT_m_dt = K_mp * (T_p - T_m) + K_ms * (T_sat - T_m);

    % Saturation pressure dynamics
    dP_sat_dt = (K_sm * (T_m - T_sat) - m_cs0 * (U_v - c_pi * T_fi)) / K_s;

    % Pack derivatives
    dydt = zeros(10,1);
    dydt(1) = dphi_dt;
    dydt(2) = dC_dt;
    dydt(3) = dT_f_dt;
    dydt(4) = dT_c1_dt;
    dydt(5) = dT_c2_dt;
    dydt(6) = dT_HL_dt;
    dydt(7) = dT_CL_dt;
    dydt(8) = dT_p_dt;
    dydt(9) = dT_m_dt;
    dydt(10) = dP_sat_dt;  % P_sat in MPa
end
