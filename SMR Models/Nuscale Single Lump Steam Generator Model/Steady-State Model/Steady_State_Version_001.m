% MATLAB Script to Solve Steady-State Equations of a Reactor System
% Author: [Your Name]
% Date: [Current Date]

% Clear workspace and command window
clear; clc;

% Add XSteam path if necessary
% addpath('path_to_XSteam');  % Uncomment and set the correct path if needed

%% Known Parameters (Provided Values)

% Neutron kinetics parameters
beta = 0.007;            % Total delayed neutron fraction (dimensionless)
Lambda = 2e-5;           % Neutron generation time (s)
lambda = 0.1;            % Decay constant of delayed neutron precursors (1/s)
alpha_f = -2.16e-5;      % Fuel temperature coefficient of reactivity (1/°C)
alpha_c = -1.8e-4;       % Coolant temperature coefficient of reactivity (1/°C)
rho_ext = 0;             % External reactivity insertion (dimensionless)

% Reactor and primary system parameters
P0 = 160e6;              % Rated thermal power (W)
tau = 0.97;              % Fraction of power deposited in fuel (dimensionless)
h_fc = 1135;             % Heat transfer coefficient between fuel and coolant (W/m²·°C)
A_fc = 583;              % Heat transfer area between fuel and coolant (m²)
m_f = 11252;             % Mass of fuel (kg)
c_pf = 467;              % Specific heat capacity of fuel (J/kg·°C)
h_pm = 20391;            % Heat transfer coefficient between primary coolant and SG metal (W/m²·°C)
A_pm = 1123;             % Heat transfer area between primary coolant and SG metal (m²)
h_ms = 4950;             % Heat transfer coefficient between SG metal and secondary coolant (W/m²·°C)
A_ms = 1214;             % Heat transfer area between SG metal and secondary coolant (m²)
m_m = 7869;              % Mass of SG metal lump (kg)
c_m = 450;               % Specific heat capacity of SG metal lump (J/kg·°C)
T_fi = 170;              % Feedwater inlet temperature (°C)
mdot_cp_rated = 587;     % Rated mass flow rate of primary coolant (kg/s)
V_core = 1.879;          % Volume of core region (m³)
V_SG_primary = 3.564;    % Volume of primary coolant in SG region (m³)
V_HL = 9.7;              % Volume of hot leg (m³)
V_CL = 26.8;             % Volume of cold leg (m³)
P_primary = 127.6;       % Primary system pressure (bar)
eta_T = 0.6;             % Turbine efficiency (dimensionless)
P_mech_fixed = 45e6;     % Mechanical power output (W)
T_p_fixed = 252.5;         % Fixed primary coolant temperature in SG region (°C)

% Thermodynamic constants
% No need to define c_pi here since it depends on T_fi and will be calculated using XSteam

%% Initial Guesses for Unknown Variables

% Initial guesses for unknown temperatures (°C)
T_f_guess = 600;         % Fuel temperature
T_c1_guess = 300;        % Primary coolant temperature at node 1
T_c2_guess = 320;        % Primary coolant temperature at node 2
T_m_guess = 250;         % SG metal lump temperature
T_sat_guess = 240;       % Saturation temperature in SG
mdot_cs_guess = 100;     % Mass flow rate of secondary coolant (kg/s)
phi_guess = 0.7;         % Neutron flux (dimensionless)

% Combine initial guesses into a vector
x0 = [T_f_guess; T_c1_guess; T_c2_guess; T_m_guess; T_sat_guess; mdot_cs_guess; phi_guess];

%% Solve the System of Equations Using fsolve

% Define options for fsolve
options = optimoptions('fsolve', 'Display', 'iter', 'TolFun', 1e-6);

% Call fsolve to solve the system of equations
[x_sol, fval, exitflag, output] = fsolve(@steady_state_equations, x0, options);

% Extract the solutions
T_f = x_sol(1);         % Fuel temperature (°C)
T_c1 = x_sol(2);        % Primary coolant temperature at node 1 (°C)
T_c2 = x_sol(3);        % Primary coolant temperature at node 2 (°C)
T_m = x_sol(4);         % SG metal lump temperature (°C)
T_sat = x_sol(5);       % Saturation temperature in SG (°C)
mdot_cs = x_sol(6);     % Mass flow rate of secondary coolant (kg/s)
phi = x_sol(7);         % Neutron flux (dimensionless)

% Calculate dependent variables
P_th = P0 * phi;        % Thermal power (W)
mdot_cp = mdot_cp_rated * phi^(1/3);   % Mass flow rate of primary coolant (kg/s)

% Additional calculations for missing parameters

% Hot Leg Temperature
T_HL = T_c2;                             % Hot leg temperature (°C)

% Cold Leg Temperature
T_p = T_p_fixed;                         % Fixed primary coolant temperature in SG region (°C)
T_CL = 2 * T_p - T_HL;                   % Cold leg temperature (°C)

% Specific heat capacity of primary coolant in core region (J/kg·°C)
T_avg_core = (T_c1 + T_c2) / 2;          % Average temperature in core (°C)
c_pc = XSteam('Cp_pT', P_primary, T_avg_core) * 1000;  % Convert kJ/kg·°C to J/kg·°C

% Specific volume of primary coolant in SG region (m³/kg)
v_p = XSteam('v_pT', P_primary, T_p);    % m³/kg

% Mass of primary coolant in SG region (kg)
m_p = V_SG_primary / v_p;

% Heat transfer coefficients
K_HL = mdot_cp / m_p;                                    % (1/s)
K_m = (h_pm * A_pm) / (m_p * c_pc);                      % (1/s)
K_mp = (h_pm * A_pm) / (m_m * c_m);                      % (1/s)
K_ms = (h_ms * A_ms) / (m_m * c_m);                      % (1/s)
K_sm = h_ms * A_ms;                                      % (W/°C)

% Specific heat capacity of feedwater (J/kg·°C)
c_pi = XSteam('CpL_T', T_fi) * 1000;                     % Convert kJ/kg·°C to J/kg·°C

% Internal energy of saturated vapor at saturation temperature (kJ/kg)
U_v = XSteam('uV_T', T_sat) * 1000;                             % kJ/kg

% Enthalpy difference (latent heat of vaporization) at saturation temperature (kJ/kg)
h_fg = XSteam('hV_T', T_sat) - XSteam('hL_T', T_sat);    % kJ/kg

% Saturation pressure corresponding to T_sat (bar)
P_sat = XSteam('psat_T', T_sat);                         % Saturation pressure (bar)

% Mass of primary coolant in core region (kg)
v_core = XSteam('v_pT', P_primary, T_avg_core);          % Specific volume in core region (m³/kg)
m_c = V_core / v_core;                                   % Mass in core region (kg)


%% Display results
fprintf('\nSteady-State Solution:\n');
fprintf('-----------------------\n');

% Known Parameters
fprintf('\nKnown Parameters:\n');
fprintf('----------------\n');
fprintf('beta (Total delayed neutron fraction): %.4f\n', beta);
fprintf('Lambda (Neutron generation time): %.2e s\n', Lambda);
fprintf('lambda (Decay constant of delayed neutron precursors): %.2f 1/s\n', lambda);
fprintf('alpha_f (Fuel temperature coefficient of reactivity): %.2e 1/°C\n', alpha_f);
fprintf('alpha_c (Coolant temperature coefficient of reactivity): %.2e 1/°C\n', alpha_c);
fprintf('rho_ext (External reactivity insertion): %.2f\n', rho_ext);
fprintf('P0 (Rated thermal power): %.2f MW\n', P0/1e6);
fprintf('tau (Fraction of power deposited in fuel): %.2f\n', tau);
fprintf('h_fc (Heat transfer coefficient between fuel and coolant): %.2f W/m²·°C\n', h_fc);
fprintf('A_fc (Heat transfer area between fuel and coolant): %.2f m²\n', A_fc);
fprintf('m_f (Mass of fuel): %.2f kg\n', m_f);
fprintf('c_pf (Specific heat capacity of fuel): %.2f J/kg·°C\n', c_pf);
fprintf('h_pm (Heat transfer coefficient between primary coolant and SG metal): %.2f W/m²·°C\n', h_pm);
fprintf('A_pm (Heat transfer area between primary coolant and SG metal): %.2f m²\n', A_pm);
fprintf('h_ms (Heat transfer coefficient between SG metal and secondary coolant): %.2f W/m²·°C\n', h_ms);
fprintf('A_ms (Heat transfer area between SG metal and secondary coolant): %.2f m²\n', A_ms);
fprintf('m_m (Mass of SG metal lump): %.2f kg\n', m_m);
fprintf('c_m (Specific heat capacity of SG metal lump): %.2f J/kg·°C\n', c_m);
fprintf('T_fi (Feedwater inlet temperature): %.2f °C\n', T_fi);
fprintf('mdot_cp_rated (Rated mass flow rate of primary coolant): %.2f kg/s\n', mdot_cp_rated);
fprintf('V_core (Volume of core region): %.3f m³\n', V_core);
fprintf('V_SG_primary (Volume of primary coolant in SG region): %.3f m³\n', V_SG_primary);
fprintf('V_HL (Volume of hot leg): %.3f m³\n', V_HL);
fprintf('V_CL (Volume of cold leg): %.3f m³\n', V_CL);
fprintf('P_primary (Primary system pressure): %.2f bar\n', P_primary);
fprintf('eta_T (Turbine efficiency): %.2f\n', eta_T);
fprintf('P_mech_fixed (Mechanical power output): %.2f MW\n', P_mech_fixed/1e6);
fprintf('T_p_fixed (Fixed primary coolant temperature in SG region): %.2f °C\n', T_p_fixed);

% Unknown Parameters (Calculated Results)
fprintf('\nUnknown Parameters (Calculated Results):\n');
fprintf('----------------------------------------\n');
fprintf('Neutron Flux (phi): %.4f\n', phi);
fprintf('Thermal Power (P_th): %.2f MW\n', P_th / 1e6);
fprintf('Fuel Temperature (T_f): %.2f °C\n', T_f);
fprintf('Primary Coolant Temperature at Node 1 (T_c1): %.2f °C\n', T_c1);
fprintf('Primary Coolant Temperature at Node 2 (T_c2): %.2f °C\n', T_c2);
fprintf('Hot Leg Temperature (T_HL): %.2f °C\n', T_HL);
fprintf('Cold Leg Temperature (T_CL): %.2f °C\n', T_CL);
fprintf('SG Metal Lump Temperature (T_m): %.2f °C\n', T_m);
fprintf('Saturation Temperature in SG (T_sat): %.2f °C\n', T_sat);
fprintf('Saturation Pressure in SG (P_sat): %.2f bar\n', P_sat);
fprintf('Mass Flow Rate of Primary Coolant (mdot_cp): %.2f kg/s\n', mdot_cp);
fprintf('Mass Flow Rate of Secondary Coolant (mdot_cs): %.2f kg/s\n', mdot_cs);
fprintf('Specific Heat Capacity of Primary Coolant in Core Region (c_pc): %.2f J/kg·°C\n', c_pc);
fprintf('Specific Volume of Primary Coolant in SG (v_p): %.6f m³/kg\n', v_p);
fprintf('Mass of Primary Coolant in SG Region (m_p): %.2f kg\n', m_p);
fprintf('Internal Energy of Saturated Vapor (U_v): %.2f kJ/kg\n', U_v);
fprintf('Enthalpy Difference (h_fg): %.2f kJ/kg\n', h_fg);
fprintf('Heat Transfer Coefficient K_HL: %.6f 1/s\n', K_HL);
fprintf('Heat Transfer Coefficient K_m: %.6f 1/s\n', K_m);
fprintf('Heat Transfer Coefficient K_mp: %.6f 1/s\n', K_mp);
fprintf('Heat Transfer Coefficient K_ms: %.6f 1/s\n', K_ms);
fprintf('Heat Transfer Coefficient K_sm: %.2f W/°C\n', K_sm);
fprintf('Specific Heat Capacity of Feedwater (c_pi): %.2f J/kg·°C\n', c_pi);
fprintf('Mass of Primary Coolant in Core Region (m_c): %.2f kg\n', m_c);

%% Function Defining the System of Steady-State Equations

function F = steady_state_equations(x)
    % Unpack the unknown variables
    T_f = x(1);       % Fuel temperature (°C)
    T_c1 = x(2);      % Primary coolant temperature at node 1 (°C)
    T_c2 = x(3);      % Primary coolant temperature at node 2 (°C)
    T_m = x(4);       % SG metal lump temperature (°C)
    T_sat = x(5);     % Saturation temperature in SG (°C)
    mdot_cs = x(6);   % Mass flow rate of secondary coolant (kg/s)
    phi = x(7);       % Neutron flux (dimensionless)

    % Known parameters (from the main script)
    P0 = 160e6;              % Rated thermal power (W)
    tau = 0.97;              % Fraction of power deposited in fuel
    h_fc = 1135;             % Heat transfer coefficient between fuel and coolant (W/m²·°C)
    A_fc = 583;              % Heat transfer area between fuel and coolant (m²)
    h_pm = 20391;            % Heat transfer coefficient between primary coolant and SG metal (W/m²·°C)
    A_pm = 1123;             % Heat transfer area between primary coolant and SG metal (m²)
    h_ms = 4950;             % Heat transfer coefficient between SG metal and secondary coolant (W/m²·°C)
    A_ms = 1214;             % Heat transfer area between SG metal and secondary coolant (m²)
    m_m = 7869;              % Mass of SG metal lump (kg)
    c_m = 450;               % Specific heat capacity of SG metal lump (J/kg·°C)
    mdot_cp_rated = 587;     % Rated mass flow rate of primary coolant (kg/s)
    P_primary = 127.6;       % Primary system pressure (bar)
    T_p = 252.5;               % Fixed primary coolant temperature in SG region (°C)
    V_SG_primary = 3.564;    % Volume of primary coolant in SG region (m³)
    V_core = 1.879;          % Volume of core region (m³)
    T_fi = 170;              % Feedwater inlet temperature (°C)
    eta_T = 0.69;             % Turbine efficiency (dimensionless)
    P_mech = 45e6;           % Mechanical power output (W)
    c_pf = 467;              % Specific heat capacity of fuel (J/kg·°C)
    m_f = 11252;             % Mass of fuel (kg)
    V_HL = 9.7;              % Volume of hot leg (m³)
    V_CL = 26.8;             % Volume of cold leg (m³)
    m_m = 7869;              % Mass of SG metal lump (kg)
    c_m = 450;               % Specific heat capacity of SG metal lump (J/kg·°C)
    h_ms = 4950;             % Heat transfer coefficient between SG metal and secondary coolant (W/m²·°C)
    A_ms = 1214;             % Heat transfer area between SG metal and secondary coolant (m²)
    h_pm = 20391;            % Heat transfer coefficient between primary coolant and SG metal (W/m²·°C)
    A_pm = 1123;             % Heat transfer area between primary coolant and SG metal (m²)
    h_fc = 1135;             % Heat transfer coefficient between fuel and coolant (W/m²·°C)
    A_fc = 583;              % Heat transfer area between fuel and coolant (m²)

    % Mass flow rate of primary coolant (kg/s)
    mdot_cp = mdot_cp_rated * phi^(1/3);

    % Thermal power (W)
    P_th = P0 * phi;

    % Specific heat capacity of primary coolant in core region (J/kg·°C)
    T_avg_core = (T_c1 + T_c2) / 2;    % Average temperature in core (°C)
    c_pc = XSteam('Cp_pT', P_primary, T_c1) * 1000;  % Convert kJ/kg·°C to J/kg·°C
    %c_pc = 4900;

    % Specific volume of primary coolant in core region (m³/kg)
    v_core = XSteam('v_pT', P_primary, T_avg_core);  % m³/kg

    % Mass of primary coolant in core region (kg)
    m_c = V_core / v_core;

    % Specific heat capacity of primary coolant in SG region (J/kg·°C)
    c_p = XSteam('Cp_pT', P_primary, T_p) * 1000;  % Convert kJ/kg·°C to J/kg·°C
    %c_p = 4775;

    % Specific volume of primary coolant in SG region (m³/kg)
    v_p = XSteam('v_pT', P_primary, T_p);  % m³/kg

    % Mass of primary coolant in SG region (kg)
    m_p = V_SG_primary / v_p;

    % Heat transfer coefficients
    K_HL = mdot_cp / m_p;    % (1/s)
    K_m = (h_pm * A_pm) / (m_p * c_p);  % (1/s)
    K_mp = (h_pm * A_pm) / (m_m * c_m);  % (1/s)
    K_ms = (h_ms * A_ms) / (m_m * c_m);  % (1/s)
    K_sm = h_ms * A_ms;    % (W/°C)

    % Specific heat capacity of feedwater (J/kg·°C)
    %c_pi = XSteam('CpL_T', T_fi) * 1000;  % Convert kJ/kg·°C to J/kg·°C
    %c_pi = XSteam('Cv_pT',(XSteam('psat_T', T_sat)), T_fi) * 1000;  % Convert kJ/kg·°C to J/kg·°C
    c_pi = 4362;

    % Internal energy of saturated vapor at saturation temperature (kJ/kg)
    U_v = XSteam('uV_T', T_sat) * 1000;  % kJ/kg

    h_v = XSteam('hV_T', T_sat) * 1000;  % kJ/kg

    % Enthalpy difference (latent heat of vaporization) at saturation temperature (kJ/kg)
    S_Turbine = XSteam('sV_T', T_sat);
    h_fg = XSteam('hV_T', T_sat) - XSteam('h_ps', 0.08, S_Turbine);  % kJ/kg

    
    % Equation residuals
    F = zeros(8,1);

    % Equation 1: Fuel temperature equation
    F(1) = tau * P_th + h_fc * A_fc * (T_c1 - T_f);

    % Equation 2: Primary coolant temperature at node 1
    T_cl = 2 * T_p - T_c2;  % Cold leg temperature (°C)
    F(2) = (1 - tau) * P_th + h_fc * A_fc * (T_f - T_c1) + 2 * mdot_cp * c_pc * (T_cl - T_c1);

    % Equation 3: Primary coolant temperature at node 2
    F(3) = (1 - tau) * P_th + h_fc * A_fc * (T_f - T_c1) + 2 * mdot_cp * c_pc * (T_c1 - T_c2);

    % Equation 4: SG primary side heat transfer
    F(4) = K_HL * (T_c2 - T_p) + (K_m * (T_m - T_p))/2;

    % Equation 5: SG metal lump temperature
    F(5) = K_mp * (T_p - T_m) + K_ms * (T_sat - T_m);

    % Equation 6: Steady-state temperature difference relation (primary side)
    F(6) = (P0 * phi) / (2 * m_p * K_m * c_pc) - (T_p - T_m);

    % Equation 7: Steady-state temperature difference relation (secondary side)
    F(7) = (T_p - T_m) - (K_ms * mdot_cs * (h_v - c_pi * T_fi)) / (K_mp * K_sm); % replaced by 1e6 

    % Equation 8: Mechanical power equation
    F(8) = P_mech - eta_T * mdot_cs * h_fg * 1000;  % Convert kJ/kg to J/kg

    % Ensure that the number of equations matches the number of unknowns
end


