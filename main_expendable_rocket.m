%% Description
% Gravity turn (cf. SD2900 - Rocket Dynamics, slide 25)
    %5 phases:
    % 1. Let the 1st stage go vertical for a short time to pick up speed to
    % avoid a too large dgamma/dt
    % 2. Set the flight path angle to a small non-zero value to start the
    % gravity turn
    % 3. At the burnout of the 1st stage, stop the time integration and restart
    % the integration with data for 2nd stage
    % 4. At the burnout of 2nd stage, stop the time integration and restart the
    % integration with data for 3rd stage.
    % 5. Use the tangent steering law for the 3rd stage to enforce gamma=0 at
    % the correct altitude. Note that gamma is NOT a variable for 3rd stage is
    % used.
    % Final conditions: gamma=0, H=H*, V=V*
    % 6. Once we are at this parking orbit, we begin the Hohmann transfer to reach the target orbit.   

 %% Simulating the powered descent landing on the Moon
clear,clc,clf

    
%% Constant values

mu_E=3.968e14; % gravitational parameter of Earth (m^3s^-2)
Re=6378e3; % mean radius of Earth (m)
g0=9.80665; % gravity of Earth at sea level (m/s^2)

m0 = 1000e3; % initiale mass of the rocket (kg)
t0 = 0; % ignition time of the 1st stage (s)

%% Phases

% 1.
% Initial condition*

% On a paper we got 12s.

%We need to estimate the thrust delivered by the engines
%[t, y] = ode45(@(t, y) ascent_dynamicsODE(T, y0), [t0, tf1], y0);


% 6. Hohmann transfer : using 3rd stage
Re=6378e3; % mean radius of Earth (m)
h= Re + 200e3; % Assuming a parking orbit of 200km.
mu = 3.986004418e14; % constant, assuming M+m is approximately M and constant
h_target = 29599.8e3; % target orbit
m_init=2100; %initial mass for phase 6. 2100kg for the moment, will be y(5) in the end
Isp= 467; %Isp stage 3
m_star=732.8;

g0=9.80665; 

Vc1 = sqrt(mu/h);% current speed on parking orbit
Vc2 = sqrt(mu/h_target); % to be speed on MEO
a = (h+h_target)/2; %semi major axis of the transfer orbit
V_perigee = sqrt(mu*(2/h-1/a)); %velocity at the perigee of the transfer orbit 
V_apogee = sqrt(mu*(2/h_target-1/a)); %velocity at the apogee of the transfer orbit
Delta_V_perigee = V_perigee - Vc1;
Delta_V_apogee = Vc2 - V_apogee;
Delta_V = Delta_V_perigee + Delta_V_apogee; % cost of the total transfer 
Delta_t = pi*sqrt(a^3/h_target); % transfer time

g_perigee=mu/(h^2);
m_final_perigee=m_init*exp(-Delta_V_perigee/(Isp*g_perigee)); % final mass after first thrust (at perigee)
g_apogee=mu/(h_target^2);
m_final_apogee=m_final_perigee*exp(-Delta_V_apogee/(Isp*g_apogee)); % final mass after second thrust (at apogee)
Delta_m= m_init-m_final_apogee; %total propellant cost of the transfer


