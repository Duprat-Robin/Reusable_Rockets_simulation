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
Delta_t = pi*sqrt(a^3/mu); % transfer time

Delta_m = m_init-m_init*exp(-Delta_V/(Isp*g0));
m_s=Delta_m-m_star; %in general m_s = 1/7 m_p

% 7. reentry of a stage 
% reentry heat and velocity + try do decrease speed : rocket + maybe use a
% parachute/other (take into account the m_s added, compared to the burned
% mp) to help braking
% approx 100m ground : hoovering with high thrusts to reduce speed to 0
% when touching V gamma x h m
Re=6378e3; % mean radius of Earth (m)
Cd = 0.85; %Drag coefficient. 1st assumption: the rocket is a cylinder (cf. Wikipedia Drag Coefficient)
A = 1; %Surface of the rocket in contact with the airflow (m^2)
stage=1;
T = [11e4, 11e4, 180e3]; %stages' thrust (N)

Isp = 266; %
Initial_speed= 6000e3/3600 ;% for example 8000km/h ie, will be y(1) as soon as possible
Initial_gamma= pi/4; %Assuming a first stage separation at around pi/4, will be y(2) as soon as possible.
Initial_x = 10000; % Assuming a first stage separation at around 10km, will be y(3) as soon as possible.
Initial_h= 100e3; % Assuming a first stage separation at around 100km, will be y(4) as soon as possible.
Initial_m = 25000; % 20000kg structural mass + 5000kg remaining propellant mass, will be y(5)- M_other_stages - M_payload

options = odeset('RelTol',1e-10,'AbsTol',1e-9);


