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
clear,clc,clf,close all

    
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
% reentry heat and velocity + try do decrease speed : rocket +
% parachute/other (take into account the m_s added, compared to the burned
% mp) to help braking
% approx 100m ground : hoovering with high thrusts to reduce speed to 0
% when touching
Re=6378e3; % mean radius of Earth (m)
Cd = 0.85; %Drag coefficient. 1st assumption: the rocket is a cylinder (cf. Wikipedia Drag Coefficient)
A = 1; %Surface of the rocket in contact with the airflow (m^2)
stage=1;
T = [11e4, 11e4, 180e3]; %stages' thrust (N)

Isp = 266; %
Initial_speed= 6000e3/3600 ;% for example 8000km/h ie, will be y(1) as soon as possible
Initial_gamma= pi/4; %Assuming a first stage separation at around pi/4, will be y(2) as soon as possible.
Initial_x = 10000; % Assuming a first stage separation at around 10km, will be y(3) as soon as possible.
Initial_h= 60e3; % Assuming a first stage separation at around 100km, will be y(4) as soon as possible.
Initial_m = 25000; % 20000kg structural mass + 5000kg remaining propellant mass, will be y(5)- M_other_stages - M_payload
Initial_alpha=0;

options = odeset('RelTol',1e-10,'AbsTol',1e-9);

%7.1 - turning phase : thursting to turn to opposite gamma
end_gamma=pi+5*pi/180-Initial_alpha;%TBD, 5° entry for the moment
phase=7.1;
param = [Isp, Cd, A, stage, phase];

y0_reentry_1 = [Initial_speed Initial_gamma Initial_h Initial_x Initial_m Initial_alpha]; % Initial state vector, 1 line

ti_reentry_1 = 600; % Initial time for phase 7.1 (s). TBD, will be end of burnout of phase 1
tf_reentry_1 = 10+ti_reentry_1; % Final time for phase 7.1 (s). TBD


[t_reentry_1, y_reentry_1] = ode45(@(t, y) reentry_dynamicsODE(t, 0, y, param, [Initial_alpha, end_gamma], [ti_reentry_1, tf_reentry_1]), (ti_reentry_1:0.05:tf_reentry_1), y0_reentry_1, options);

%7.2 going in the opposite direction
phase=7.2;
param = [Isp, Cd, A, stage, phase];

y0_reentry_2 = [y_reentry_1(end, 1) y_reentry_1(end,2) y_reentry_1(end,3) y_reentry_1(end,4) y_reentry_1(end,5) y_reentry_1(end,6)]; % Initial state vector, 1 line

ti_reentry_2 = tf_reentry_1; % Initial time for phase 7.2, end of phase 7.1 (s)
tf_reentry_2 = 60+ti_reentry_2; % Final time for phase 7.2 (s). TBD, for the moment 1 minute

[t_reentry_2, y_reentry_2] = ode45(@(t, y) reentry_dynamicsODE(t, T(stage), y, param), (ti_reentry_2:0.05:tf_reentry_2), y0_reentry_2, options);

%7.3 turning to gamma =pi/2
phase=7.3;
%end_gamma=3*pi/2;
end_gamma=pi/2-y_reentry_2(end,2);
param = [Isp, Cd, A, stage, phase];

y0_reentry_3 = [y_reentry_2(end, 1) y_reentry_2(end,2) y_reentry_2(end,3) y_reentry_2(end,4) y_reentry_2(end,5) y_reentry_2(end,6)]; % Initial state vector, 1 line

ti_reentry_3 = tf_reentry_2; % Initial time for phase 7.3, end of phase 7.2 (s)
tf_reentry_3 = 60+ti_reentry_3; % Final time for phase 7.3 (s). TBD, for the moment 10 s under : [y0_reentry_3(2), end_gamma]

[t_reentry_3, y_reentry_3] = ode45(@(t, y) reentry_dynamicsODE(t, T(stage), y, param, [y0_reentry_3(6), end_gamma], [ti_reentry_3, tf_reentry_3]), (ti_reentry_3:0.05:tf_reentry_3), y0_reentry_3, options);

%7.4 - braking phase: parachute to brake : to be deployed after 340 m/s
%only
phase=7.4;
Cd = 1.5; %Drag coefficient for a space parachute (kevlar/nylon stuff, deployed at approx 10km to avoid it to brake)
A = 20; %Surface of the 3 parachutes in contact with the airflow (m^2) TBD
param = [Isp, Cd, A, stage, phase];

y0_reentry_4 = [y_reentry_3(end, 1) y_reentry_3(end,2) y_reentry_3(end,3) y_reentry_3(end,4) y_reentry_3(end,5) y_reentry_3(end,6)]; % Initial state vector, 1 line

ti_reentry_4 = tf_reentry_3; % Initial time for phase 7.4, end of phase 7.3 (s)
tf_reentry_4 = 60*4+ti_reentry_4; % Final time for phase 7.4 (s). TBD, for the moment 3 minutes


[t_reentry_4, y_reentry_4] = ode45(@(t, y) reentry_dynamicsODE(t, 0 , y, param), (ti_reentry_4:0.05:tf_reentry_4), y0_reentry_4, options);

%7.5 - controlled descent phase : landing (high thrust). Maybe remove
%parachute ?
phase=7.5;

y0_reentry_5 = [y_reentry_4(end, 1) y_reentry_4(end,2) y_reentry_4(end,3) y_reentry_4(end,4) y_reentry_4(end,5) y_reentry_4(end,6)]; % Initial state vector, 1 line

ti_reentry_5 = tf_reentry_4; % Initial time for phase 7.5, end of phase 7.4(s)
tf_reentry_5 = 20+ti_reentry_5; % Final time for phase 7.5 (s). TBD, for the moment 20s


[t_reentry_5, y_reentry_5] = ode45(@(t, y) reentry_dynamicsODE(t, -T(stage), y, param), (ti_reentry_5:0.05:tf_reentry_5), y0_reentry_5, options);



%% Ploting phase
figure(1); hold on;
plot(t_reentry_1,y_reentry_1(:,3)/1e3,'r','LineWidth',2);
plot(t_reentry_2,y_reentry_2(:,3)/1e3,'g','LineWidth',2);
plot(t_reentry_3,y_reentry_3(:,3)/1e3,'b','LineWidth',2);
plot(t_reentry_4,y_reentry_4(:,3)/1e3,'c','LineWidth',2);
plot(t_reentry_5,y_reentry_5(:,3)/1e3,'m','LineWidth',2);
title('Altitude change');
xlabel('Time (s)');
ylabel('Altitude (km)');
grid;

figure(2); hold on;
plot(t_reentry_1,y_reentry_1(:,2)*180/pi,'r','LineWidth',2);
plot(t_reentry_2,y_reentry_2(:,2)*180/pi,'g','LineWidth',2);
plot(t_reentry_3,y_reentry_3(:,2)*180/pi,'b','LineWidth',2);
plot(t_reentry_4,y_reentry_4(:,2)*180/pi,'c','LineWidth',2);
plot(t_reentry_5,y_reentry_5(:,2)*180/pi,'m','LineWidth',2);
title('Flight path angle change');
xlabel('Time (s)');
ylabel('Flight path angle (deg)');
grid;

figure(3); hold on;
plot(y_reentry_1(:,4)/1e3,y_reentry_1(:,3)/1e3,'r','LineWidth',2);
plot(y_reentry_2(:,4)/1e3,y_reentry_2(:,3)/1e3,'g','LineWidth',2);
plot(y_reentry_3(:,4)/1e3,y_reentry_3(:,3)/1e3,'b','LineWidth',2);
plot(y_reentry_4(:,4)/1e3,y_reentry_4(:,3)/1e3,'c','LineWidth',2);
plot(y_reentry_5(:,4)/1e3,y_reentry_5(:,3)/1e3,'m','LineWidth',2);
title('Altitude change');
xlabel('X position (km)');
ylabel('Altitude (km)');
grid;

figure(4); hold on;
plot(t_reentry_1,y_reentry_1(:,5),'r','LineWidth',2);
plot(t_reentry_2,y_reentry_2(:,5),'g','LineWidth',2);
plot(t_reentry_3,y_reentry_3(:,5),'b','LineWidth',2);
plot(t_reentry_4,y_reentry_4(:,5),'c','LineWidth',2);
plot(t_reentry_5,y_reentry_5(:,5),'m','LineWidth',2);
title('Mass change');
xlabel('Time (s)');
ylabel('Mass (kg)');
grid;

figure(5); hold on;
plot(t_reentry_1,y_reentry_1(:,1),'r','LineWidth',2);
plot(t_reentry_2,y_reentry_2(:,1),'g','LineWidth',2);
plot(t_reentry_3,y_reentry_3(:,1),'b','LineWidth',2);
plot(t_reentry_4,y_reentry_4(:,1),'c','LineWidth',2);
plot(t_reentry_5,y_reentry_5(:,1),'m','LineWidth',2);
title('velocity change');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid;

figure(6); hold on;
plot(t_reentry_1,y_reentry_1(:,6)*180/pi,'r','LineWidth',2);
plot(t_reentry_2,y_reentry_2(:,6)*180/pi,'g','LineWidth',2);
plot(t_reentry_3,y_reentry_3(:,6)*180/pi,'b','LineWidth',2);
plot(t_reentry_4,y_reentry_4(:,6)*180/pi,'c','LineWidth',2);
plot(t_reentry_5,y_reentry_5(:,6)*180/pi,'m','LineWidth',2);
title('alpha angle change');
xlabel('Time (s)');
ylabel('alpha angle (deg)');
grid;
