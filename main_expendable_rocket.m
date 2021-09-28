%% Description
% Gravity turn (cf. SD2900 - Rocket Dynamics, slide 25)
    %Assuming Lift = 0
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
    % Final conditions (for the rocket and not the mission): gamma=0, H=H*, V=V*

%% Constant values

mu_E=3.968e14; % gravitational parameter of Earth (m^3s^-2)
Re=6378e3; % mean radius of Earth (m)
g0=9.80665; % gravity of Earth at sea level (m/s^2)
Hf = 200e3; %parking orbit altitude (m)
Xf = 200e3; %distance from launch pad at the end of ascent phase (m). Usefull for BVP

t0 = 0; % ignition time of the 1st stage (s)

options = odeset('RelTol',1e-10,'AbsTol',1e-9);
iV = 1; igamma = 2; ih = 3; ix = 4; im = 5;
%% Rocket parameters
Isp = [378, 359, 467]; %Isp (s) for 1st stage. TBD
Cd = 0.85; %Drag coefficient. 1st assumption: the rocket is a cylinder (cf. Wikipedia Drag Coefficient)
A = pi*3.66^2; %Surface of the rocket in contact with the airflow (m^2)
nb_engines = [3, 1, 1];
T = [2205000, 533000, 180000].*nb_engines; %stages' trhrust (N)
ms = [16e3, 4e3, 1.5e3]; %stages' strucutal mass (kg)
mp = [418647.6, 31134, 3751.57]; %stages' propellant mass (kg)
m_p = 732.8; %Mass of the payload at launch (kg)
m0 = sum(ms) + sum(mp) + m_p; %Total mass of the rocket at lift-off (kg)

V0 = 0; % (m/s)
gamma0 = pi/2; % (rad)
gamma1 = gamma0 - 0.1*pi/180;
x0 = 0; % (m)
h0 = 0; % (m) Sea level, to adpat in function of the sarting point
%% Phases

% 1.
% Initial condition
phase = 1; %Current phase
stage = 1; %Stage currently use by the rocket
param = [Isp(stage), Cd, A, stage, phase];

y01 = [V0 gamma0 h0 x0 m0]; % Initial state vector, 1 line

ti1 = 10; % Final time for phase 1 (s). Must be a short time (i for intermediate)
% On a paper we got 12s.
[t1, y1] = ode45(@(t, y) ascent_dynamicsODE(t, T(stage), y, param, [gamma0, gamma1], [t0, ti1]), [t0 ti1], y01, options);

%2.
phase = 2;
param = [Isp(stage), Cd, A, stage, phase];

y02 = [y1(end,iV) y1(end,igamma) y1(end,ih) y1(end,ix) y1(end,im)];

tb1 = mp(stage)*g0*Isp(stage)/T(stage); %burnout time of 1st stage (s)
tf1 = tb1; %end of the 1st stage phase
[t2, y2] = ode45(@(t, y) ascent_dynamicsODE(t, T(stage), y, param), [ti1 tf1], y02, options);

%3.
phase = 3;
stage = 2;
param = [Isp(stage), Cd, A, stage, phase];

y03 = [y2(end,iV) y2(end,igamma) y2(end,ih) y2(end,ix) y2(end,im)];
y03(end,im) = y03(end,im)-ms(stage-1); %1st stage removal

tb2 = mp(stage)*g0*Isp(stage)/T(stage); %burnout time of 2nd stage (s)
tf2 = tb2 + tf1; %end of the 2nd stage phase
[t3, y3] = ode45(@(t, y) ascent_dynamicsODE(t, T(stage), y, param), [tf1 tf2], y03, options);

%4.
phase = 4;
stage = 3;
param = [Isp(stage), Cd, A, stage, phase];

y04 = [y3(end,iV) y3(end,igamma) y3(end,ih) y3(end,ix) y3(end,im)];
y04(end,im) = y04(end,im)-ms(stage-1); %2nd stage removal

tb3 = mp(stage)*g0*Isp(stage)/T(stage); %burnout time of 2nd stage (s)
tf3 = tb3 + tf2; %end of the 3rd stage phase
[t4, y4] = ode45(@(t, y) ascent_dynamicsODE(t, T(stage), y, param, [y3(end,igamma) 0], [tf2 tf3]), [tf2 tf3], y04, options);
for i=1:size(y4(:,2),1)
    y4(i,igamma) = atan(tan(y3(end,igamma)*(1-(t4(i)-tf2)/tb3))); %Steering law: linear tangent law
end

%% Ploting phase
figure(1); hold on;
plot(t1,y1(:,ih)/1e3,'r','LineWidth',2);
plot(t2,y2(:,ih)/1e3,'g','LineWidth',2);
plot(t3,y3(:,ih)/1e3,'b','LineWidth',2);
plot(t4,y4(:,ih)/1e3,'y','LineWidth',2);
title('Altitude change');
xlabel('Time (s)');
ylabel('Altitude (km)');
grid;

figure(2); hold on;
plot(t1,y1(:,igamma)*180/pi,'r','LineWidth',2);
plot(t2,y2(:,igamma)*180/pi,'g','LineWidth',2);
plot(t3,y3(:,igamma)*180/pi,'b','LineWidth',2);
plot(t4,y4(:,igamma)*180/pi,'y','LineWidth',2);
title('Flight path angle change');
xlabel('Time (s)');
ylabel('Flight path angle (deg)');
grid;

figure(3); hold on;
plot(y1(:,ix)/1e3,y1(:,ih)/1e3,'r','LineWidth',2);
plot(y2(:,ix)/1e3,y2(:,ih)/1e3,'g','LineWidth',2);
plot(y3(:,ix)/1e3,y3(:,ih)/1e3,'b','LineWidth',2);
plot(y4(:,ix)/1e3,y4(:,ih)/1e3,'y','LineWidth',2);
title('Altitude change');
xlabel('X position (km)');
ylabel('Altitude (km)');
grid;

figure(4); hold on;
plot(t1,y1(:,iV)/1e3,'r','LineWidth',2);
plot(t2,y2(:,iV)/1e3,'g','LineWidth',2);
plot(t3,y3(:,iV)/1e3,'b','LineWidth',2);
plot(t4,y4(:,iV)/1e3,'y','LineWidth',2);
title('Speed change');
xlabel('Time (s)');
ylabel('Speed (km/s)');
grid;