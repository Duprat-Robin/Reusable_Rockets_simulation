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

t0 = 0; % ignition time of the 1st stage (s)

options = odeset('RelTol',1e-10,'AbsTol',1e-9);
%% Rocket parameters
Isp = [400, 400, 467]; %Isp (s) for 1st stage. TBD
Cd = 0.85; %Drag coefficient. 1st assumption: the rocket is a cylinder (cf. Wikipedia Drag Coefficient)
A = 1; %Surface of the rocket in contact with the airflow (m^2)

T = [11e6, 11e6, 180e3]; %stages' trhrust (N)
ms = [100, 100, 100]; %stages' strucutal mass (kg)
mp = [500, 500, 500]; %stages' propellant mass (kg)

m_p = 732.8; %Mass of the payload at launch (kg)
m0 = sum(ms) + sum(mp) + m_p; %Total mass of the rocket at lift-off (kg)
%% Phases

% 1.
% Initial condition
stage = 1; %Stage currently use by the rocket
V = 0; % (m/s)
gamma = 0; % (rad)
x = 0; % (m)
h = 0; % (m) Sea level, to adpat in function of the sarting point
m = m0;
y0 = [V; gamma; x; h; m]; % State vector, 1 column

ti1 = 10; % Final time for phase 1 (s). Must be a short time (i for intermediate)
% On a paper we got 12s.
param = [Isp(stage), Cd, A, stage];
%We need to estimate the thrust delivered by the engines
[t, y] = ode45(@(t, y) ascent_dynamicsODE(T(stage), y, param), [t0, ti1], y0, options);

%2.
gamma = 0.05*pi/180; % (rad) Non-zero value to start gravity turn
y01 = y;
y01(2) = gamma;
tb1 = mp(stage)*g0*Isp(stage)/T(stage); %burnout time of 1st stage (s)
tf1 = tb1 - ti1; %end of the 1st stage phase
[t1, y1] = ode45(@(t, y) ascent_dynamicsODE(T(stage), y, param), [ti1, tf1], y01, options);

%3.
stage = 2;
y02 = y;
y02(5) = y02(5)-ms(stage-1); %1st stage removal
param(1) = Isp(stage); %Isp update for 2nd stage
param(5) = stage;
tb2 = mp(stage)*g0*Isp/T(stage); %burnout time of 2nd stage (s)
tf2 = tb2 - tf1; %end of the 2nd stage phase
[t2, y2] = ode45(@(t, y) ascent_dynamicsODE(T(stage), y, param), [tf1, tf2], y02, options);

%4.
stage = 3;
y03 = y;
y03(5) = y03(5)-ms(stage-1); %2nd stage removal
param(1) = Isp(stage); %Isp update for 2nd stage
param(5) = stage;
tb3 = mp(stage)*g0*Isp(stage)/T(stage); %burnout time of 2nd stage (s)
tf3 = tb3 - tf2; %end of the 2nd stage phase
[t3, y3] = ode45(@(t, y) ascent_dynamicsODE(T(stage), y, param), [tf2, tf3], y03, options);