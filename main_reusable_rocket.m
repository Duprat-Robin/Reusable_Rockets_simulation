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
    % Final conditions (for the rocket and not the mission): gamma=0, H=H*, V=V*

%% Constant values

mu_E=3.968e14; % gravitational parameter of Earth (m^3s^-2)
Re=6378e3; % mean radius of Earth (m)
g0=9.80665; % gravity of Earth at sea level (m/s^2)

m0 = 1000e3; % initiale mass of the rocket (kg)
t0 = O; % ignition time of the 1st stage (s)

%% Rocket parameters
Isp = [400, 467, 400]; %Isp (s) for 1st stage. TBD
Cd = 0.85; %Drag coefficient. 1st assumption: the rocket is a cylinder (cf. Wikipedia Drag Coefficient)
A = 1; %Surface of the rocket in contact with the airflow (m^2)

T = [11e6, 11e6]; %stages trhrust (N)

%% Phases

% 1.
% Initial condition
V = 0; % (m/s)
gamma = 0; % (rad)
x = 0; % (m)
h = 0; % (m) Sea level, to adpat in function of the sarting point
m = m0;
y0 = [V; gamma; x; h; m]; % State vector, 1 column

tf1 = 10; % Final time for phase 1 (s). Must be a short time
% On a paper we got 12s.
param = [Isp(1), Cd, A];
%We need to estimate the thrust delivered by the engines
[t, y] = ode45(@(t, y) ascent_dynamicsODE(T(1), y0, param), [to, tf1], y0);

%2.
gamma = 0.05*pi/180; % (rad) Non-zero value to start gravity turn
y1 = y;
y1(2) = gamma;
tb1 = mp*g0*Isp/T1; %burnout time of 1st stage (s)

[t, y] = ode45(@(t, y) ascent_dynamicsODE(T(2), y1, param), [tf1, tb1-tf1], y1);
