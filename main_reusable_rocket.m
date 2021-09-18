%% Description
% Gravity turn (cf. SD2900 - Rocket Dynamics, slide 25)
    %5 phases:
    % 1. Let the 1st staghe go vertical for a short time to pick up speed to
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
    
%% Constant values

mu_E=3.968e14; % gravitational parameter of Earth (m^3s^-2)
R_E=6378e3; % mean radius of Earth (m)
g_E=9.80665; % gravity of Earth at sea level (m/s^2)

m0 = 1000e3; % initiale mass of the rocket (kg)
t0 = O; % ignition time of the 1st stage (s)

%% Phases

% 1.
% Initial condition*
V = 0; % (m/s)
gamma = 0; % (rad)
x = 0; % (m)
h = 0; % (m) Sea level, to adpat in function of the sarting point
m = m0;
y = [V; gamma; x; H; m]; % State vector, 1 column
