function dy = ascent_dynamicsODE(t, T, y, param, gammas, tf)
%ASCENT_DYNAMICS: This function provide the differential equations for the
%dynamics of the rockets during the ascention phase.
%   The function return a column vector dy corresponding to [acceleration,
%   flight path angle, ground distance rate from lift-off, altitude
%   rate, mass flow rate].
%   Units: [a = m/s^2, dgamma = 1/s, dh = m/s, dx = m/s, dm = kg/s]. 
%   We assume that%   the rocket doesn't generate lift. 
%   Forces at stake: thrust T, drag D, weight m*g.
%   y is the state vector, dy is the time derivative of y

if ~exist('gammas', 'var')
    gammas = [0, 0];
end

if ~exist('tf', 'var')
    tf = [0, 10];
end

%% Earth parameters
Re = 6371e3; %Earth radius (m). To adapt with the Launch point (average 6371km)
g0=9.80665; %gravitational acceleration on Earth at sea level (m/s^2)
mu_E = 3.986e14; %Earth's gravitationnal constant (m^3/s^2)

%% Rocket parameters
Isp = param(1); %Isp (s). TBD
Cd = param(2); %Drag coefficient. 1st assumption: the rocket is a cylinder (cf. Wikipedia Drag Coefficient)
A = param(3); %Surface of the rocket in contact with the airflow (m^2)
stage = param(4);
phase = param(5);
%% Rocket dynamaic equation in the Troposphere
dy = zeros(5,1);
g = mu_E/((Re+y(3))^2); %Earth model: gravitational accelration in function of the alitude

[Temp, sound_vel, P, rho] = atmoscoesa(y(2)); %Matlab atmospheric model

D = 0.5*A*rho*Cd*y(1)^2; % Drag (N)

dy(1) = (T-D)/y(5) - g*sin(y(2)); %acceleration (m/s^2)
if phase == 1
    dy(2) = (gammas(2)-gammas(1))/(tf(2)-tf(1)); %Linear progression for dgamma during 1st phase
else
    dy(2) = -1/y(1) * (g-(y(1)^2)/(Re+y(4)))*cos(y(2)); %flight path angle (1/s)
end
% if stage ~= 3
%     dy(2) = -1/y(1) * (g-(y(1)^2)/(Re+y(4)))*cos(y(2)); %flight path angle (1/s)
% else
%     dy(2) = 0; %in this phase, gamma is no longer a variable because of steering law
% end
dy(3) = Re*y(1)*cos(y(2))/(Re+y(4)); %ground distance rate (m/s)
dy(4) = y(1)*sin(y(2)); %altitude rate (m/s)
dy(5) = -abs(T)/Isp/g0; %mass flow rate (kg/s)

end

