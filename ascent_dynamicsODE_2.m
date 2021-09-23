function dy = ascent_dynamicsODE_2(t, T, y, param, gammas, tf)
%ASCENT_DYNAMICS: This function provide the differential equations for the
%dynamics of the rockets during the ascention phase.
%   The function return a column vector dy corresponding to [acceleration,
%   flight path angle, ground distance rate from lift-off, altitude
%   rate, mass flow rate].
%   Units: [dx = m/s, ddx = m/s^2, dh = m/s, ddh = m/s^2, dm = kg/s]. 
%   We assume that%   the rocket doesn't generate lift. 
%   Forces at stake: thrust T, drag D, weight m*g.
%   y is the state vector, dy is the time derivative of y

if ~exist('gammas', 'var')
    gammas = [0, 0];
end

if ~exist('tf', 'var')
    tf = [0, 10];
end
ix = 1; idx = 2; ih = 3; idh = 4; im = 5;
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
g = mu_E/((Re+y(ih))^2); %Earth model: gravitational acceleration in function of the alitude

[Temp, sound_vel, P, rho] = atmoscoesa(y(ih), 'None'); %Matlab atmospheric model

D = 0.5*A*rho*Cd*(y(idx)^2+y(idh)^2); % Drag (N)
q = 0.5*rho*(y(idx)^2+y(idh)^2); %dynamic pressure

if phase == 1
    gamma = gammas;
elseif stage ~= 3 && phase ~= 4
    gamma = acos(y(idx)/(sqrt(y(idx)^2+y(idh)^2))); %flight path angle
elseif phase == 4
    gamma = atan(tan(gammas)*(1-(t-tf(1))/tf(2))); %in this phase, gamma is no longer a variable because of steering law
end

dy(ix) = y(idx);
dy(idx) = (T+D)*cos(gamma)/y(im) - y(idx)*y(idh)/(Re+y(ih));
dy(ih) = y(idh);
dy(idh) = (T+D)*sin(gamma)/y(im) - (g-y(idx)^2/(Re+y(ih)));
dy(im) = -abs(T)/Isp/g0; %mass flow rate (kg/s)

end

