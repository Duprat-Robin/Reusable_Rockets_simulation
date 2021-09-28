function dy = reentry_dynamicsODE(t, T, y, param, gammas, tf)
%ASCENT_DYNAMICS: This function provide the differential equations for the
%dynamics of the stages during the reentry phase.
%   The function return a column vector dy corresponding to [acceleration,
%   flight path angle, ground distance rate from lift-off, altitude
%   rate, mass flow rate].
%   Units: [a = m/s^2, dgamma = 1/s, dh = m/s, dx = m/s, dm = kg/s]. 
%   the stages might generate lift due to the parachutes, or it will be modelised by an added thrust/higher drag. 
%   Forces at stake: thrust T, drag D, weight m*g.
%   y is the state vector, dy is the time derivative of y


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
g = mu_E/((Re+y(3))^2); %Earth model: gravitational acceleration in function of the alitude

[Temp, sound_vel, P, rho] = atmoscoesa(y(3), 'None'); %Matlab atmospheric model

D = 0.5*A*rho*Cd*y(1)^2; % Drag (N)

dy(1) = (T-D)/y(5) - g*sin(y(2)); %acceleration (m/s^2)
if phase == 7.1 || phase==7.3 %turning phase
    dy(2) = (gammas(2)-gammas(1))/(tf(2)-tf(1)); %Linear progression for dgamma during 1st phase
else %going down phase
    dy(2) = 0; %flight path angle fixed to be pi/2 (1/s)
end
%dy(2)=(gamma(2)-gamma(1))/(tf(2)-tf(1)); %% flight path angle (1/s)
dy(3) = y(1)*sin(y(2)); %altitude rate (m/s)
dy(4) = Re*y(1)*cos(y(2))/(Re+y(3)); %ground distance rate (m/s)
dy(5) = -abs(T)/Isp/g0; %mass flow rate (kg/s)



end