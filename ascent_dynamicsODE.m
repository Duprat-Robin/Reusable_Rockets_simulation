function [dy] = ascent_dynamicsODE(T, y)
%ASCENT_DYNAMICS: This function provide the differential equations for the
%dynamics of the rockets during the ascention phase.
%   The function return a column vector dy corresponding to [acceleration,
%   flight path angle, ground distance rate from lift-off, altitude
%   rate, mass flow rate].
%   Units: [a = m/s^2, dgamma = 1/s, dh = m/s, dx = m/s, dm = kg/s]. 
%   We assume that%   the rocket doesn't generate lift. 
%   Forces at stake: thrust T, drag D, weight m*g.
%   y is the state vector, dy is the time derivative of y

%% Earth parameters
Re = 6371e3; %Earth radius (m). To adapt with the Launch point (average 6371km)
g0=9.80665; %gravitational acceleration on Earth at sea level (m/s^2)

%% Standard atmosphere
p0 = 101325; %Pressure at sea level (Pa)
T0 = 288.15; %Temperature at sea level (K)
R = 287.058; %Specific gas constant for dry air (J/(kg.K))
M = 2.89652e-2; %Molar mass of dry air (kg/mol)
L = 6.e-3; %(K/m) linear evolution of temperature with atm in the Troposphere
%% Rocket parameters
Isp = 400; %Isp (s). TBD
Cd = 0.85; %Drag coefficient. 1st assumption: the rocket is a cylinder (cf. Wikipedia Drag Coefficient)
A = 1; %Surface of the rocket in contact with the airflow (m^2)
%% Rocket dynamaic equation in the Troposphere
dy = zeros(5,1);
g = g0*(y(4)/(Re+y(4)))^2; %Earth model: gravitational accelration in function of the alitude

rho = (p0*M/(R*T0)) * (1-L*y(3)/T0)^(g*M/(R*L)-1); %air density at altitude h=y(3) (kg/m^3) in the Troposphere
D = 0.5*rho*A*Cd*y(1)^2; % Drag (N)

dy(1) = (T-D)/y(5) - g*sin(y(2)); %acceleration (m/s^2)
dy(2) = -1/y(1) * (g-(y(1)^2)/(Re+y(4)))*cos(y(2)); %flight path angle (1/s)
dy(3) = Re*y(1)*cos(y(2))/(Re+y(4)); %ground distance rate (m/s)
dy(4) = y(1)*sin(y(2)); %altitude rate (m/s)
dy(5) = -abs(T)/Isp/g0; %mass flow rate (kg/s)

end

