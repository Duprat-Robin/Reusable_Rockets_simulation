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

%Taking into account the Earth rotating speed V_E: The speed given in the
%state vector is related to the Earth rotating referential. So y(iV) is in
%regards to the launching point. Thuse, we have to substract V_E to the
%orbital speed. Then we need a smaller DeltaV to reach orbt. This DeltaV is
%provided by the Earth itself. That is why we have an interest with launch
%sites close to the equator.

if ~exist('gammas', 'var')
    gammas = [0, 0];
end

if ~exist('tf', 'var')
    tf = [0, 10];
end
iV = 1; igamma = 2; ih = 3; ix = 4; im = 5; %state vector index
%% Earth parameters
Re = 6371e3; %Earth radius (m). To adapt with the Launch point (average 6371km)
g0=9.80665; %gravitational acceleration on Earth at sea level (m/s^2)
mu_E = 3.986e14; %Earth's gravitationnal constant (m^3/s^2)
i0 = 5.2308333*pi/180; %Kourou inlcination (deg)
CirE = 2*pi*Re*cos(i0);
V_E = CirE/(24*3600); %m/s Earth rotation velocity at Kourou, ~461.8897 m/s

%% Rocket parameters
Isp = param(1); %Isp (s). TBD
Cd = param(2); %Drag coefficient. 1st assumption: the rocket is a cylinder (cf. Wikipedia Drag Coefficient)
A = param(3); %Surface of the rocket in contact with the airflow (m^2)
stage = param(4);
phase = param(5);

%% Final conditions
Hf = 200e3; %altitude of the parking orbit
Vf = sqrt(mu_E/(Re+Hf))-V_E; %Speed related to the launch site (moving with the Earth at a velocity V_E)
%% Rocket dynamaic equation in the Troposphere
dy = zeros(5,1);
g = mu_E/((Re+y(ih))^2); %Earth model: gravitational acceleration in function of the alitude

[Temp, sound_vel, P, rho] = atmoscoesa(y(ih), 'None'); %Matlab atmospheric model

if isnan(rho)
    %disp(y(ih)) %to check where this equation can't provide a value for rho
    rho = 0;
end

D = 0.5*A*rho*Cd*y(iV)^2; % Drag (N)

if phase == 4 || y(igamma) <= 0
    y(igamma) = atan(tan(gammas(1))*(1-(t-tf(1))/(tf(2)-tf(1)))); %Steering law: linear tangent law
end

if y(ih) >= Hf && y(iV) >= Vf
    T = 0; %On the parking with orbital speed, we stop the thrust
end

dy(iV) = (T-D)/y(im) - (g-(y(iV)^2)/(Re+y(ih)))*sin(y(igamma)); %acceleration (m/s^2)

if phase == 2 && stage == 1
    dy(igamma) = (-cos(y(igamma))/y(iV))*(g-(y(iV)^2)/(Re+y(ih))); %during phase 2, natural gravity turn
elseif ~(phase == 4 || y(igamma) <= 0)
    dy(igamma) = (gammas(2)-gammas(1))/(tf(2)-tf(1)); %Linear progression for dgamma during 2nd and 3rd phases
else
     dy(igamma) = 0; %in this phase, gamma is no longer a variable because of steering law
end

dy(ih) = y(iV)*sin(y(igamma)); %altitude rate (m/s)

%In an ideal world, we want gamma = 0 just before reaching Vf and Hf
if y(ih) >= Hf
     dy(ih) = 0;
end

if y(iV) >= Vf
     dy(iV) = 0;
end

dy(ix) = Re*y(iV)*cos(y(igamma))/(Re+y(ih)); %ground distance rate (m/s)
dy(im) = -abs(T)/Isp/g0; %mass flow rate (kg/s)

end

