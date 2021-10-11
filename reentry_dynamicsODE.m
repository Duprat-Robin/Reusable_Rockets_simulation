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
global landing_printed flight_time

%% Rocket parameters
Isp = param(1); %Isp (s). TBD
Cd = param(2); %Drag coefficient. 1st assumption: the rocket is a cylinder (cf. Wikipedia Drag Coefficient)
A = param(3); %Surface of the rocket in contact with the airflow (m^2)
stage = param(4);
phase = param(5);

iV = 1; igamma = 2; ih = 3; ix = 4; im = 5;
%% Rocket dynamaic equation in the Troposphere
dy = zeros(5,1);
g = mu_E/((Re+y(ih))^2); %Earth model: gravitational acceleration in function of the alitude

[Temp, sound_vel, P, rho] = atmoscoesa(y(ih), 'None'); %Matlab atmospheric model

if isnan(rho)
    rho = 0;
end

D = 0.5*A*rho*Cd*y(iV)^2; % Drag (N)

h_hoovering=149; % thrusting before landing
coef_thurst=0.99999;
if stage==2
    h_hoovering=87;
    coef_thurst=0.99;
end
if y(ih)<h_hoovering
    T=-coef_thurst*y(im)*g; %high thrust with the remaining propellant! (as soon as it does not exceed the maximal one.
end

if y(ih)<=0 %landed!
    dy(iV) = 0;
    dy(ih) = 0; %altitude rate (m/s)
    dy(ix) = 0; %ground distance rate (m/s)
    dy(im) = 0; %mass flow rate (kg/s)
    dy(igamma) =0;
    if ~landing_printed
        landing_printed=true;
        flight_time= t;
        disp('=======================================================================');
        disp(['Rocket stage ' num2str(stage) ' landed safely, landing parameters:']);
        disp(['Landing velocity: ' num2str(y(iV)) ' m/s']);
        disp(['Flight path angle: ' num2str(y(igamma)*180/pi) '°']);
        disp(['height: ' num2str(y(ih)) ' m']);
        disp(['Landing ground distance from launching site: ' num2str(y(ix)) ' m']);
        disp(['Final mass: ' num2str(y(im)) ' kg']);
        disp(['Total flight time: ' num2str(flight_time) ' s']);

    end
    return
else  
    dy(iV) = (T-D)/y(im) - g*sin(y(igamma)); %acceleration (m/s^2)
    if phase == 7.1 || phase==7.3 %turning phase
        dy(igamma) = (gammas(2)-gammas(1))/(tf(2)-tf(1)); %Linear progression for dgamma during 1st phase
    elseif phase==7.4 || phase==7.5 %going down phase
        dy(igamma) = 0; %flight path angle fixed to be pi/2 (1/s)
    else
        dy(igamma) = -1/y(iV) * (g-(y(iV)^2)/(Re+y(ih)))*cos(y(igamma));
    end
    dy(ih) = y(iV)*sin(y(igamma)); %altitude rate (m/s)
    dy(ix) = Re*y(iV)*cos(y(igamma))/(Re+y(ih)); %ground distance rate (m/s)
    dy(im) = -abs(T)/Isp/g0; %mass flow rate (kg/s)
end



end