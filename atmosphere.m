function rho = atmosphere(h)
%ATMOSPHERE Summary of this function goes here
%   Detailed explanation goes here
%% Standard atmosphere
rho0 = 1.225; %air density at sea level (kg/m^3)

%% Models
H0 = 8.5e3; %atmosphere scale height

rho = rho0*exp(-h/H0);
end

