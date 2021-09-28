function rho = atmosphere(h)
%ATMOSPHERE Summary of this function goes here
%   Detailed explanation goes here
%% Standard atmosphere
p0 = 101325; %Pressure at sea level (Pa)
T0 = 288.15; %Temperature at sea level (K)
R = 287.058; %Specific gas constant for dry air (J/(kg.K))
M = 2.89652e-2; %Molar mass of dry air (kg/mol)
L = 6.e-3; %(K/m) linear evolution of temperature with atm in the Troposphere

H_tropopause = 11e3; %altitude of the tropopause (m)

%% Models

if h<H_tropopause
    rho = (p0*M/(R*T0)) * (1-L*y(3)/T0)^(g*M/(R*L)-1); %air density at altitude h=y(3) (kg/m^3) in the Troposphere
end

