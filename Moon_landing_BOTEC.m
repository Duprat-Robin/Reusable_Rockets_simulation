%% Gunnar Tibert, KTH, 2020-04-20

%% Moon landing BOTEC
clear, clc

mu_M=4.9048695e12; %% gravitational parameter of the Moon (m^3s^-2)
R_M=1737e3; %% mean radius of the Moon (m)
g_M=1.625; %% gravity of the Moon (m/s^2)
h1=100e3; %% altitude 1 (m)
h2=20e3; %% altitude 2 (m)

disp('=======================================================================');
disp('Hohmann transfer from 100 km to 20 km altitude');
r1=R_M+h1; %% radius at apolune
r2=R_M+h2; %% radius at perilune
at=(r1+r2)/2; %% semi major axis of transfer ellipse

disp('Specific impulse of the spacecraft propulsion system (s)');
Isp_sc=321;
Isp_sc

disp('Orbital velocity (km/s) of 100 km circular orbit');
V1=sqrt(mu_M/r1);
V1/1e3

disp('Orbital velocity (km/s) at 100 km of elliptical orbit with 100 km apolune and and 20 km perilune');
V2=sqrt(mu_M*(2/r1-1/at));
V2/1e3

disp('Hohmann transfer DV (m/s) to lower altitude from 100 km circular to 20 km perilune');
DV_Hohmann=V1-V2;
DV_Hohmann

disp('Required DV (m/s) with 100% margin');
DV_Hohmann*2

disp('Spacecraft mass (kg) at start of Hohmann transfer');
m_sc0=573

disp('Fuel mass (kg) for the Hohmann transfer with 100% DV margin');
m_fuel_Hohmann=m_sc0*(1-exp(-2*DV_Hohmann/Isp_sc/9.81))


disp('=======================================================================');
disp('Back-of-the-Envelope-Analysis for the powered descent from 20 km perilune');

disp('Assumed total time for the powered descent (min)');
t_des=13*60;
t_des/60

disp('Orbital velocity (km/s) at 20 km perilune');
DV_hor=sqrt(mu_M*(2/r2-1/at));
DV_hor/1e3

disp('Vertical speed gain (km/s) for a constant mass falling from 20 km in the Moons gravity field (constant g_Moon)');
DV_ver=sqrt(2*g_M*h2);
DV_ver/1e3

%disp('Total speed (km/s) to brake away is approximated as sqrt(DV_hor^2+DV_ver^2)');
%DV_tot=sqrt(DV_hor^2+DV_ver^2);
%DV_tot/1e3

disp('Total speed (km/s) to brake away is approximated as DV_hor+DV_ver)');
DV_tot=DV_hor+DV_ver;
DV_tot/1e3


disp('Required constant acceleration (m/s^2) to brake from orbital speed to zero speed');
a_des=DV_tot/t_des;
a_des

disp('Assumed spacecraft mass (kg) at the start of the powered descent');
m_sc=566;
m_sc

disp('Required fuel mass (kg) to brake away all speed with 5% DV margin');
m_fuel_des=m_sc*(1-exp(-DV_tot*1.05/Isp_sc/9.81));
m_fuel_des

disp('Constant mass flow (kg/s) during the powered descent');
dmdt_fuel_des=m_fuel_des/t_des;
dmdt_fuel_des

disp('Required constant thrust(N) for the average mass flow');
F_des=dmdt_fuel_des*Isp_sc*9.81;
F_des

disp('=======================================================================');
disp('Cost of final vertical descent/hovering at 100 m altitude');

disp('Hovering time (min)');
t_hov=3*60;
t_hov/60

disp('Mass of spacecraft (kg) at start of hovering');
m_hov=m_sc-m_fuel_des;
m_hov

disp('Maximum thrust (N) to hover at constant altitude');
F_hov=m_hov*g_M;
F_hov

disp('Maximum required mass flow (kg/s) to maintain hovering thrust');
dmdt_fuel_hov=F_hov/Isp_sc/9.81;
dmdt_fuel_hov

disp('Maximum required fuel mass (kg) to hover');
m_fuel_hov=dmdt_fuel_hov*t_hov

disp('Total DV (m/s) for the hovering');
DV_hov=-Isp_sc*9.81*log((m_hov-m_fuel_hov)/m_hov)

disp('Note that the required thrust is lower at the end of the hovering');
disp('phase as the mass is lower');
average_mass_ratio=(m_hov-m_fuel_hov+m_hov)/2/m_hov

disp('Required fuel mass (kg) to hover with average mass and 5% mass margin');
m_fuel_hov*average_mass_ratio*1.05











