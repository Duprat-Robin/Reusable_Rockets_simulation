%% Gunnar Tibert, KTH, 2020-04-20

%% Simulating the powered descent landing on the Moon
clear,clc,clf

%% Initial speed and altitude
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

disp('Orbital velocity (km/s) at 20 km of elliptical orbit with 100 km apolune and and 20 km perilune');
V0=sqrt(mu_M*(2/r2-1/at));

disp('Initial mass at start of powered descent');
m0=566

options = odeset('RelTol',1e-10,'AbsTol',1e-9);


%% time integration of the braking phase
tbraking=[0 10*60]; %% braking during the first 10 min
T1=-1200; %% braking thrust (N)
gamma1=[0 -4.8]*pi/180; %% start and end flight path angles

tf1=tbraking; %% start and end times
y0braking=[V0 0 h2 0 m0]; %% starting values for speed-flight path angle-altitude-ground distance-mass
[t1,y1]=ode45(@(t,y) MoonDescentODE(t,y,T1,gamma1,tf1),tbraking,y0braking,options);

%% time integration of the controlled descent phase
tdescent=[10*60 13*60]; %% controlled descent during 3 min
m0_des=y1(end,5) %% starting mass for the controlled descent phase
T2_est=-m0_des*g_M %% estimated thrust to brake the speed = m*g
scalefactor=0.955 %% scaling down the constant thrust
T2=scalefactor*T2_est; %% braking thrust (N)
gamma2=[y1(end,2)*180/pi -89.9]*pi/180; %% start and end flight path angles
tf2=tdescent; %% start and end times

y0descent=[y1(end,1) y1(end,2) y1(end,3) y1(end,4) y1(end,5)]; %% starting values for speed-flight path angle-altitude-ground distance-mass
[t2,y2]=ode45(@(t,y) MoonDescentODE(t,y,T2,gamma2,tf2),tdescent,y0descent,options);


%% Computing the Delta-V from the fuel consumption
Isp=321; %% specific impulse of engine (s)
g0=9.81;
m0
y1(end,5)
y2(end,5)
DVbrake=-Isp*g0*log(y1(end,5)/m0)
DVdescent=-Isp*g0*log(y2(end,5)/y1(end,5))

%% plotting the results for the total powered descent

figure(1); hold on;
plot(t1,y1(:,1),'r','LineWidth',2);
plot(t2,y2(:,1),'g','LineWidth',2);
title('Speed change during the powered descent');
xlabel('Time (s)');
ylabel('Speed (m/s)');
grid;
legend('Braking phase','Descent phase');

figure(2); hold on;
plot(t1,y1(:,2)*180/pi,'r','LineWidth',2);
plot(t2,y2(:,2)*180/pi,'g','LineWidth',2);
title('Flight path angle change during the powered descent');
xlabel('Time (s)');
ylabel('Flight path angle (deg)');
grid;
legend('Braking phase','Descent phase');

figure(3); hold on;
plot(t1,y1(:,3)/1e3,'r','LineWidth',2);
plot(t2,y2(:,3)/1e3,'g','LineWidth',2);
title('Altitude change during the powered descent');
xlabel('Time (s)');
ylabel('Altitude (km)');
grid;
legend('Braking phase','Descent phase');

figure(4); hold on;
plot(t1,y1(:,5),'r','LineWidth',2);
plot(t2,y2(:,5),'g','LineWidth',2);
title('Spacecraft mass change during the powered descent');
xlabel('Time (s)');
ylabel('Mass (kg)');
grid;
legend('Braking phase','Descent phase');

figure(5); hold on;
plot(y1(:,4)/1e3,y1(:,3)/1e3,'r','LineWidth',2);
plot(y2(:,4)/1e3,y2(:,3)/1e3,'g','LineWidth',2);
title('Altitude vs ground distance during the powered descent');
xlabel('Ground distance (km)');
ylabel('Altitude (km)');
grid;
legend('Braking phase','Descent phase');

figure(6); hold on;
plot(t1,abs(T1)*ones(size(t1)),'r','LineWidth',2);
plot(t2,abs(T2)*ones(size(t2)),'g','LineWidth',2);
title('Braking thrust for the powered descent');
ylabel('Braking thrust (N)');
xlabel('Time (s)');
grid;
legend('Braking phase','Descent phase');






