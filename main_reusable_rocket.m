%% Description
% Gravity turn (cf. SD2900 - Rocket Dynamics, slide 25)
    %Assuming Lift = 0
    %5 phases:
    % 1. Let the 1st stage go vertical for a short time to pick up speed to
    % avoid a too large dgamma/dt
    % 2. Set the flight path angle to a small non-zero value to start the
    % gravity turn
    % 3. At the burnout of the 1st stage, stop the time integration and restart
    % the integration with data for 2nd stage
    % 4. At the burnout of 2nd stage, stop the time integration and restart the
    % integration with data for 3rd stage.
    % 5. Use the tangent steering law for the 3rd stage to enforce gamma=0 at
    % the correct altitude. Note that gamma is NOT a variable for 3rd stage is
    % used.
    % Final conditions (for the rocket and not the mission): gamma=0, H=H*, V=V*
    % 6. Once we are at this parking orbit, we begin the Hohmann transfer to reach the target orbit.   
    % 7. We compute the reentry phases using the previous elements.  
    
clear,clc,clf, close all

%% Constant values

mu_E=3.968e14; % gravitational parameter of Earth (m^3s^-2)
Re=6378e3; % mean radius of Earth (m)
g0=9.80665; % gravity of Earth at sea level (m/s^2)
Hf = 200e3; %parking orbit altitude (m)
Xf = 200e3; %distance from launch pad at the end of ascent phase (m). Usefull for BVP

t0 = 0; % ignition time of the 1st stage (s)

options = odeset('RelTol',1e-10,'AbsTol',1e-9);
iV = 1; igamma = 2; ih = 3; ix = 4; im = 5;
%% Rocket parameters
Isp = [378, 359, 467]; %Isp (s) for 1st stage. TBD
Cd = 0.85; %Drag coefficient. 1st assumption: the rocket is a cylinder (cf. Wikipedia Drag Coefficient)
A = pi*3.66^2; %Surface of the rocket in contact with the airflow (m^2)
nb_engines = [3, 1, 1];

T = [2205000, 533000, 180000].*nb_engines; %stages' thrust (N)
ms = [16e3, 4e3, 1.5e3]; %stages' strucutal mass (kg)
mp = [418647.6, 31134, 3751.57]; %stages' propellant mass (kg)
m_fuel_reentry = [3750, 45, 0]; %stages' propellant mass need for reentry, includeded in mp. TBD(kg)

m_p = 732.8; %Mass of the payload at launch (kg)
m0 = sum(ms) + sum(mp) + m_p; %Total mass of the rocket at lift-off (kg)

V0 = 0; % (m/s)
gamma0 = pi/2; % (rad)
gamma1 = gamma0 - 0.1*pi/180;
x0 = 0; % (m)
h0 = 0; % (m) Sea level, to adpat in function of the sarting point
%% Phases

% 1.
% Initial condition
phase = 1; %Current phase
stage = 1; %Stage currently use by the rocket
param = [Isp(stage), Cd, A, stage, phase];

y01 = [V0 gamma0 h0 x0 m0]; % Initial state vector, 1 line

ti1 = 10; % Final time for phase 1 (s). Must be a short time (i for intermediate)
% On a paper we got 12s.
[t1, y1] = ode45(@(t, y) ascent_dynamicsODE(t, T(stage), y, param, [gamma0, gamma1], [t0, ti1]), [t0 ti1], y01, options);

%2.
phase = 2;
param = [Isp(stage), Cd, A, stage, phase];

y02 = [y1(end,iV) y1(end,igamma) y1(end,ih) y1(end,ix) y1(end,im)];

tb1 = (mp(stage)-m_fuel_reentry(stage))*g0*Isp(stage)/T(stage); %burnout time of 1st stage (s)
tf1 = tb1; %end of the 1st stage phase
[t2, y2] = ode45(@(t, y) ascent_dynamicsODE(t, T(stage), y, param), [ti1 tf1], y02, options);

%3.
phase = 3;
stage = 2;
param = [Isp(stage), Cd, A, stage, phase];

y03 = [y2(end,iV) y2(end,igamma) y2(end,ih) y2(end,ix) y2(end,im)];
y03(end,im) = y03(end,im)-ms(stage-1)-m_fuel_reentry(stage-1); %1st stage removal

tb2 = (mp(stage)-m_fuel_reentry(stage))*g0*Isp(stage)/T(stage); %burnout time of 2nd stage (s)
tf2 = tb2 + tf1; %end of the 2nd stage phase
[t3, y3] = ode45(@(t, y) ascent_dynamicsODE(t, T(stage), y, param), [tf1 tf2], y03, options);

%4.
phase = 4;
stage = 3;
param = [Isp(stage), Cd, A, stage, phase];

y04 = [y3(end,iV) y3(end,igamma) y3(end,ih) y3(end,ix) y3(end,im)];
y04(end,im) = y04(end,im)-ms(stage-1)-m_fuel_reentry(stage-1); %2nd stage removal

tb3 = (mp(stage)-m_fuel_reentry(stage))*g0*Isp(stage)/T(stage); %burnout time of 2nd stage (s)
tf3 = tb3 + tf2; %end of the 3rd stage phase
[t4, y4] = ode45(@(t, y) ascent_dynamicsODE(t, T(stage), y, param), [tf2 tf3], y04, options);
for i=1:size(y4(:,2),1)
    y4(i,igamma) = atan(tan(y3(end,igamma)*(1-(t4(i)-tf2)/tb3))); %Steering law: linear tangent law
end

% 6. Hohmann transfer : using 3rd stage
%Re=6378e3; % mean radius of Earth (m)
h= Re + y4(end,ih); % Assuming a parking orbit of 200km.
%mu = 3.986004418e14; % constant, assuming M+m is approximately M and constant
h_target = 29599.8e3; % target orbit
m_init=y4(end,im); %initial mass for phase 6. 2100kg for the moment, will be y(5) in the end

%m_star=732.8;

%g0=9.80665; 

Vc1 = sqrt(mu_E/h);% current speed on parking orbit
Vc2 = sqrt(mu_E/h_target); % to be speed on MEO
a = (h+h_target)/2; %semi major axis of the transfer orbit
V_perigee = sqrt(mu_E*(2/h-1/a)); %velocity at the perigee of the transfer orbit 
V_apogee = sqrt(mu_E*(2/h_target-1/a)); %velocity at the apogee of the transfer orbit
Delta_V_perigee = V_perigee - Vc1;
Delta_V_apogee = Vc2 - V_apogee;
Delta_V = Delta_V_perigee + Delta_V_apogee; % cost of the total transfer 
Delta_t = pi*sqrt(a^3/mu_E); % transfer time

Delta_m = m_init-m_init*exp(-Delta_V/(Isp(stage)*g0));
disp(Delta_m);
m_s=Delta_m-m_p; %in general m_s = 1/7 m_p

%% Ploting phase for payload
plot_payload = false;

if plot_payload
    figure(1); hold on;
    plot(t1,y1(:,ih)/1e3,'r','LineWidth',2);
    plot(t2,y2(:,ih)/1e3,'g','LineWidth',2);
    plot(t3,y3(:,ih)/1e3,'b','LineWidth',2);
    plot(t4,y4(:,ih)/1e3,'y','LineWidth',2);
    title('Altitude change for payload');
    xlabel('Time (s)');
    ylabel('Altitude (km)');
    grid;

    figure(2); hold on;
    plot(t1,y1(:,igamma)*180/pi,'r','LineWidth',2);
    plot(t2,y2(:,igamma)*180/pi,'g','LineWidth',2);
    plot(t3,y3(:,igamma)*180/pi,'b','LineWidth',2);
    plot(t4,y4(:,igamma)*180/pi,'y','LineWidth',2);
    title('Flight path angle change for payload');
    xlabel('Time (s)');
    ylabel('Flight path angle (deg)');
    grid;

    figure(3); hold on;
    plot(y1(:,ix)/1e3,y1(:,ih)/1e3,'r','LineWidth',2);
    plot(y2(:,ix)/1e3,y2(:,ih)/1e3,'g','LineWidth',2);
    plot(y3(:,ix)/1e3,y3(:,ih)/1e3,'b','LineWidth',2);
    plot(y4(:,ix)/1e3,y4(:,ih)/1e3,'y','LineWidth',2);
    title('Altitude change for payload');
    xlabel('X position (km)');
    ylabel('Altitude (km)');
    grid;

    figure(4); hold on;
    plot(t1,y1(:,iV)/1e3,'r','LineWidth',2);
    plot(t2,y2(:,iV)/1e3,'g','LineWidth',2);
    plot(t3,y3(:,iV)/1e3,'b','LineWidth',2);
    plot(t4,y4(:,iV)/1e3,'y','LineWidth',2);
    title('Speed change for payload');
    xlabel('Time (s)');
    ylabel('Speed (km/s)');
    grid;
end

 %% Simulating the powered descent landing of the reusable rocket

    
%% Constant values

% mu_E=3.968e14; % gravitational parameter of Earth (m^3s^-2)
% Re=6378e3; % mean radius of Earth (m)
% g0=9.80665; % gravity of Earth at sea level (m/s^2)
% 
% m0 = 1000e3; % initiale mass of the rocket (kg)
% t0 = 0; % ignition time of the 1st stage (s)

% 7. reentry of a stage 
% reentry heat and velocity + decrease speed : rocket + a
% parachute (take into account the m_s added, compared to the burned
% mp) to help braking
% approx 100m ground : hoovering with high thrusts to reduce speed to 0
% when touching V gamma x h m
global landing_printed
landing_printed = false;
reentry_stage=1;

Initial_speed= y03(iV) ;% for example 8000km/h ie, will be y(1) as soon as possible
Initial_gamma= y03(igamma); %Assuming a first stage separation at around pi/4, will be y(2) as soon as possible.
Initial_x = y03(ix); % Assuming a first stage separation at around 10km, will be y(4) as soon as possible.
Initial_h= y03(ih); % Assuming a first stage separation at around 100km, will be y(3) as soon as possible.
Initial_m = m_fuel_reentry(reentry_stage)+ms(reentry_stage); % 20000kg structural mass + 5000kg remaining propellant mass, will be y(5)- M_other_stages - M_payload
A=pi*3^2; %total rocket : pi*3.66^2
Cd=0.85;

%7.0 waiting phase
phase=7.0;
param = [Isp(reentry_stage), Cd, A, reentry_stage, phase];

y0_reentry_0 = [Initial_speed Initial_gamma Initial_h Initial_x Initial_m]; % Initial state vector, 1 line

ti_reentry_0 = tf1; % Initial time for phase 7.1 (s) ie end of burnout of phase 1
tf_reentry_0 = 13+ti_reentry_0; % Final time for phase 7.1 (s). TBD

[t_reentry_0, y_reentry_0] = ode45(@(t, y) reentry_dynamicsODE(t, 0, y, param), (ti_reentry_0:0.05:tf_reentry_0), y0_reentry_0, options);

%7.1 - turning phase : thursting to turn to opposite gamma
end_gamma=pi+74.09*pi/180;%TBD, 5° entry for the moment
phase=7.1;

param = [Isp(reentry_stage), Cd, A, reentry_stage, phase];

y0_reentry_1 = [y_reentry_0(end, iV) y_reentry_0(end,igamma) y_reentry_0(end,ih) y_reentry_0(end,ix) y_reentry_0(end,im)]; % Initial state vector, 1 line

ti_reentry_1 = tf_reentry_0; % Initial time for phase 7.1 (s) ie end of burnout of phase 1
tf_reentry_1 = 20+ti_reentry_1; % Final time for phase 7.1 (s). TBD

[t_reentry_1, y_reentry_1] = ode45(@(t, y) reentry_dynamicsODE(t, -T(reentry_stage)/100, y, param, [Initial_gamma, end_gamma], [ti_reentry_1, tf_reentry_1]), (ti_reentry_1:0.05:tf_reentry_1), y0_reentry_1, options);

%7.2 going in the opposite direction
phase=7.2;
param = [Isp(reentry_stage), Cd, A, reentry_stage, phase];

y0_reentry_2 = [y_reentry_1(end, iV) y_reentry_1(end,igamma) y_reentry_1(end,ih) y_reentry_1(end,ix) y_reentry_1(end,im)]; % Initial state vector, 1 line

ti_reentry_2 = tf_reentry_1; % Initial time for phase 7.2, end of phase 7.1 (s)
tf_reentry_2 = 57.2+ti_reentry_2; % Final time for phase 7.2 (s). TBD, for the moment 1 minute

[t_reentry_2, y_reentry_2 ] = ode45(@(t, y) reentry_dynamicsODE(t, -T(reentry_stage)/55 , y, param), (ti_reentry_2:0.05:tf_reentry_2), y0_reentry_2, options);

%7.3 turning to gamma =pi/2
phase=7.3;
end_gamma=3*pi/2;
param = [Isp(reentry_stage), Cd, A, reentry_stage, phase];

y0_reentry_3 = [y_reentry_2(end, iV) y_reentry_2(end,igamma) y_reentry_2(end,ih) y_reentry_2(end,ix) y_reentry_2(end,im)]; % Initial state vector, 1 line

ti_reentry_3 = tf_reentry_2; % Initial time for phase 7.3, end of phase 7.2 (s)
tf_reentry_3 = 15+ti_reentry_3; % Final time for phase 7.3 (s). TBD, for the moment 10 s

[t_reentry_3, y_reentry_3] = ode45(@(t, y) reentry_dynamicsODE(t, -T(reentry_stage)/80, y, param, [y0_reentry_3(igamma), end_gamma], [ti_reentry_3, tf_reentry_3]), (ti_reentry_3:0.05:tf_reentry_3), y0_reentry_3, options);

%7.4 - braking phase: parachute to brake
phase=7.4;
Cd_parachute = 1.5; %Drag coefficient for a space parachute (kevlar/nylon stuff, deployed at approx 10km to avoid it to brake)
A_parachute = 250; %Surface of the 3 parachutes in contact with the airflow (m^2) TBD
param = [Isp(reentry_stage), Cd_parachute, A_parachute, reentry_stage, phase];

y0_reentry_4 = [y_reentry_3(end, iV) y_reentry_3(end,igamma) y_reentry_3(end,ih) y_reentry_3(end,ix) y_reentry_3(end,im)]; % Initial state vector, 1 line

ti_reentry_4 = tf_reentry_3; % Initial time for phase 7.4, end of phase 7.3 (s)
tf_reentry_4 = 60*4+ti_reentry_4; % Final time for phase 7.4 (s). TBD, for the moment 3 minutes


[t_reentry_4, y_reentry_4] = ode45(@(t, y) reentry_dynamicsODE(t, 0 , y, param), (ti_reentry_4:0.05:tf_reentry_4), y0_reentry_4, options);

%7.5 - controlled descent phase : landing (high thrust). Maybe remove
%parachute ?
%phase=7.5;
%param = [Isp(reentry_stage), Cd, A, reentry_stage, phase];

%y0_reentry_5 = [y_reentry_4(end, iV) y_reentry_4(end,igamma) y_reentry_4(end,ih) y_reentry_4(end,ix) y_reentry_4(end,im)]; % Initial state vector, 1 line

%ti_reentry_5 = tf_reentry_4; % Initial time for phase 7.5, end of phase 7.4(s)
%tf_reentry_5 = 2+ti_reentry_5; % Final time for phase 7.5 (s). TBD, for the moment 20s


%[t_reentry_5, y_reentry_5] = ode45(@(t, y) reentry_dynamicsODE(t, -T(reentry_stage)/10, y, param), (ti_reentry_5:0.05:tf_reentry_5), y0_reentry_5, options);



%% Ploting phase
plot_stage1 = true;
add_ascent_plot=true;
lengend_list=["Waiting phase before reentry", "Turning phase to revert speed", "Reentry inside atmosphere phase", ...
    "Turning phase to land at 90°", "Deploy parachutes and thrusting before landing"];

if plot_stage1
    figure(5); hold on;
    if add_ascent_plot
        plot(t1,y1(:,ih)/1e3,'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',2);
        plot(t2,y2(:,ih)/1e3,'MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',2);
        lengend_list=["Ascent phase before gravity turn", "Thrusting phase during gravity turn", "Waiting phase before reentry",...
            "Turning phase to revert speed", "Reentry inside atmosphere phase", ...
            "Turning phase to land at 90°", "Deploy parachutes and thrusting before landing"] ;
    end
    plot(t_reentry_0,y_reentry_0(:,ih)/1e3,'m','LineWidth',2);
    plot(t_reentry_1,y_reentry_1(:,ih)/1e3,'r','LineWidth',2);
    plot(t_reentry_2,y_reentry_2(:,ih)/1e3,'MarkerFaceColor',[0.9290 0.6940 0.1250],'LineWidth',2);
    plot(t_reentry_3,y_reentry_3(:,ih)/1e3,'b','LineWidth',2);
    plot(t_reentry_4,y_reentry_4(:,ih)/1e3,'c','LineWidth',2);
    %plot(t_reentry_5,y_reentry_5(:,ih)/1e3,'m','LineWidth',2);
    title(['Altitude change for reentry of stage ' num2str(reentry_stage)]);
    xlabel('Time (s)');
    ylabel('Altitude (km)');
    set(findall(gcf,'-property','FontSize'),'FontSize',12)
    legend(lengend_list)
    grid;

    % figure(6); hold on;
    % plot(t_reentry_1,y_reentry_1(:,igamma)*180/pi,'r','LineWidth',2);
    % plot(t_reentry_2,y_reentry_2(:,igamma)*180/pi,'g','LineWidth',2);
    % plot(t_reentry_3,y_reentry_3(:,igamma)*180/pi,'b','LineWidth',2);
    % plot(t_reentry_4,y_reentry_4(:,igamma)*180/pi,'c','LineWidth',2);
    % %plot(t_reentry_5,y_reentry_5(:,igamma)*180/pi,'m','LineWidth',2);
    % title(['Flight path angle change for reentry of stage ' num2str(reentry_stage)]);
    % xlabel('Time (s)');
    % ylabel('Flight path angle (deg)');
    % grid;

    figure(7); hold on;
    if add_ascent_plot
        plot(y1(:,ix)/1e3,y1(:,ih)/1e3,'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',2);
        plot(y2(:,ix)/1e3,y2(:,ih)/1e3,'MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',2);
    end
    plot(y_reentry_0(:,ix)/1e3,y_reentry_0(:,ih)/1e3,'m','LineWidth',2);
    plot(y_reentry_1(:,ix)/1e3,y_reentry_1(:,ih)/1e3,'r','LineWidth',2);
    plot(y_reentry_2(:,ix)/1e3,y_reentry_2(:,ih)/1e3,'MarkerFaceColor',[0.9290 0.6940 0.1250],'LineWidth',2);
    plot(y_reentry_3(:,ix)/1e3,y_reentry_3(:,ih)/1e3,'b','LineWidth',2);
    plot(y_reentry_4(:,ix)/1e3,y_reentry_4(:,ih)/1e3,'c','LineWidth',2);
    %plot(y_reentry_5(:,ix)/1e3,y_reentry_5(:,ih)/1e3,'m','LineWidth',2);
    title(['Altitude change for reentry of stage ' num2str(reentry_stage)]);
    xlabel('X position (km)');
    ylabel('Altitude (km)');
    grid;

    figure(8); hold on;
    if add_ascent_plot
        plot(t1,y1(:,im),'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',2);
        plot(t2,y2(:,im),'MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',2);
    end
    plot(t_reentry_0,y_reentry_0(:,im),'m','LineWidth',2);
    plot(t_reentry_1,y_reentry_1(:,im),'r','LineWidth',2);
    plot(t_reentry_2,y_reentry_2(:,im),'MarkerFaceColor',[0.9290 0.6940 0.1250],'LineWidth',2);
    plot(t_reentry_3,y_reentry_3(:,im),'b','LineWidth',2);
    plot(t_reentry_4,y_reentry_4(:,im),'c','LineWidth',2);
    %plot(t_reentry_5,y_reentry_5(:,im),'m','LineWidth',2);
    yline(ms(reentry_stage),'-','Reentry stage structural mass');
    ylim([0.99*ms(reentry_stage),1.01*y_reentry_0(1,im)])
    title(['mass change for reentry of stage ' num2str(reentry_stage)]);
    xlabel('Time (s)'); 
    ylabel('Mass (kg)');
    grid;

    figure(9); hold on;
    if add_ascent_plot
        plot(t1,y1(:,iV),'MarkerFaceColor','#0072BD','LineWidth',2);
        plot(t2,y2(:,iV),'MarkerFaceColor','#D95319','LineWidth',2);
    end
    plot(t_reentry_0,y_reentry_0(:,iV),'m','LineWidth',2);
    plot(t_reentry_1,y_reentry_1(:,iV),'r','LineWidth',2);
    plot(t_reentry_2,y_reentry_2(:,iV),'MarkerFaceColor','#77AC30','LineWidth',2);
    plot(t_reentry_3,y_reentry_3(:,iV),'b','LineWidth',2);
    plot(t_reentry_4,y_reentry_4(:,iV),'c','LineWidth',2);
    %plot(t_reentry_5,y_reentry_5(:,iV),'m','LineWidth',2);
    title(['velocity change for reentry of stage '  num2str(reentry_stage)]);
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    grid;
end

 %% Simulating the powered descent landing of the reusable rocket

    
%% Constant values

% mu_E=3.968e14; % gravitational parameter of Earth (m^3s^-2)
% Re=6378e3; % mean radius of Earth (m)
% g0=9.80665; % gravity of Earth at sea level (m/s^2)
% 
% m0 = 1000e3; % initiale mass of the rocket (kg)
% t0 = 0; % ignition time of the 1st stage (s)

% 7. reentry of a stage 
% reentry heat and velocity + decrease speed : rocket + a
% parachute (take into account the m_s added, compared to the burned
% mp) to help braking
% approx 100m ground : hoovering with high thrusts to reduce speed to 0
% when touching V gamma x h m
landing_printed = false;
reentry_stage=2;

Initial_speed= y04(iV) ;% for example 8000km/h ie, will be y(1) as soon as possible
Initial_gamma= y04(igamma); %Assuming a first stage separation at around pi/4, will be y(2) as soon as possible.
Initial_x = y04(ix); % Assuming a first stage separation at around 10km, will be y(4) as soon as possible.
Initial_h= y04(ih); % Assuming a first stage separation at around 100km, will be y(3) as soon as possible.
Initial_m = m_fuel_reentry(reentry_stage)+ms(reentry_stage); % 20000kg structural mass + 5000kg remaining propellant mass, will be y(5)- M_other_stages - M_payload
A=pi*1^2; %total rocket : pi*3.66^2
Cd=0.85;

%7.0 waiting phase
phase=7.0;
param = [Isp(reentry_stage), Cd, A, reentry_stage, phase];

y0_reentry_0 = [Initial_speed Initial_gamma Initial_h Initial_x Initial_m]; % Initial state vector, 1 line

ti_reentry_0 = tf1; % Initial time for phase 7.1 (s) ie end of burnout of phase 1
tf_reentry_0 = 16+ti_reentry_0; % Final time for phase 7.1 (s). TBD

[t_reentry_0, y_reentry_0] = ode45(@(t, y) reentry_dynamicsODE(t, 0, y, param), (ti_reentry_0:0.05:tf_reentry_0), y0_reentry_0, options);

%7.1 - turning phase : thursting to turn to opposite gamma
end_gamma=pi+27.3*pi/180;%TBD, 5° entry for the moment
phase=7.1;

param = [Isp(reentry_stage), Cd, A, reentry_stage, phase];

y0_reentry_1 = [y_reentry_0(end, iV) y_reentry_0(end,igamma) y_reentry_0(end,ih) y_reentry_0(end,ix) y_reentry_0(end,im)]; % Initial state vector, 1 line

ti_reentry_1 = tf_reentry_0; % Initial time for phase 7.1 (s) ie end of burnout of phase 1
tf_reentry_1 = 13.1+ti_reentry_1; % Final time for phase 7.1 (s). TBD

[t_reentry_1, y_reentry_1] = ode45(@(t, y) reentry_dynamicsODE(t, -T(reentry_stage)/200, y, param, [Initial_gamma, end_gamma], [ti_reentry_1, tf_reentry_1]), (ti_reentry_1:0.05:tf_reentry_1), y0_reentry_1, options);

%7.2 going in the opposite direction
phase=7.2;
param = [Isp(reentry_stage), Cd, A, reentry_stage, phase];

y0_reentry_2 = [y_reentry_1(end, iV) y_reentry_1(end,igamma) y_reentry_1(end,ih) y_reentry_1(end,ix) y_reentry_1(end,im)]; % Initial state vector, 1 line

ti_reentry_2 = tf_reentry_1; % Initial time for phase 7.2, end of phase 7.1 (s)
tf_reentry_2 = 119.7+ti_reentry_2; % Final time for phase 7.2 (s). TBD, for the moment 1 minute

[t_reentry_2, y_reentry_2 ] = ode45(@(t, y) reentry_dynamicsODE(t, 0 , y, param), (ti_reentry_2:0.05:tf_reentry_2), y0_reentry_2, options);

%7.3 turning to gamma =pi/2
phase=7.3;
end_gamma=3*pi/2;
param = [Isp(reentry_stage), Cd, A, reentry_stage, phase];

y0_reentry_3 = [y_reentry_2(end, iV) y_reentry_2(end,igamma) y_reentry_2(end,ih) y_reentry_2(end,ix) y_reentry_2(end,im)]; % Initial state vector, 1 line

ti_reentry_3 = tf_reentry_2; % Initial time for phase 7.3, end of phase 7.2 (s)
tf_reentry_3 = 20+ti_reentry_3; % Final time for phase 7.3 (s). TBD, for the moment 10 s

[t_reentry_3, y_reentry_3] = ode45(@(t, y) reentry_dynamicsODE(t, -T(reentry_stage)/200, y, param, [y0_reentry_3(igamma), end_gamma], [ti_reentry_3, tf_reentry_3]), (ti_reentry_3:0.05:tf_reentry_3), y0_reentry_3, options);

%7.4 - braking phase: parachute to brake
phase=7.4;
Cd_parachute = 1.5; %Drag coefficient for a space parachute (kevlar/nylon stuff, deployed at approx 10km to avoid it to brake)
A_parachute = 250; %Surface of the 3 parachutes in contact with the airflow (m^2) TBD
param = [Isp(reentry_stage), Cd_parachute, A_parachute, reentry_stage, phase];

y0_reentry_4 = [y_reentry_3(end, iV) y_reentry_3(end,igamma) y_reentry_3(end,ih) y_reentry_3(end,ix) y_reentry_3(end,im)]; % Initial state vector, 1 line

ti_reentry_4 = tf_reentry_3; % Initial time for phase 7.4, end of phase 7.3 (s)
tf_reentry_4 = 60*4+ti_reentry_4; % Final time for phase 7.4 (s). TBD, for the moment 3 minutes


[t_reentry_4, y_reentry_4] = ode45(@(t, y) reentry_dynamicsODE(t, 0 , y, param), (ti_reentry_4:0.05:tf_reentry_4), y0_reentry_4, options);


%% Ploting phase
figure(10); hold on;
if add_ascent_plot
    plot(t1,y1(:,ih)/1e3,'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',2);
    plot(t2,y2(:,ih)/1e3,'MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',2);
    plot(t3,y3(:,ih)/1e3,'MarkerFaceColor','#EDB120','LineWidth',2);
end
plot(t_reentry_0,y_reentry_0(:,ih)/1e3,'m','LineWidth',2);
plot(t_reentry_1,y_reentry_1(:,ih)/1e3,'r','LineWidth',2);
plot(t_reentry_2,y_reentry_2(:,ih)/1e3,'MarkerFaceColor','#77AC30','LineWidth',2);
plot(t_reentry_3,y_reentry_3(:,ih)/1e3,'b','LineWidth',2);
plot(t_reentry_4,y_reentry_4(:,ih)/1e3,'c','LineWidth',2);
%plot(t_reentry_5,y_reentry_5(:,ih)/1e3,'m','LineWidth',2);
title(['Altitude change for reentry of stage ' num2str(reentry_stage)]);
xlabel('Time (s)');
ylabel('Altitude (km)');
grid;

% figure(11); hold on;
% plot(t_reentry_1,y_reentry_1(:,igamma)*180/pi,'r','LineWidth',2);
% plot(t_reentry_2,y_reentry_2(:,igamma)*180/pi,'g','LineWidth',2);
% plot(t_reentry_3,y_reentry_3(:,igamma)*180/pi,'b','LineWidth',2);
% plot(t_reentry_4,y_reentry_4(:,igamma)*180/pi,'c','LineWidth',2);
% %plot(t_reentry_5,y_reentry_5(:,igamma)*180/pi,'m','LineWidth',2);
% title(['Flight path angle change for reentry of stage ' num2str(reentry_stage)]);
% xlabel('Time (s)');
% ylabel('Flight path angle (deg)');
% grid;

figure(12); hold on;
if add_ascent_plot
    plot(y1(:,ix)/1e3,y1(:,ih)/1e3,'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',2);
    plot(y2(:,ix)/1e3,y2(:,ih)/1e3,'MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',2);
    plot(y3(:,ix)/1e3,y3(:,ih)/1e3,'MarkerFaceColor','#EDB120','LineWidth',2);
end
plot(y_reentry_0(:,ix)/1e3,y_reentry_0(:,ih)/1e3,'m','LineWidth',2);
plot(y_reentry_1(:,ix)/1e3,y_reentry_1(:,ih)/1e3,'r','LineWidth',2);
plot(y_reentry_2(:,ix)/1e3,y_reentry_2(:,ih)/1e3,'MarkerFaceColor','#77AC30','LineWidth',2);
plot(y_reentry_3(:,ix)/1e3,y_reentry_3(:,ih)/1e3,'b','LineWidth',2);
plot(y_reentry_4(:,ix)/1e3,y_reentry_4(:,ih)/1e3,'c','LineWidth',2);
%plot(y_reentry_5(:,ix)/1e3,y_reentry_5(:,ih)/1e3,'m','LineWidth',2);
title(['Altitude change for reentry of stage ' num2str(reentry_stage)]);
xlabel('X position (km)');
ylabel('Altitude (km)');
grid;

figure(13); hold on;
if add_ascent_plot
    plot(t1,y1(:,im),'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',2);
    plot(t2,y2(:,im),'MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',2);
    plot(t3,y3(:,im),'MarkerFaceColor','#EDB120','LineWidth',2);
end
plot(t_reentry_0,y_reentry_0(:,im),'m','LineWidth',2);
plot(t_reentry_1,y_reentry_1(:,im),'r','LineWidth',2);
plot(t_reentry_2,y_reentry_2(:,im),'MarkerFaceColor','#77AC30','LineWidth',2);
plot(t_reentry_3,y_reentry_3(:,im),'b','LineWidth',2);
plot(t_reentry_4,y_reentry_4(:,im),'c','LineWidth',2);
%plot(t_reentry_5,y_reentry_5(:,im),'m','LineWidth',2);
yline(ms(reentry_stage),'-','Reentry stage structural mass');
ylim([0.99*ms(reentry_stage),1.01*y_reentry_0(1,im)])
title(['mass change for reentry of stage ' num2str(reentry_stage)]);
xlabel('Time (s)'); 
ylabel('Mass (kg)');
grid;

figure(14); hold on;
if add_ascent_plot
    plot(t1,y1(:,iV),'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',2);
    plot(t2,y2(:,iV),'MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',2);
    plot(t3,y3(:,iV),'MarkerFaceColor','#EDB120','LineWidth',2);
end
plot(t_reentry_0,y_reentry_0(:,iV),'m','LineWidth',2);
plot(t_reentry_1,y_reentry_1(:,iV),'r','LineWidth',2);
plot(t_reentry_2,y_reentry_2(:,iV),'MarkerFaceColor','#77AC30','LineWidth',2);
plot(t_reentry_3,y_reentry_3(:,iV),'b','LineWidth',2);
plot(t_reentry_4,y_reentry_4(:,iV),'c','LineWidth',2);
%plot(t_reentry_5,y_reentry_5(:,iV),'m','LineWidth',2);
title(['velocity change for reentry of stage '  num2str(reentry_stage)]);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid;