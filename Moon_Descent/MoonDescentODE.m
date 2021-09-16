%% Gunnar Tibert, KTH, 2020-04-20

function dy = MoonDescentODE(t,y,T,gamma,tf)
dy = zeros(5,1);    % a column vector

R_M=1737e3; %% mean radius of the Moon (m)
Isp=321; %% specific impulse of engine (s)
g=1.625; %% gravitational acceleration on Moon surface (m/s^2)
g0=9.81; %% gravitational acceleration on Earth surface (m/s^2)

dy(1)=T/y(5)-g*sin(y(2)); %% acceleration in spacecraft system (m/s^2)
dy(2)=(gamma(2)-gamma(1))/(tf(2)-tf(1)); %% flight path angle (1/s)
dy(3)=y(1)*sin(y(2)); %% altitude rate (m/s)
dy(4)=y(1)*cos(y(2))*R_M/(R_M+y(3)); %% ground distance rate (m/s)
dy(5)=-abs(T)/Isp/g0; %% Propellant mass flow rate (kg/s)


