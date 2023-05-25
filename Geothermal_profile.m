function [Temperatures,Solidification_time]=Geothermal_profile(T,Calculation_times,F_atm,F_mo)
% This function was created by Miguel Angel Carapia Perez for the paper 
% "Atmospheric influence on lithosphere formation during cooling of a 
% global magma ocean" in 2023.
% This function evaluates the heat equation in two dimensions for an Earth
% symmetric split in half. This Earth has thermal layers similar to today's
% Earth.
% 
% The input elements are:
% T:        A row vector that includes the temperatures of the [Core,Inner_Mantle,
%           Outer_Mantle, Atmosphere] Temperatures are entered in degrees Kelvins
% 
% Calculation_times: It is a vector that includes the calculation times for
%           the model. These calculation times are entered in millions of years. 
% 
% F_atm:    Heat flux that leaves the atmosphere.
% 
% F_mo:     Heat flux thal leaves the magma ocean.
% 
% The output elements is:
% Temperatures: They are the temperatures that are obtained directly 
%             from the calculations.
%
% Solidification_time: It is the time (in Myrs) in which the temperature 
%             drops below the solidification temperature of the peridotite 
%             (approx 1400 K).
% 
%-----------------------EXAMPLE------------------------------

% This is just an example that is returned by default in case the input 
% data has not been entered. It can easily be blocked or commented out 
% to work without it.
if nargin==0
close all
clear all
clc

disp('"This is only an example"')

% Calculation times 
Calculation_times=[0,0.1:0.1:0.9,1:1,9,10:10:90,100:50:900,1000:100:9000,...
    10000:500:90000,100000:10000:30000000];
Calculation_times=Calculation_times.*(1*10^-6);

% Temperatures of the Earth´s layers
TN1=5000; TM1=6000; TM2=8000; TA1=4000;
T=[TN1 TM1 TM2 TA1];

% Heat flux from magma ocean to atmosphere
F_mo=100000;

% Percentage of heat flux leaving the atmosphere as a function of heat flux 
% from the magma ocean, 1 is for 100%, 0.9 is for 90%,... 0.1 for 10%, etc. 
% The retention percentage is the remaining percentage, that is, if the 
% atmosphere lets 40% pass (q_atm=0.4) The percentage that the atmosphere 
% retains is 60%.
F_atm=1;

end
%-----------------------EXAMPLE------------------------------




%-----------------------GEOMETRY------------------------------

% Radii of Earth's layers (in meters)
Co= 3500000;    % Core
IM= 5500000;    % Inner mantle
OM= 6370000;    % Outer mantle 
A=  6800000;    % Atmosphere

% Geometry of the Earth with concentric spheres and a box behind forming the space
C=10000000;
cuadro1 = [2 
    4
    -C
    C
    C
    -C
    C-1000000
    C-1000000
    -C+1000000
    -C+1000000];
cuadro2 = [2 
    4
    -C
    0
    0
    -C
    C
    C
    -C 
    -C];
C1=[1
    0
    0
    Co];
C2=[1
    0
    0
    IM];
C3=[1
    0
    0
    OM];
C4=[1
    0
    0
    A];

C1 = [C1;zeros(length(cuadro1) - length(C1),1)];
C2 = [C2;zeros(length(cuadro1) - length(C2),1)];
C3 = [C3;zeros(length(cuadro1) - length(C3),1)];
C4 = [C4;zeros(length(cuadro1) - length(C4),1)];

gd = [C1,C2,C3,C4,cuadro1,cuadro2];
ns = char('C1','C2','C3','C4','cuadro1','cuadro2');
ns = ns';
sf = '(C1+C2+C3+C4)-cuadro2';

% The geometry of the model is created from the edges
[dl, bt] = decsg(gd,sf,ns);

%---------------------GEOMETRY------------------------------

%----------------------GRID---------------------------------

% Creation of the thermal model
thermalmodel = createpde('thermal','transient');
geometryFromEdges(thermalmodel,dl);

% figure(10);
% pdegplot(dl,"EdgeLabels","on","FaceLabels","on")
% axis equal

% Grid standar
mesh1=generateMesh(thermalmodel);

% The size of the grid cells is modified. These values were obtained 
% with trial and error, so I recommend not modifying them, or doing it consciously.
Tam_Max=608210;
c=ceil(Tam_Max/20)+1;

% The number of nodes at the edges of the model is highlighted, especially 
% for the edges between the outermost layers of the Earth where there is 
% the largest initial temperature difference.
mesh1=generateMesh(thermalmodel,'GeometricOrder','linear','Hedge',{[8 9 10 11 12 13 14 15],c});
disp('Mesh')
pause(1)

[p,e,t] = meshToPet(mesh1);
% figure(30)
% pdemesh(p,e,t)
% pause(3)

[p,e,t] = refinemesh(dl,p,e,t,[1 2 3 4]);
disp('First refined')
pause(1)
[p,e,t] = refinemesh(dl,p,e,t,[1 2 3 4]);
disp('Second refined')
pause(1)
[p,e,t] = refinemesh(dl,p,e,t,[2 3 4]);
disp('Third refined')
pause(1)
[p,e,t] = refinemesh(dl,p,e,t,[2 3 4]);
disp('Fourth refined')
pause(1)
[p,e,t] = refinemesh(dl,p,e,t,[2 3]);
disp('Fifth refined')
pause(1)

%-----------------------------GRID------------------------------



%---------------------THERMAL PROPERTIES--------------------------


% Earth and atmosphere surface. 
Surf_Earth=2*pi*OM;
Surf_Atm=2*pi*A;
% Heat flux from the magma ocean distributed on the surface of the atmosphere.
F_atm2=(F_mo*Surf_Earth)/Surf_Atm;
k_space=1*10^8;


% This was done in order to model conduction cooling. 
% In this case the heat flow from the atmosphere is no longer a function 
% of the magma ocean and the cooling of the atmosphere can be seen more clearly.
if F_mo==0
k_eq_mo=3.6;
k_eq_atm=F_atm*(A-OM)/T(4);
else
k_eq_mo=F_mo*(A-OM)/T(4);
k_eq_atm=F_atm*F_atm2*(A-OM)/T(4);

if k_eq_mo>k_space
k_eq_mo=k_space;
else
end

if k_eq_atm>k_space
k_eq_atm=k_space;
else
end

end


% Core
thermalProperties(thermalmodel,'Face',4,'ThermalConductivity',60, ...   % Melchior(1986) 
                                        'MassDensity',11000, ...
                                        'SpecificHeat',1400);
% Inner mantle
thermalProperties(thermalmodel,'Face',1,'ThermalConductivity',k_eq_mo, ... %%%%%%%%%%%%%%%%
                                        'MassDensity',3500, ...
                                        'SpecificHeat',700);
% Outer mantle
thermalProperties(thermalmodel,'Face',2,'ThermalConductivity',k_eq_mo, ...%%%%%%%%%%%%%%%%%%%%
                                        'MassDensity',3500, ...
                                        'SpecificHeat',700);
% Atmosphere
thermalProperties(thermalmodel,'Face',3,'ThermalConductivity',k_eq_atm, ...%%%%%%%%%%%%%%%%%%
                                        'MassDensity',1.293, ...
                                        'SpecificHeat',1004); % ATMOSPHERE, OCEAN, AND CLIMATE DYNAMICS: AN INTRODUCTORY TEXT

thermalmodel.StefanBoltzmannConstant = 5.670367e-8;

% % Temperatures
% These temperatures are entered before running the function.
TC=T(1); TIM=T(2); TOM=T(3); TA=T(4);

% Initial thermal conditions of each layer 
thermalIC(thermalmodel,TC,'Face',4)
thermalIC(thermalmodel,TIM,'Face',1)
thermalIC(thermalmodel,TOM,'Face',2)
thermalIC(thermalmodel,TA,'Face',3)

% Thermal boundary conditions that establish the value of 0°K for outer space.
thermalBC(thermalmodel,'Edge',[14 15],'Temperature',0);

% Radiogenic heat source
internalHeatSource(thermalmodel,@heatSource,'Face',[1,2]);

% Calculation times (in seconds) in the function are entered in millions of years
My_to_s=1000000*365*24*60*60;
tlist=Calculation_times.*My_to_s;
TMa=Calculation_times;

x = linspace(0,C,1000);
y = zeros(1,length(x));

pause(1)
disp('Started calculations')

% Solution of the thermal model
R = solve(thermalmodel,tlist);  

% Interpolation for the thermal profile of the whole Earth
Tintrp = interpolateTemperature(R,x,y,1:length(tlist));

pause(1)
disp('finished calculations')


% Time of solidification at the surface
for i=1:length(Calculation_times) 

    if Tintrp(637,i)<1400
        Solid_time(i)=TMa(i);
    else
        Solid_time(i)=NaN;
        n(i)=length(Solid_time);
    end

end

% Time of solidification at the surface
Solidification_time1=min(Solid_time);
%---------------------------Graphics--------------------------

Prof=0:1:2870;
[T_sol,T_liq]=T_sol_liq(Prof);
Profu=6370-Prof;

NER=length(Tintrp(:,1));
deepT=linspace(0,C/1000,NER);
f1=figure('color','white');
hold on 
for i=1:1:length(n)+1
axis([3500 7000 0 10000])
grid on
plot(Profu,T_sol,'k-',LineWidth=2)
plot(Profu,T_liq,'k-',LineWidth=2)
plot(deepT,Tintrp(:,i))
xlabel('Earth`s radio(Km)','FontSize',12),ylabel('Temperature (K)','FontSize',12)
title(['Temperature at Time ' num2str(TMa(i)) ' Myrs'],'FontSize',16)
%    pause(1)
end

%---------------------------Graphics--------------------------

%------------------------Data Ouput---------------------------

Solidification_time=Solidification_time1; % In million years
Temperatures=Tintrp; % In Kelvin
disp('All completed')
%------------------------Data Ouput---------------------------

%------------------------Extra Functions---------------------------

function Q = heatSource(location,state)
  Q = zeros(1,numel(location.x));
t0=(4.5*10^9)*365*24*60*60;
s_a=3.168808781402895*10^-8;
%CONCENTRATION  
       U238_C=0.9928*0.20;     U235_C=0.0072*0.20;     Th232_C=0.069;
       K40_C=(1.17*10^-4)*270;      AL26_C=(5*10^-5)*8650;
%HEAT PRODUCTION
       U238_H=9.17*10^-5;     U235_H=5.75*10^-4;     Th232_H=2.56*10^-5;
        K40_H=2.97*10^-5;      AL26_H=3.54*10^-1;
%LAMBDA        
        U238=(4.19*10^-18);     U235=(3.12*10^-17);     Th232=(1.56*10^-18);
        K40=(1.72*10^-17);      AL26=(3.06*10^-14);

if(isnan(state.time))
% Returning a NaN when time=NaN tells the solver that the heat source is a function of time.
  Q(1,:) = NaN;
  return
end

if state.time > 1
  Q(1,:) = ((U238_C.*U238_H.*exp(-U238*((state.time*s_a)-t0))...
        +U235_C*U235_H.*exp(-U235.*((state.time*s_a)-t0))...
         +Th232_C.*Th232_H.*exp(-Th232.*((state.time*s_a)-t0))...
         +K40_C.*K40_H.*exp(-K40.*((state.time*s_a)-t0))) ...
         +AL26_C.*AL26_H.*exp(-AL26.*state.time))*10^-6;
end

end
%------------------------Extra Functions---------------------------


 
end
