% This script was created by Miguel Angel Carapia Perez for the paper 
% "Atmospheric influence on lithosphere formation during cooling of a 
% global magma ocean" in 2023.
% This small code was made to reproduce the figures of the manuscript
% using the Geothermal_profile function. 
% This is only as an example to show how the function can be used

close all
clear all
clc

% Calculation times 
Calculation_times=[0,0.1:0.1:0.9,1:1,9,10:10:90,100:50:900,1000:100:9000,...
    10000:500:90000,100000:10000:30000000]; % Years
Calculation_times=Calculation_times.*(1*10^-6); % Million years

% Percentage of heat flux leaving the atmosphere as a function of heat flux 
% from the magma ocean, 1 is for 100%, 0.9 is for 90%,... 0.1 for 10%, etc. 
% The retention percentage is the remaining percentage, that is, if the 
% atmosphere lets 40% pass (q_atm=0.4) The percentage that the atmosphere 
% retains is 60%.
q_atm = [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1];

% Magma ocean heat flux
% q_mo  = 0; % This case is for heat flux with only conduction  
q_mo  = 100000;

% Temperatures
t1=[5000 6000 8000 4000 0]; % Profile of maximum temperatures
t2=[4000 3500 5000 2000 0]; % Profile of minimum temperatures
T=t1; % Choice of temperature for current modeling

% Empty matrices to facilitate calculations
T_Lithos=zeros(length(q_mo),length(q_atm));
Temper_1=zeros(length(Calculation_times),length(q_mo));

% Calculation cycles that match each heat flow from the atmosphere with 
% each heat flow from the magma ocean
for i=1:length(q_mo)
    for j=1:length(q_atm)
q_mo_t=q_mo(i)
q_atm_t=q_atm(j)

[Temperatures_1,T_Lithos_1]=Geothermal_profile(T,Calculation_times,q_atm_t,q_mo_t);

    if length(q_mo) < length(q_atm)
Temper_1(:,j)=Temperatures_1(637,:);
T_Lithos(i,j)=T_Lithos_1(:);
    else
Temper_1(:,i)=Temperatures_1(637,:);
T_Lithos(i,j)=T_Lithos_1(:);
    end
% % %     % Uncomment for figure 6b
T_0y(j,:)=Temperatures_1(:,1);
T_7300y(j,:)=Temperatures_1(:,102);
T_8600y(j,:)=Temperatures_1(:,115);
T_10500y(j,:)=Temperatures_1(:,121);
T_13000y(j,:)=Temperatures_1(:,126);
T_16500y(j,:)=Temperatures_1(:,133);
T_21500y(j,:)=Temperatures_1(:,143);
T_28500y(j,:)=Temperatures_1(:,157);
T_40000y(j,:)=Temperatures_1(:,180);
T_61000y(j,:)=Temperatures_1(:,222);
T_130000y(j,:)=Temperatures_1(:,284);

% % %     % Uncomment for figure 6a
% T_0y(j,:)=Temperatures_1(:,1);
% T_1500y(j,:)=Temperatures_1(:,44);
% T_1700y(j,:)=Temperatures_1(:,46);
% T_2100y(j,:)=Temperatures_1(:,50);
% T_2600y(j,:)=Temperatures_1(:,55);
% T_3200y(j,:)=Temperatures_1(:,61);
% T_4300y(j,:)=Temperatures_1(:,72);
% T_6100y(j,:)=Temperatures_1(:,90);
% T_10000y(j,:)=Temperatures_1(:,120);
% T_17000y(j,:)=Temperatures_1(:,134);
% T_37500y(j,:)=Temperatures_1(:,175);

    end
end

% Melting curves, 
Tsolid=ones(1,length(Calculation_times)).*1400;
Prof=0:1:2870;
[T_sol,T_liq]=T_sol_liq(Prof);
Profu=6370-Prof;
deep=linspace(0,10000,1000);
Solidification_times=T_Lithos.*1000000


%%%%%% Figure 3a, 3b
% f1=figure('Color','white');
% axis([0 30 0 7000])
% hold on
% grid on
% plot(Calculation_times,Temper_1(:,1),'k-',LineWidth=2)
% plot(Calculation_times,Temper_1(:,2),'k:',LineWidth=2)
% plot(Calculation_times,Temper_1(:,3),'k.',LineWidth=2)
% plot(Calculation_times,Temper_1(:,4),'k-.',LineWidth=2)
% plot(Calculation_times,Temper_1(:,5),'g',LineWidth=2)
% plot(Calculation_times,Temper_1(:,6),'b',LineWidth=2)
% plot(Calculation_times,Temper_1(:,7),'r',LineWidth=2)
% plot(Calculation_times,Tsolid,'m--',LineWidth=2)
% ylabel('Surface temperature (Kº)','FontSize',12)
% xlabel('Time (Myrs)','FontSize',12)
% legend('0.1 W/m^2','0.2 W/m^2','0.5 W/m^2','0.7 W/m^2','1 W/m^2',...
%     '10 W/m^2','100 W/m^2','Temperature solidus','location','northeast')
% title('Conduction cooling (Max Temperatures)','FontSize',16) 

%%%%%% Figure 4a, 4b
% f2=figure('Color','white');
% axis([0 1000000 0 7000])
% hold on
% grid on
% plot(Calculation_times*1000000,Temper_1(:,1),'r',LineWidth=2)
% plot(Calculation_times*1000000,Temper_1(:,2),'g',LineWidth=2)
% plot(Calculation_times*1000000,Temper_1(:,3),'b',LineWidth=2)
% plot(Calculation_times*1000000,Temper_1(:,4),'c',LineWidth=2)
% plot(Calculation_times*1000000,Tsolid,'k--',LineWidth=2)
% ylabel('Surface temperature (Kº)','FontSize',12)
% xlabel('Time (Years)','FontSize',12)
% legend('1,000 W/m^2','10,000 W/m^2','100,000 W/m^2','1,000,000 W/m^2',...
%     'Temperature solidus','location','northeast')
% title('Heat flux from magma ocean convection','FontSize',16) 
% 

%%%%% Figure 6a
% f6=figure('Color','white');
% axis([3500 7000 0 8000])
% hold on
% grid on
% plot(deep,T_0y(1,:),'k',LineWidth=2);
% plot(deep,T_1500y(1,:),'b:',LineWidth=2);
% plot(deep,T_1700y(2,:),'b--',LineWidth=2);
% plot(deep,T_2100y(3,:),'b-.',LineWidth=2);
% plot(deep,T_2600y(4,:),'b-',LineWidth=2);
% plot(deep,T_3200y(5,:),'g:',LineWidth=2);
% plot(deep,T_4300y(6,:),'g--',LineWidth=2);
% plot(deep,T_6100y(7,:),'g',LineWidth=2);
% plot(deep,T_10000y(8,:),'r--',LineWidth=2);
% plot(deep,T_17000y(9,:),'r-.',LineWidth=2);
% plot(deep,T_37500y(10,:),'r-',LineWidth=2);
% plot(Profu,T_sol,'k',LineWidth=2)
% plot(Profu,T_liq,'k',LineWidth=2)
% ylabel('Temperature (Kº)','FontSize',12)
% xlabel('Earth`s radio','FontSize',12)
% legend('Initial conditions',...
%     '1,500 yrs (0 %)','1,700 yrs (10 %)','2,100 yrs (20 %)','2,600 yrs (30 %)',...
%     '3,200 yrs (40 %)','4,300 yrs (50 %)','6,100 yrs (60 %)','10,000 yrs (70 %)',...
%     '17,000 yrs (80 %)','37,500 yrs (90 %)',...
%     'Temperature solidus','Temperature liquidus','location','northeast')

%%%% Figure 6b 
f5=figure('Color','white');
axis([3500 7000 0 8000])
hold on
grid on
plot(deep,T_0y(1,:),'k',LineWidth=2);
plot(deep,T_7300y(1,:),'b:',LineWidth=2);
plot(deep,T_8600y(2,:),'b--',LineWidth=2);
plot(deep,T_10500y(3,:),'b-.',LineWidth=2);
plot(deep,T_13000y(4,:),'b-',LineWidth=2);
plot(deep,T_16500y(5,:),'b:',LineWidth=2);
plot(deep,T_21500y(6,:),'g--',LineWidth=2);
plot(deep,T_28500y(7,:),'g',LineWidth=2);
plot(deep,T_40000y(8,:),'r--',LineWidth=2);
plot(deep,T_61000y(9,:),'r-.',LineWidth=2);
plot(deep,T_130000y(10,:),'r-',LineWidth=2);
plot(Profu,T_sol,'k',LineWidth=2)
plot(Profu,T_liq,'k',LineWidth=2)
ylabel('Temperature (Kº)','FontSize',12)
xlabel('Earth`s radio','FontSize',12)
legend('Initial conditions',...
    '7,300 yrs (0 %)','8,600 yrs (10 %)','10,500 yrs (20 %)','13,000 yrs (30 %)',...
    '16,500 yrs (40 %)','21,500 yrs (50 %)','28,500 yrs (60 %)','40,000 yrs (70 %)',...
    '61,000 yrs (80 %)','130,000 yrs (90 %)',...
    'Temperature solidus','Temperature liquidus','location','northeast')
% title('Temperature at the formation of lithosphere','FontSize',16) 



