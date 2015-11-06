%Battery Simulation Test
close all
clear all
clc

Battery_Capacity=100;
%Simulate Production=Consumption
%Should give a flat line at the batteries starting capacity or charge
%level
for i=1:10
    Power_Production_Vector(i)=10;
    Power_Consumption_Vector(i)=10;
end
%Simulate no production and consumption enough to drain the battery
%should give a linearly decreasing charge or capacity level until the
%battery becomes empty where it should show that all power is being drawn
%from the grid
for i=10:20
    Power_Production_Vector(i)=0;
    Power_Consumption_Vector(i)=30;
end
%Simulate battery charging with no load enough to fully charge the battery
%Should give a linearyly increasing charge of capacity level untill the
%battery becomes full where it should give a constant charge or capacity
%level at the batteries full capacity and show that the amount of power
%lost is increasing linearly
for i=20:35
    Power_Production_Vector(i)=30;
    Power_Consumption_Vector(i)=0;
end

for i=35:40
    Power_Production_Vector(i)=30;
    Power_Consumption_Vector(i)=15;
end

[ Energy_Level,Power_Lost,Power_Grid,Energy_Lost,Energy_Grid,Energy_Conv ] = Bat_Sim( Power_Production_Vector,Power_Consumption_Vector,Battery_Capacity );

figure
subplot(2,1,1)
plot(1:40,Power_Production_Vector,'g')
hold on
plot(1:40,Power_Consumption_Vector,'r')
plot(1:40,Power_Lost,'y','LineWidth',3)
plot(1:40,Power_Grid,'k','LineWidth',3)
title('Power Production and Consumption')
xlabel('Time[Days]')
ylabel('Power[W]')
legend('Power Production','Power Consumption','Power Lost','Power Grid')
subplot(2,1,2)
plot(1:39,Energy_Level)
title('Battery Level')
xlabel('Time[Days]')
ylabel('Battery Level[Wh]')
legend('Battery Level','Dumped Energy','Grid Energy')