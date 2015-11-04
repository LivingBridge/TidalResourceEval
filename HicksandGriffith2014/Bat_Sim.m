function [ Charge_Level ] = Bat_Sim( Power_Production_Vector,Power_Consumption_Vector,Battery_Capacity )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

Battery_Start_Capacity=0.5*Battery_Capacity;    %[Wh] Assumed starting capacity
Circuit_Voltage=120;                            %[V-AC] Assumed...this could be different in reality depending on circuit design. And could also be different coming from turbine vs coming off alternator to the load.


%Detailed Battery Parameters
%Battery Profile(This assumes linear profile between capacity[Ah] and charge level[V], I know this is not realistic)
N=100;%Number of profile points
Max_Voltage=450;%[V] Battery full voltage found on http://www.teslamotors.com/powerwall
Min_Voltage=350;%[V] Battery empty voltage found on http://www.teslamotors.com/powerwall
Charge_Prof=linspace(Min_Voltage,Max_Voltage,N);%[V]
Capacity_Prof=linspace(0,Battery_Capacity,N);   %[Wh]


Energy_Consumption_Vector=Power_Consumption_Vector/4;%The amount of [Wh] consumed in a vector where each element represents 15 minutes
Current_Load_Vector=Energy_Consumption_Vector/Circuit_Voltage;%This is the amount of current[A] that the load is drawing assuming everything is running on a 120 V circuit.
Energy_Production_Vector=Power_Production_Vector/4;%The amount of [Wh] produced in a vector where each element represents 15 minutes
Current_Source_Vector=Energy_Production_Vector/Circuit_Voltage;%This is the amount of current[A] that the turbine is producing assuming everything is running on a 120 V circuit.

j=1;
k=1;
l=1;
% Capacity_Level=zeros(length(Power_Production_Vector));
% Charge_Level
for i=2:length(Power_Production_Vector-1)
    Capacity_Level(1)=Battery_Start_Capacity; %[Wh]
    Charge_Level(1)=interp1(Capacity_Prof,Charge_Prof,Battery_Start_Capacity); %[V]
    Capacity_Level(i)=Capacity_Level(i-1)+(Current_Source_Vector(i)-Current_Load_Vector(i));%subract out Ah consumed by load add in Ah produced by turbine
    Charge_Level(i)=interp1(Capacity_Prof,Charge_Prof,Capacity_Level(i));
    %Battery becomes full ATS=Open, Break=On
    if Charge_Level(i)>Max_Voltage
        Current_Source_Vector(i+1)=0;%Break is on, thus no current comes from turbine to battery
        Power_Lost(i+1)=Power_Production_Vector(i+1);%The amount of power lost is equal to the amount of power that could have been generated
        j=j+1;
    end
    %Battery becomes empty ATS=Closed, Break=Off
    if Charge_Level(i)<Min_Voltage
        Current_Load_Vector(i+1)=0;%Current draw from the battery goes to zero because the load is drawing from the grid
        Power_Grid(i+1)=Power_Consumption_Vector(i+1);%All power consumption comes off the grid
        k=k+1;
    end
    %Battery is in between empty and full ATS=Open, Break=Off
    if Charge_Level(i)>Min_Voltage && Charge_Level(i)<Max_Voltage
        Power_Grid(i+1)=0;%No power is drawn from the grid
        Power_Lost(i+1)=0;%No power is wasted as the turbine spins and all possible power is stored in the battery and then consumed
        l=l+1;
    end

end

end
