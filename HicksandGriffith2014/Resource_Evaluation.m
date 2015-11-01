clear all; close all; clc
A=4.5;              %Turbine Area[m^2](1 for power density)
S=29.8;             %Average Salinity[PPT] over a year(2004) at Great Bay Coastal Marine Lab (greatbaydata.org)
T=8.749;            %Average Temperature[C] over a year(2004) at Great Bay Coastal Marine Lab (greatbaydata.org)
sample_period=15;   %minutes(check this)
cut_in=0.7;         %[m/s]
Cp=0.35;            %Efficiency of turbine
bin_space=0.3;      %[m]

Data=load('ADCP.mat');
Other_Data=load('MEMBRENG.CNV');
ADCP_Depth=Other_Data(:,11)/10;
c=size(Data.adcp);
for i=1:c(2)
    Time{i}=Data.adcp(:,i,7);
    Velocity{i}=Data.adcp(:,i,9);
    Speed{i}=Data.adcp(:,i,13);
    for j=2592:5472 %Cut out month where constuction barge was above from day 27 to day 57
        Speed{i}(j)=NaN;
    end
    Direction{i}=degtorad(Data.adcp(:,i,14)+16.25);
    Bin_Height{i}=Data.adcp(:,i,8);    %approximate depth below water surface[m]
    rho(i)=seadens(T,S,ADCP_Depth(i));
    Power{i} = (1/2)*A*rho(i)*Speed{i}.^3;
    for j=1:length(Speed{i})
        if Speed{i}(j) < cut_in
            Power_lim{i}(j,1)=0;
        else
            Power_lim{i}(j,1)=Cp*Power{i}(j,1);
        end
    end
end
for i=1:c(2)
    max_current_bin(i)=max(Speed{i});
end
max_current=max(max_current_bin);
Time_Increment=sample_period*linspace(1,c(1),c(1));
Time_Increment_Days=Time_Increment/1440;
% Determine peaks and troughs of power generation to descretize tidal cycles

for i=1:c(2)
    [Speed_Smoothed{i},L] = wsmooth(Speed{i},Time_Increment,2);
    [maxtab{i}, mintab{i}]=peakdet(Speed_Smoothed{i}, .1 , Time_Increment);
end

%Finding zero crossings to determine amount of time tides flood and ebb for
for i=5:c(2)
    k=1;
    for j=1:length(Velocity{1,i})-1
        if Velocity{1,i}(j,1)*Velocity{1,i}(j+1,1)<0
            Zero_Crossing_Indexs{1,i}(k)=j*15; %These are the times(in minutes(15)) of slack water
            k=k+1;
        end
    end
    for l=1:length(Zero_Crossing_Indexs{1,i})-1
        Time_Btwn_Tides{1,i}(l)=Zero_Crossing_Indexs{1,i}(1,l+1)-Zero_Crossing_Indexs{1,i}(1,l);
        ebb_durations{i}=Time_Btwn_Tides{i}(1:2:end);
        flood_durations{i}=Time_Btwn_Tides{i}(2:2:end);
    end
end

Mean_Ebb_Duration=mean(ebb_durations{c(2)});
Mean_Flood_Duration=mean(flood_durations{c(2)});

%Use Figure 1 in "Histograms of Energy in Tidal Cycle" to check that peaks
%and troughs are determined correctly
%Statistics about "Index_Diff"(the times between what the program believes
%to be one half cycle) can also be used for this purpose

for i=5:c(2)
    k=1;
    for j=1:length(Time_Increment)
        if Time_Increment(j)==mintab{i}(k,1)
            index{i}(k)=j;
            k=k+1;
            if k>length(mintab{i})
                break
            end
        end
    end
    for l=1:length(index{i})-1
        %energy in Watt-minutes
        energy{i}(l)=trapz(Time_Increment(index{i}(l):index{i}(l+1)),Power{i}(index{i}(l):index{i}(l+1)));
        energy{i}(l)=energy{i}(l)/60;%convert to Watt-hr
    end
    ebb_energy{i}=energy{i}(1:2:end);
    flood_energy{i}=energy{i}(2:2:end);
end

for i=5:c(2)
    for j=1:length(index)-1
        Index_Diff{i}(j)=(index{i}(j+1)-index{i}(j))*15/60;
    end
end



% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% ::::::::::::::::::: Code Used for Generation of Plots :::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


% -------------------------------------------------------------------------
% --------------           Speed vs Time Plots               --------------
% -------------------------------------------------------------------------


% % for i=1:c(2)
% %     figure(i)
% %     plot(Time_Increment_Days,Speed{i})
% %     datetick('x','mmm-dd','keepticks')
% %     NumTicks = 300;
% %     xlabel('Time[]')
% %     ylabel('Speed[m/s]')
% %     title(strcat('Current Speeds at:',num2str(Bin_Height{1,i}(1,1)),'m Above Seafloor'))
% %     pause
% %     close all
% % end

% -------------------------------------------------------------------------
% --------------    Flow Power vs Turbine Power Plots        --------------
% -------------------------------------------------------------------------

% % %29.52 days from start increment to end increment
% % start_increment=5472;
% % end_increment=8306;
% % for i=5:c(2)
% %     figure(i)
% %     plot(Time_Increment_Days(start_increment:end_increment),Power{i}(start_increment:end_increment))
% %     hold on
% %     plot(Time_Increment_Days(start_increment:end_increment),Power_lim{i}(start_increment:end_increment),'-g')
% % %     datetick('x','mmm-dd','keepticks')
% % %     NumTicks = 300;
% %     xlabel('Time[Days]')
% %     ylabel('Power Density[W/m^2]')
% %     legend('Current','Theoretical Turbine')
% %     title(strcat('Flow Power Density vs Turbine Power Density at:',num2str(height(i)),'m Above Seafloor'))
% %     pause
% %     close all
% % end

% -------------------------------------------------------------------------
% --------------      Histogram of Speed Measurements        --------------
% -------------------------------------------------------------------------


% % for i=1:c(2)
% %     figure(i)
% %     bin_centers=[1000:1000:12000];
% %     hist(Speed{i},2*(length(Speed{i})^(1/3)));
% %     title(strcat('Current Speeds at:',num2str(height(i)),'m Above Seafloor'))
% %     ylabel('Frequency')
% %     xlabel('Speed[m/s]')
% %     xlim([0,2.5])
% %     ylim([0,700])
% %     pause
% %     close all
% % end

% -------------------------------------------------------------------------
% --------------     Histograms of Energy in Tidal Cycle      -------------
% -------------------------------------------------------------------------

% % figure(1)
% % plot(Time_Increment,Speed_Smoothed{5})
% % plot(maxtab{5}(:,1), maxtab{5}(:,2), 'r*');
% % plot(mintab{5}(:,1), mintab{5}(:,2), 'g*');
% % xlabel('Time[min]')
% % ylabel('Speed[m/s]')
% % title('Cycle Descretization')
% % legend('Speed','Peak','Trough')
% % 
% % for i=1:c(2)
% %     figure(i+1)
% %     title(strcat('Depth:',num2str(height(i))))
% %     subplot(3,1,1);
% %     bin_centers=[1000:1000:16000];
% %     hist(energy{i},bin_centers);
% %     title(strcat('Energy per Half Tidal Cycle 1m^2 Turbine at:',num2str(height(i)),'m Above Seafloor'))
% %     ylabel('Frequency')
% %     xlabel('Energy per Cycle[Wh]')
% %     xmin=0;
% %     xmax=16000;
% %     ymin=0;
% %     ymax=150;
% %     xlim([xmin,xmax])
% %     ylim([ymin,ymax])
% %     
% %     subplot(3,1,2);
% %     hist(ebb_energy{i},bin_centers);
% %     title(strcat('Energy per Ebb Tidal Cycle 1m^2 Turbine at:',num2str(height(i)),'m Above Seafloor'))
% %     ylabel('Frequency')
% %     xlabel('Energy per Cycle[Wh]')
% %     xlim([xmin,xmax])
% %     ylim([ymin,ymax])
% %     
% %     subplot(3,1,3);
% %     hist(flood_energy{i},bin_centers);
% %     title(strcat('Energy per Flood Tidal Cycle 1m^2 Turbine at:',num2str(height(i)),'m Above Seafloor'))
% %     ylabel('Frequency')
% %     xlabel('Energy per Cycle[Wh]')
% %     xlim([xmin,xmax])
% %     ylim([ymin,ymax])
% %     pause
% %     close all
% % end

% -------------------------------------------------------------------------
% ------------   Polar Plots of Speed and Current Direction    ------------
% -------------------------------------------------------------------------

% % for i=1:c(2)
% %     figure
% %     t = 0 : .01 : 2 * pi;
% %     P = polar(t, 2.5 * ones(size(t)));
% %     set(P, 'Visible', 'off')
% %     hold on
% %     polar(Direction{i},Speed{i},'.b')
% %     set(gca,'color','none')
% %     view([90 -90])
% %     title({strcat('Current Speed[m/s] and Direction at Approximately',num2str(height(i)),'m Above Seafloor'),''})
% %     pause
% %     close all
% % end

% -------------------------------------------------------------------------
% --Plot Comparing Energy Yeild vs Turbine Size, Efficieny, and Cut-In ----
% -------------------------------------------------------------------------

% % Diameter=4;             %[m]
% % cut_in=0.4:0.1:1.5;     %[m/s]
% % Turbine_Length=2;       %[m]
% % for i=1:length(cut_in)
% %     [ Thirty_Day_Energy(i) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in(i) );
% % end
% % plot(cut_in,Thirty_Day_Energy)
% % hold on
% % Turbine_Length=3;       %[m]
% % for i=1:length(cut_in)
% %     [ Thirty_Day_Energy(i) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in(i) );
% % end
% % plot(cut_in,Thirty_Day_Energy)
% % Turbine_Length=4;       %[m]
% % for i=1:length(cut_in)
% %     [ Thirty_Day_Energy(i) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in(i) );
% % end
% % plot(cut_in,Thirty_Day_Energy)
% % xlabel('Cut In[m/s]')
% % ylabel('30 Day Energy Production[Wh]')
% % legend('8 m^2 Turbine Cp=0.35','12 m^2 Turbine Cp=0.35','16 m^2 Turbine Cp=0.35') 


% -------------------------------------------------------------------------
% --Plot Comparing Energy Yeild vs Turbine Size, Efficieny, and Cut-In ----
% -------------------------------------------------------------------------

% % Diameter=4;         %[m]
% % Cp=0.45;
% % cut_in=1;           %[m/s]
% % for Turbine_Length=2:5
% %     [ Thirty_Day_Energy(Turbine_Length-1) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % figure(1)
% % hold on
% % plot((Diameter*(2:5)),Thirty_Day_Energy,'--m')
% % cut_in=0.9;
% % for Turbine_Length=2:5
% %     [ Thirty_Day_Energy(Turbine_Length-1) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(2:5)),Thirty_Day_Energy,':m')
% % cut_in=0.5;
% % for Turbine_Length=2:5
% %     [ Thirty_Day_Energy(Turbine_Length-1) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(2:5)),Thirty_Day_Energy,'-m')
% % Cp=0.35;
% % cut_in=1;
% % for Turbine_Length=2:5
% %     [ Thirty_Day_Energy(Turbine_Length-1) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(2:5)),Thirty_Day_Energy,'--g')
% % cut_in=0.9;
% % for Turbine_Length=2:5
% %     [ Thirty_Day_Energy(Turbine_Length-1) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(2:5)),Thirty_Day_Energy,':g')
% % cut_in=0.5;
% % for Turbine_Length=2:5
% %     [ Thirty_Day_Energy(Turbine_Length-1) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(2:5)),Thirty_Day_Energy,'-g')
% % Cp=0.25;
% % cut_in=1;
% % for Turbine_Length=2:5
% %     [ Thirty_Day_Energy(Turbine_Length-1) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(2:5)),Thirty_Day_Energy,'--k')
% % cut_in=0.9;
% % for Turbine_Length=2:5
% %     [ Thirty_Day_Energy(Turbine_Length-1) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(2:5)),Thirty_Day_Energy,':k')
% % cut_in=0.5;
% % for Turbine_Length=2:5
% %     [ Thirty_Day_Energy(Turbine_Length-1) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(2:5)),Thirty_Day_Energy,'-k')
% % xlabel('Turbine Swept Area[m^2]')
% % ylabel('Energy Produced in 30 Days[Wh]')
% % legend('Cp=0.45, Cut In=1m/s','Cp=0.45, Cut In=0.9m/s','Cp=0.45, Cut In=0.5m/s','Cp=0.35, Cut In=1m/s','Cp=0.35, Cut In=0.9m/s','Cp=0.35, Cut In=0.5m/s','Cp=0.25, Cut In=1m/s','Cp=0.25, Cut In=0.9m/s','Cp=0.25, Cut In=0.5m/s','Location','NorthWest')

% -------------------------------------------------------------------------
% --------  Theoretical Energy Production and Consumption Plot  -----------
% -------------------------------------------------------------------------

% % start_index=5473;                       %Start at end of barge over ADCP bad data
% % A=7;                                   %[m^2]
% % E_bat=10000;                           %[Wh]
% % Energy_Consumption_Per_Month=900000;    %[Wh]
% % Power_Consumption=Energy_Consumption_Per_Month/(24*30); %[W]
% % Power_Consumption_Vector=Power_Consumption*ones(1,length(Power_lim{10}(start_index:end)));
% % Energy_Consumption_Vector=Power_Consumption_Vector/4;
% % figure
% % subplot(2,1,1)
% % plot(Time_Increment_Days(start_index:end),A*Power_lim{10}(start_index:end))
% % hold on
% % plot(Time_Increment_Days(start_index:end),Power_Consumption_Vector,'r')
% % title('Power Production and Consumption')
% % xlabel('Time[Days]')
% % ylabel('Power[W]')
% % legend('Power Production','Power Consumption')
% % j=1;
% % k=1;
% % l=1;
% % %This might be ok or it might mess things up
% % Under_indeces=zeros(length(Time_Increment_Days(start_index:end)));
% % Over_indeces=zeros(length(Time_Increment_Days(start_index:end)));
% % Middle_indeces=zeros(length(Time_Increment_Days(start_index:end)));
% % %To here
% % for i=2:length(Power_lim{10}(start_index:end))-1
% %     E(1)=E_bat;
% %     E_test(1)=E_bat;
% %     E_test(i)=E_test(i-1)+(-Energy_Consumption_Vector(i)+(A*Power_lim{10}(start_index+i)/4));
% %     E(i)=E(i-1)+(-Energy_Consumption_Vector(i)+(A*Power_lim{10}(start_index+i)/4));
% %     if E(i)>E_bat;
% %         E_over(i)=A*Power_lim{10}(start_index+i)/4;
% %         E(i)=E(i-1)-Energy_Consumption_Vector(i);
% %         Over_indeces(j)=i;
% %         j=j+1;
% %     end
% %     if E(i)<0;
% %         E_under(i)=-E(i);
% %         E_over(i)=0;
% %         Under_indeces(k)=i;
% %         k=k+1;
% %     end
% %     if E(i)<E_bat && E(i)>0
% %         E_under(i)=0;
% %         Middle_indeces(l)=i;
% %         l=l+1;
% %         E_over(i)=0;
% %     end
% % end
% % for k=1:length(Under_indeces)
% %     E(Under_indeces(k))=0;
% % end
% % for i=1:2:length(Over_indeces)-1
% %     Dump(i)=trapz(Time_Increment(Over_indeces(i):Over_indeces(i+1)),(A*Power_lim{10}(start_index+Over_indeces(i):start_index+Over_indeces(i+1))/60));
% % end
% % Dumped_Energy=sum(Dump);
% % for i=1:2:length(Under_indeces)-1
% %     Grid(i)=trapz(Time_Increment(Under_indeces(i):Under_indeces(i+1)),(Power_Consumption_Vector(Under_indeces(i):Under_indeces(i+1))/60));
% % end
% % Grid_Energy=sum(Grid);
% % for i=1:2:length(Middle_indeces)-1
% %     Converted(i)=trapz(Time_Increment(Middle_indeces(i):Middle_indeces(i+1)),(A*Power_lim{10}(start_index+Middle_indeces(i):start_index+Middle_indeces(i+1))/60));
% % end
% % Converted_Energy=sum(Converted);
% % 
% % subplot(2,1,2)
% % plot(Time_Increment_Days(start_index+1:end),E)
% % hold on
% % plot(Time_Increment_Days(start_index+1:end),E_over,'g')
% % plot(Time_Increment_Days(start_index+1:end),E_under,'r')
% % % plot(Time_Increment_Days(start_index+1:end),E_test,'k')
% % title('Battery Level')
% % xlabel('Time[Days]')
% % ylabel('Battery Level[Wh]')
% % legend('Battery Level','Dumped Energy','Grid Energy')

% -------------------------------------------------------------------------
% --------  Theoretical Energy Production and Consumption Plot 2  ---------
% -------------------------------------------------------------------------

start_index=5473;                       %Start at end of barge over ADCP bad data
Bin_Rep=5;                              %Bin{5} was chosen as a representative bin height it is approximately 12m off the seabed. It's also the highest bin that we always have reasonable data for.
Battery_Capacity=10000;                 %[Wh] 1 Tesla Power Wall
Battery_Start_Capacity=5000;            %[Wh] Assumed starting capacity
Energy_Consumption_Per_Month=900000;    %[Wh] Based off TECH 797 Report
Circuit_Voltage=120;                    %[V-AC] Assumed...this could be different in reality depending on circuit design. And could also be different coming from turbine vs coming off alternator to the load.
Power_Consumption=Energy_Consumption_Per_Month/(24*30); %[W] How many watts are consumed per hour assuming constant load
Power_Consumption_Vector=Power_Consumption*ones(1,length(Power_lim{Bin_Rep}(start_index:end))); %A vector of constant power consumption
Energy_Consumption_Vector=Power_Consumption_Vector/4;%The amount of [Wh] consumed in a vector where each element represents 15 minutes
Current_Load_Vector=Energy_Consumption_Vector/Circuit_Voltage;%This is the amount of current[A] that the load is drawing assuming everything is running on a 120 V circuit.
Power_Production_Vector=Power_lim{Bin_Rep}(start_index:end);
Energy_Production_Vector=Power_Production_Vector/4;%The amount of [Wh] produced in a vector where each element represents 15 minutes
Current_Source_Vector=Energy_Production_Vector/Circuit_Voltage;%This is the amount of current[A] that the turbine is producing assuming everything is running on a 120 V circuit.
%Battery Profile(This assumes linear profile between capacity[Ah] and charge level[V], I know this is not realistic)
N=100;%Number of profile points
Max_Voltage=450;%[V] Battery full voltage found on http://www.teslamotors.com/powerwall
Min_Voltage=350;%[V] Battery empty voltage found on http://www.teslamotors.com/powerwall
Charge_Prof=linspace(Min_Voltage,Max_Voltage,N);
Capacity_Prof=linspace(0,Battery_Capacity,N);
figure
subplot(2,1,1)
plot(Time_Increment_Days(start_index:end),Power_Production_Vector)
hold on
plot(Time_Increment_Days(start_index:end),Power_Consumption_Vector,'r')
title('Power Production and Consumption')
xlabel('Time[Days]')
ylabel('Power[W]')
legend('Power Production','Power Consumption')
j=1;
k=1;
l=1;
% Capacity_Level=zeros(length(Power_Production_Vector));
% Charge_Level
for i=2:length(Power_Production_Vector-1)
    Capacity_Level(1)=Battery_Start_Capacity;
    Charge_Level(1)=interp1(Capacity_Prof,Charge_Prof,Battery_Start_Capacity);
    Capacity_Level(i)=Capacity_Level(i-1)+Current_Source_Vector(i)-Current_Load_Vector(i);%subract out Ah consumed by load add in Ah produced by turbine
    Charge_Level(i)=interp1(Capacity_Prof,Charge_Prof,Capacity_Level(i));
    %Battery becomes full ATS=Open, Break=On
    if Charge_Level(i)>Max_Voltage
        Current_Source_Vector(i)=0;%Break is on thus no current comes from turbine to battery
        Power_Lost(i)=Power_Production_Vector(i);%The amount of power lost is equal to the amount of power that could have been generated
        j=j+1;
    end
    %Battery becomes empty ATS=Closed, Break=Off
    if Charge_Level(i)<Min_Voltage
        Current_Load_Vector(i)=0;%Current draw from the battery goes to zero because the load is drawing from the grid
        Power_Grid(i)=Power_Consumption_Vector(i);%All power consumption comes off the grid
        k=k+1;
    end
    %Battery is in between empty and full ATS=Open, Break=Off
    if Charge_Level(i)>Min_Voltage && Charge_Level(i)<Max_Voltage
        Power_Grid(i)=0;%No power is drawn from the grid
        Power_Lost(i)=0;%No power is wasted as the turbine spins and all possible power is stored in the battery and then consumed
        l=l+1;
    end

end
% %This might be ok or it might mess things up
% Under_indeces=zeros(length(Time_Increment_Days(start_index:end)));
% Over_indeces=zeros(length(Time_Increment_Days(start_index:end)));
% Middle_indeces=zeros(length(Time_Increment_Days(start_index:end)));
% %To here
% for i=2:length(Power_lim{10}(start_index:end))-1
%     Battery_Charge(1)=E_bat;
%     E_test(1)=E_bat;
%     E_test(i)=E_test(i-1)+(-Energy_Consumption_Vector(i)+(A*Power_lim{10}(start_index+i)/4));
%     E(i)=E(i-1)+(-Energy_Consumption_Vector(i)+(A*Power_lim{10}(start_index+i)/4));
%     if E(i)>E_bat;
%         E_over(i)=A*Power_lim{10}(start_index+i)/4;
%         E(i)=E(i-1)-Energy_Consumption_Vector(i);
%         Over_indeces(j)=i;
%         j=j+1;
%     end
%     if E(i)<0;
%         E_under(i)=-E(i);
%         E_over(i)=0;
%         Under_indeces(k)=i;
%         k=k+1;
%     end
%     if E(i)<E_bat && E(i)>0
%         E_under(i)=0;
%         Middle_indeces(l)=i;
%         l=l+1;
%         E_over(i)=0;
%     end
% end
% for k=1:length(Under_indeces)
%     E(Under_indeces(k))=0;
% end
% for i=1:2:length(Over_indeces)-1
%     Dump(i)=trapz(Time_Increment(Over_indeces(i):Over_indeces(i+1)),(A*Power_lim{10}(start_index+Over_indeces(i):start_index+Over_indeces(i+1))/60));
% end
% Dumped_Energy=sum(Dump);
% for i=1:2:length(Under_indeces)-1
%     Grid(i)=trapz(Time_Increment(Under_indeces(i):Under_indeces(i+1)),(Power_Consumption_Vector(Under_indeces(i):Under_indeces(i+1))/60));
% end
% Grid_Energy=sum(Grid);
% for i=1:2:length(Middle_indeces)-1
%     Converted(i)=trapz(Time_Increment(Middle_indeces(i):Middle_indeces(i+1)),(A*Power_lim{10}(start_index+Middle_indeces(i):start_index+Middle_indeces(i+1))/60));
% end
% Converted_Energy=sum(Converted);

subplot(2,1,2)
plot(Time_Increment_Days(start_index:end),Charge_Level)
hold on
% plot(Time_Increment_Days(start_index:end),E_over,'g')
% plot(Time_Increment_Days(start_index:end),E_under,'r')
% plot(Time_Increment_Days(start_index:end),E_test,'k')
title('Battery Level')
xlabel('Time[Days]')
ylabel('Battery Level[Wh]')
legend('Battery Level','Dumped Energy','Grid Energy')

% -------------------------------------------------------------------------
% --------       Calculate Energy Density in Top Bins         -----------
% -------------------------------------------------------------------------

% % for bins=5:50   %range over the top 1.75m where the turbine will be
% %     [ binE(bins) ] = onemonthbinE( Speed,Time_Increment,bins );
% % end
% % 
% % plot(1:50,binE)
% % % ylim([500000,600000])

% -------------------------------------------------------------------------
% --------  Theoretical Energy Production and Consumption Plot 2 -----------
% -------------------------------------------------------------------------
% A=[4.5,6.75,9];                                    %[m^2]
% Battery_Capacity=20000:20000:100000;                   %[Wh]
% A=7;                                    %[m^2]
% Battery_Capacity=100000;                   %[Wh]
% start_index=5473;
% Energy_Consumption_Per_Month=900000;    %[Wh]
% Power_Consumption=Energy_Consumption_Per_Month/(24*30); %[W]
% Power_Consumption_Vector=Power_Consumption*ones(1,length(Power_lim{10}(start_index:end))); %[W] check this
% Energy_Consumption_Vector=Power_Consumption_Vector/4; %[Wh]
% Power_Production_Vector=A*Power_lim{10}(start_index:end); %[W]
% Energy_Production_Vector=Power_Production_Vector/4; %[Wh]
% Total_Energy_Production=sum(Energy_Production_Vector);
% Total_Energy_Production_Check=trapz(Time_Increment(start_index:end),Power_Production_Vector)/60;
% Total_Energy_Consumption=sum(Energy_Consumption_Vector);
% Total_Energy_Consumption_Check=trapz(Time_Increment(start_index:end),Power_Consumption_Vector)/60;
% figure
% subplot(2,1,1)
% plot(Time_Increment_Days(start_index:end),Power_Production_Vector)
% hold on
% plot(Time_Increment_Days(start_index:end),Power_Consumption_Vector,'r')
% title('Power Production and Consumption')
% xlabel('Time[Days]')
% ylabel('Power[W]')
% legend('Power Production','Power Consumption')
% j=1;
% k=1;
% l=1;
% for i=2:length(Power_lim{10}(start_index:end))-1
%     Battery_Level(1)=Battery_Capacity/2;%Start the battery at half a charge
%     %The battery level at each index is equal to whatever the previous
%     %indexes state of the battery was minus the amount of energy consumed
%     %and plus the amount of energy converted
%     Battery_Level(i)=Battery_Level(i-1)+(-Energy_Consumption_Vector(i)+Energy_Production_Vector(i));
%     %If the battery is neither empty or over full no energy is drawn from
%     %the grid or dumped.
%     if Battery_Level(i)<Battery_Capacity && Battery_Level(i-1)>0
%         E_under(i)=0;
%         E_over(i)=0;
%         Middle_indeces(l)=i;
%         l=l+1;
%     end
%     %If the battery level goes above the batteries full capacity the amount
%     %of energy overproduced will be whatever the battery level is minus the
%     %full capacity of the battery such that the battery's charge level stays
%     %at its full capacity
%     if Battery_Level(i)>=Battery_Capacity;
%         E_over(i)=Battery_Level(i)-Battery_Capacity;
%         Battery_Level(i)=Battery_Capacity;
%         E_under(i)=0;
%         Over_indeces(j)=i;
%         j=j+1;
%     end
%     %If the battery gets empty the automatic switch switches the bridge to
%     %grid power, whatever the energy demand is is pulled from the grid, and
%     %whatever the turbine converts goes into charging the battery.
%     if Battery_Level(i)<=0;
%         E_under(i)=Energy_Consumption_Vector(i);
%         Battery_Level(i)=Energy_Production_Vector(i);
%         E_over(i)=0;
%         Under_indeces(k)=i;
%         k=k+1;
%     end
% end
% 
% for i=1:2:length(Over_indeces)-1%index every other such that each space is only integrated once
%     Dump(i)=trapz(Time_Increment(Over_indeces(i):Over_indeces(i+1)),...
%         (A*Power_lim{10}(start_index+Over_indeces(i):start_index+Over_indeces(i+1))/60));%/60 because Time_Increments is in minutes
% end
% Dumped_Energy=sum(Dump);
% for i=1:2:length(Under_indeces)-1
%     Grid(i)=trapz(Time_Increment(Under_indeces(i):Under_indeces(i+1)),...
%         (Power_Consumption_Vector(Under_indeces(i):Under_indeces(i+1))/60));
% end
% Grid_Energy=sum(Grid);
% for i=1:2:length(Middle_indeces)-1
%     Converted(i)=trapz(Time_Increment(Middle_indeces(i):Middle_indeces(i+1)),...
%         (A*Power_lim{10}(start_index+Middle_indeces(i):start_index+Middle_indeces(i+1))/60));
% end
% Converted_Energy=sum(Converted);
% 
% subplot(2,1,2)
% plot(Time_Increment_Days(start_index+1:end),Battery_Level)
% hold on
% % plot(Time_Increment_Days(start_index+1:end),E_over,'g')
% % plot(Time_Increment_Days(start_index+1:end),E_under,'r')
% title('Battery Level')
% xlabel('Time[Days]')
% ylabel('Battery Level[Wh]')
% legend('Battery Level','Dumped Energy','Grid Energy')