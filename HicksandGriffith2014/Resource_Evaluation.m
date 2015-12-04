clear all; close all; clc
A=5.25;              %Turbine Area[m^2](1 for power density)
S=29.8;             %Average Salinity[PPT] over a year(2004) at Great Bay Coastal Marine Lab (greatbaydata.org)
T=8.749;            %Average Temperature[C] over a year(2004) at Great Bay Coastal Marine Lab (greatbaydata.org)
sample_period=15;   %minutes(check this)
cut_in=0.7;         %[m/s]
Cp=0.35;            %Efficiency of turbine
bin_space=0.3;      %[m]
start_index=5473;   %Start at end of barge over ADCP bad data

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
    Direction_deg{i}=Data.adcp(:,i,14)+16.25;%16.25 is magnetic variation given in MEMBRVEL.RPT file
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
% %     title(strcat('Current Speeds at:',num2str(Bin_Height{1,i}(1,1)),'m Above Seafloor'))
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

% % for i=5:c(2)
% %     l = 1;
% %     for j = 1:length(Velocity{1,i}(:,1))-1
% %         if Velocity{1,i}(j,1)>0
% %             Ebb_Speed{1,i}(j) = abs(Velocity{1,i}(j,1));
% %             Flood_Speed{1,i}(j) = 0;
% %         end
% %         if Velocity{1,i}(j,1)<0
% %             Flood_Speed{1,i}(j) = abs(Velocity{1,i}(j,1));
% %             Ebb_Speed{1,i}(j) = 0;
% %         end
% %         if isnan(Velocity{1,i}(j,1)) == 1 || Velocity{1,i}(j,1) == 0
% %             Flood_Speed{1,i}(j) = 0;
% %             Ebb_Speed{1,i}(j) = 0;
% %         end
% %         Ebb_Power_D{1,i}(j)=(1/2)*rho(i)*Ebb_Speed{1,i}(j).^3;
% %         Flood_Power_D{1,i}(j)=(1/2)*rho(i)*Flood_Speed{1,i}(j).^3;
% %         if Velocity{1,i}(j,1)*Velocity{1,i}(j+1,1)<0
% %             cycle_indices{i}(l)=j+1; 
% %             l=l+1;
% %         end
% %     end
% %     for k = 1:2:length(cycle_indices{i})-1
% %         Ebb_Cycle_Energy{1,i}(k)=trapz(Time_Increment(cycle_indices{i}(k):cycle_indices{i}(k+1)),Ebb_Power_D{1,i}(cycle_indices{i}(k):cycle_indices{i}(k+1)))/4; %[Wh]
% %     end
% %     for m= 1:length(Ebb_Cycle_Energy{1,i})-1
% %         if Ebb_Cycle_Energy{1,i}(m) == 0
% %             Ebb_Cycle_Energy{1,i}(m)=Ebb_Cycle_Energy{1,i}(m+1);
% %             Ebb_Cycle_Energy{1,i}(m+1)=0;
% %         end
% %     end
% %     for k=2:2:length(cycle_indices{i})-1
% %         Flood_Cycle_Energy{1,i}(k)=trapz(Time_Increment(cycle_indices{i}(k):cycle_indices{i}(k+1)-1),Flood_Power_D{1,i}(cycle_indices{i}(k):cycle_indices{i}(k+1)-1))/4; %[Wh]      
% %     end
% % end
% % 
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
% %     title(strcat('Depth:',num2str(mean(ADCP_Depth))))
% %     subplot(3,1,1);
% %     bin_centers=[1000:1000:16000];
% %     hist(energy{i},bin_centers);
% %     title(strcat('Energy per Half Tidal Cycle 1m^2 Turbine at:',num2str(Bin_Height{1,i}(1,1)),'m Above Seafloor'))
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
% %     hist(Ebb_Cycle_Energy{i},bin_centers);
% %     title(strcat('Energy per Ebb Tidal Cycle 1m^2 Turbine at:',num2str(Bin_Height{1,i}(1,1)),'m Above Seafloor'))
% %     ylabel('Frequency')
% %     xlabel('Energy per Cycle[Wh]')
% %     xlim([xmin,xmax])
% %     ylim([ymin,ymax])
% %     
% %     subplot(3,1,3);
% %     hist(Flood_Cycle_Energy{i},bin_centers);
% %     title(strcat('Energy per Flood Tidal Cycle 1m^2 Turbine at:',num2str(Bin_Height{1,i}(1,1)),'m Above Seafloor'))
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
% %     polar(Direction{i}(start_index:end),Speed{i}(start_index:end),'.b')
% %     set(gca,'color','none')
% %     view([90 -90])
% %     title({strcat('Current Speed[m/s] and Direction at Approximately',num2str(Bin_Height{1,i}(1,1)),'m Above Seafloor'),''})
% %     pause
% %     close all
% % end

% -------------------------------------------------------------------------
% --Plot Comparing Energy Yeild vs Turbine Size, Efficieny, and Cut-In ----
% -------------------------------------------------------------------------
% % 
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

% % Diameter=3;         %[m]
% % Cp=0.40;
% % cut_in=1;           %[m/s]
% % for Turbine_Length=1:5
% %     [ Thirty_Day_Energy(Turbine_Length) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % figure(1)
% % hold on
% % plot((Diameter*(1:5)),Thirty_Day_Energy,'--m')
% % cut_in=0.7;
% % for Turbine_Length=1:5
% %     [ Thirty_Day_Energy(Turbine_Length) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(1:5)),Thirty_Day_Energy,':m')
% % cut_in=0.5;
% % for Turbine_Length=1:5
% %     [ Thirty_Day_Energy(Turbine_Length) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(1:5)),Thirty_Day_Energy,'-m')
% % Cp=0.35;
% % cut_in=1;
% % for Turbine_Length=1:5
% %     [ Thirty_Day_Energy(Turbine_Length) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(1:5)),Thirty_Day_Energy,'--g')
% % cut_in=0.7;
% % for Turbine_Length=1:5
% %     [ Thirty_Day_Energy(Turbine_Length) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(1:5)),Thirty_Day_Energy,':g')
% % cut_in=0.5;
% % for Turbine_Length=1:5
% %     [ Thirty_Day_Energy(Turbine_Length) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(1:5)),Thirty_Day_Energy,'-g')
% % Cp=0.25;
% % cut_in=1;
% % for Turbine_Length=1:5
% %     [ Thirty_Day_Energy(Turbine_Length) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(1:5)),Thirty_Day_Energy,'--k')
% % cut_in=0.7;
% % for Turbine_Length=1:5
% %     [ Thirty_Day_Energy(Turbine_Length) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(1:5)),Thirty_Day_Energy,':k')
% % cut_in=0.5;
% % for Turbine_Length=1:5
% %     [ Thirty_Day_Energy(Turbine_Length) ] = onemonthE( Speed,bin_space,Time_Increment,Turbine_Length,Diameter,Cp,cut_in );
% % end
% % plot((Diameter*(1:5)),Thirty_Day_Energy,'-k')
% % 
% % Thirty_Day_Energy = onemonthE( Speed,bin_space,Time_Increment,1.75,3,.35,.5 );
% % plot((1.75*3),Thirty_Day_Energy,'b*')
% % 
% % xlabel('Turbine Swept Area[m^2]')
% % ylabel('Energy Converted in 30 Days[Wh]')
% % legend('Cp=0.40, Cut In=1m/s','Cp=0.40, Cut In=0.7m/s','Cp=0.40, Cut In=0.5m/s','Cp=0.35, Cut In=1m/s','Cp=0.35, Cut In=0.7m/s','Cp=0.35, Cut In=0.5m/s','Cp=0.25, Cut In=1m/s','Cp=0.25, Cut In=0.7m/s','Cp=0.25, Cut In=0.5m/s','Location','NorthWest')


% -------------------------------------------------------------------------
% --------  Theoretical Energy Production and Consumption Plot   ----------
% -------------------------------------------------------------------------

% % Bin_Rep=6;                                  %Bin{6} was chosen as a representative bin height it is approximately 12m off the seabed. It's also the highest bin that we always have reasonable data for.
% % Battery_Capacity=10000;                    %[Wh] 1 Tesla Power Wall
% % Battery_Start_Capacity=Battery_Capacity/2;  %[Wh] Assumed starting capacity
% % Energy_Consumption_Per_Month=900000;        %[Wh] Based off TECH 797 Report
% % Circuit_Voltage=120;                        %[V-AC] Assumed...this could be different in reality depending on circuit design. And could also be different coming from turbine vs coming off alternator to the load.
% % Power_Consumption=Energy_Consumption_Per_Month/(24*30); %[W] How many watts are consumed per hour assuming constant load
% % Power_Consumption_Vector=Power_Consumption*ones(1,length(Power_lim{Bin_Rep}(start_index:end))); %A vector of constant power consumption
% % Energy_Consumption_Vector=Power_Consumption_Vector/4;%The amount of [Wh] consumed in a vector where each element represents 15 minutes
% % Power_Production_Vector=Power_lim{Bin_Rep}(start_index:end);
% % Energy_Production_Vector=Power_Production_Vector/4;%The amount of [Wh] produced in a vector where each element represents 15 minutes
% % [ Energy_Level,Power_Lost,Power_Grid,Energy_Lost,Energy_Grid,Energy_Conv ] = Bat_Sim( Power_Production_Vector,Power_Consumption_Vector,Battery_Capacity );
% % 
% % figure
% % subplot(2,1,1)
% % plot(Time_Increment_Days(start_index:end),Power_Production_Vector,'b')
% % hold on
% % plot(Time_Increment_Days(start_index:end),Power_Consumption_Vector,'r')
% % plot(Time_Increment_Days(start_index:end),Power_Lost,'mo')
% % plot(Time_Increment_Days(start_index:end),Power_Grid,'k','LineWidth',3)
% % title('Power Production and Consumption')
% % xlabel('Time[Days]')
% % ylabel('Power[W]')
% % legend('Power Production','Power Consumption','Over Production','Under Production')
% % subplot(2,1,2)
% % plot(Time_Increment_Days(start_index:end-1),Energy_Level)
% % title('Battery Level')
% % xlabel('Time[Days]')
% % ylabel('Battery Level[Wh]')
% % axes('Position',[0 0 1 1],'Visible','off') % sets ax1 to current axes
% % str = {strcat('Turbine Area:',num2str(A),'m^2',', Battery Capacity:',num2str(Battery_Capacity),'Wh',', Over Producing:',num2str(100*(Energy_Lost/(Energy_Grid+Energy_Conv+Energy_Lost))), '%, Under Producing:',num2str(100*(Energy_Grid/(Energy_Grid+Energy_Conv+Energy_Lost))),'%')};
% % text(0,0.02,str(:))

% -------------------------------------------------------------------------
% ------------------  Depth vs Speed vs Time  -----------------------------
% -------------------------------------------------------------------------

% % for i = 1:length(Bin_Height)
% %     Bin_Heights(i) = Bin_Height{1,i}(1,1);
% % end
% % for i = 1:c(2)
% %     for j = start_index:c(1)
% %         Speeds(i,(j-start_index)+1)=Speed{1,i}(j,1);
% %     end
% % end
% % 
% % [time,depth]=meshgrid(Time_Increment_Days(start_index:end),Bin_Heights);
% % 
% % hFig = figure(1);
% % set(hFig, 'Position', [100 100 1000 200 ])
% % pcolor(time,depth,Speeds)
% % shading interp
% % hold on
% % plot(Time_Increment_Days(start_index:end),ADCP_Depth(start_index:end))
% % color = colorbar;
% % color.Label.String = 'Speed[m/s]';
% % ylim([0,20])
% % xlabel('Time[days]')
% % ylabel('Depth[m]')
% % 
% % for i = 1:length(Bin_Height)
% %     Bin_Heights(i) = Bin_Height{1,i}(1,1);
% % end
% % clear Speeds
% % for i = 1:c(2)
% %     for j = 6240:7200
% %         Speeds(i,(j-6240)+1)=Speed{1,i}(j,1);
% %     end
% % end
% % 
% % [time,depth]=meshgrid(Time_Increment_Days(6240:7200),Bin_Heights);
% % 
% % hFig = figure(2);
% % set(hFig, 'Position', [100 100 1000 200 ])
% % pcolor(time,depth,Speeds)
% % shading interp
% % hold on
% % plot(Time_Increment_Days(6240:7200),ADCP_Depth(6240:7200))
% % color = colorbar;
% % color.Label.String = 'Speed[m/s]';
% % ylim([0,20])
% % xlabel('Time[days]')
% % ylabel('Depth[m]')

% -------------------------------------------------------------------------
% --------       Calculate Energy Density in Top Bins         -------------
% -------------------------------------------------------------------------

% % for bins=5:50   %range over the top 1.75m where the turbine will be
% %     [ binE(bins) ] = onemonthbinE( Speed,Time_Increment,bins );
% % end
% % 
% % plot(1:50,binE)
% % ylim([500000,600000])

% -------------------------------------------------------------------------
% --------       Calculate MKPD, MKPA, Maximum Speed         --------------
% -------------------------------------------------------------------------

% % Bin_Rep=6;                                  %Bin{6} was chosen as a representative bin height it is approximately 12m off the seabed. It's also the highest bin that we always have reasonable data for.
% % MKPD = nanmean(0.5*rho(Bin_Rep)*Speed{1,Bin_Rep}(start_index:end,1).^3) %[w/m^2]
% % for i = 1:length(Velocity{1,Bin_Rep}(:,1))
% %     if Velocity{1,Bin_Rep}(i,1)>0
% %         Ebb_Speed(i) = abs(Velocity{1,Bin_Rep}(i,1));
% %         Flood_Speed(i) = NaN;
% %     end
% %     if Velocity{1,Bin_Rep}(i,1)<0
% %         Flood_Speed(i) = abs(Velocity{1,Bin_Rep}(i,1));
% %         Ebb_Speed(i) = NaN;
% %     end
% % end
% %         
% % MKPD_Ebb = nanmean(0.5*rho(Bin_Rep)*Ebb_Speed(1,start_index:end).^3);
% % MKPD_Flood = nanmean(0.5*rho(Bin_Rep)*Flood_Speed(1,start_index:end).^3);
% % MKPDA = MKPD_Ebb/MKPD_Flood
% % Max_Speed = max(abs(Velocity{1,Bin_Rep}(start_index:end,1)))

% -------------------------------------------------------------------------
% --------               Direction Analysis                  --------------
% -------------------------------------------------------------------------

% % for i = 5:c(2)
% %     k=1;
% %     l=1;
% %     for j = 1:length(Direction_deg{i})
% %         if Direction_deg{i}(j) > 180
% %             Flood_Direction{i}(k) = Direction_deg{i}(j);
% %             k=k+1;
% %         end
% %         if Direction_deg{i}(j)<180
% %             Ebb_Direction{i}(l) = Direction_deg{i}(j);
% %             l=l+1;
% %         end
% %     end
% %     avg_flood_dir_el(i)=mean(Flood_Direction{i});
% % 
% %     avg_ebb_dir_el(i)=mean(Ebb_Direction{i});
% % 
% % end
% % avg_ebb_dir=mean(avg_ebb_dir_el(5:end));
% % avg_flood_dir=mean(avg_flood_dir_el(5:end));
% % 
% % Direction_Asymmetry=(avg_flood_dir-180)-avg_ebb_dir;


% -------------------------------------------------------------------------
% --------                Pier AoA Analysis                  --------------
% -------------------------------------------------------------------------

%Bridge pier is aligned 105.05 degrees south of north. Determined using
%bridge as built drawings.

Bridge_Angle=105.05;

for i = 5:c(2)
    for j = 1:length(Direction_deg{i})
        if Speed{i}(j)<0.5
            Direction_deg{i}(j)=NaN;
        end
        if Direction_deg{i}(j) > 180
            AoA{i}(j) = (Direction_deg{i}(j)-180)-Bridge_Angle;
        end
        if Direction_deg{i}(j)<180
            AoA{i}(j) = Direction_deg{i}(j)-Bridge_Angle;
        end
    end
    plot(Time_Increment_Days(start_index:length(AoA{i})),AoA{i}(start_index:end))%Time_Increment_Days is incremented weird because AoA ends in a few NaNs
    title({strcat('Bridge Pier AoA at Approximately',num2str(Bin_Height{1,i}(1,1)),'m Above Seafloor'),''})
    xlabel('Time[days]')
    ylabel('AoA[deg]')
    pause
    close all
end

%For thin airfoils and small angles of attack Cl=2*pi*AoA (AoA in radians) https://www.grc.nasa.gov/www/k-12/airplane/incline.html
Planform_A=(80.74*75)*0.092903;%[m^2]underwater bridge pier planform area dimensions taken from bridge drawings converted to m^2
Cl=2*pi*deg2rad(AoA{6});
Lift=(1/2)*mean(rho)*Planform_A*Cl'.*Speed{6}(1:length(Cl)).^2;
plot(Time_Increment_Days(start_index:length(Lift)),Lift(start_index:end))%Time_Increment_Days is incremented weird because AoA ends in a few NaNs
title({strcat('Bridge Pier Lift Force'),''})
xlabel('Time[days]')
ylabel('Lift Force[N]')