clear; close all; clc;

% *************************************************************************
% Developed by: Hayden Hicks 01/19/2014
% Revisions: (Person) - (Short Description)
%
% Description: This .m file is designed to input ADCP data. All bins are
% automatically imported and plotted. In addition, a theoretical model is
% applied to all bins. This model can be used to describe the tidal current
% velocities over any period of time. The theoretical model can be used to
% represent velocities at any location assuming that the max velocities is
% known.
%
% *************************************************************************

% Important for acquiring data
number_of_bins = 16;                % # bins of acoustic signal to surf.
Assessment_name = 'PIR0705';        % Common name of all bins
ping_separation = 6;                % Time difference between pings (minutes)
Operating_Range = 1;                % Minimum required velocity for power
Turbine_Efficiency = 0.35;          % Approximate Gorlov Turbine Efficiency
A = 12;                             % In m^2 (when A = 1, then Pwr Density)
bin_size = 1;                       % In meters
i = 1;                              % Index for looping through bins

% Water density
ro = 1025; % kg/m^3

% Creates matching file names for loading data
for bin = 1:number_of_bins
    
    if bin < 10   % Formats file name if bin # is less than 10
        number = num2str(bin);
        bin_name = strcat('0',number);
        file_name{i} = strcat(Assessment_name,'_bin', bin_name, '.dat');
    elseif bin >= 10   % Formats file name if bin # is greater than 10
        number = num2str(bin);
        bin_name = strcat(number);
        file_name{i} = strcat(Assessment_name,'_bin', bin_name, '.dat');
    end
    
    i = i + 1;
end

% Reset index for importing and accesing data
i = 1;
bin = 1;

% Imports Header (Roll, Pitch, Temp., Depth)
Header = importdata(strcat(Assessment_name, '.hdr'),' ');

ADCP_Heading = Header.data(:,2);
ADCP_Pitch = Header.data(:,3);
ADCP_Roll = Header.data(:,4);
ADCP_Temp = Header.data(:,5);
ADCP_Depth = Header.data(:,6);


% Imports data from newly constructed file names
for bin = 1:number_of_bins
    [Bin_data{i},delimit_char,headerlinesOut] = importdata(file_name{i},' ');
    Julian_Time{i} = Bin_data{i}.data(:,1);
    Speed{i} = Bin_data{i}.data(:,2)/100;      % /100 to convert to m/s
    Direction{i} = Bin_data{i}.data(:,3);
    Vel_North{i} = Bin_data{i}.data(:,4)/100;  % /100 to convert to m/s
    Vel_East{i} = Bin_data{i}.data(:,5)/100;   % /100 to convert to m/s
    Vel_Vert{i} = Bin_data{i}.data(:,6)/100;   % /100 to convert to m/s
    Date{i} = Bin_data{i}.textdata(headerlinesOut+1:end,1);
    Time{i} = Bin_data{i}.textdata(headerlinesOut+1:end,2);
    Depth{i} = i*bin_size;
    Power{i} = (1/2)*A*ro*Speed{i}.^3;
    Turbine_Power{i} = Turbine_Efficiency.*(1/2).*A*ro.*Speed{i}.^3;
    
    % Determine the length of the data file
    Data_Length = length(Speed{i});
    
    % This loop searches for velocities that are of a magnitude that can
    % produce power. Comment when no interest (slow process).
    % ------------------------------------------------------------------
    for count = 1:Data_Length
        
        if bin == 1 
        % Create an array of the time increments
        Time_Increment{1}.data(count,1) = ping_separation*count;
        end
        
        if Bin_data{i}.data(count,2) < Operating_Range*100
            Vel_lim{i}.data(count,1) = 0;
            Power_lim{i}.data(count,1) = 0;
            Turbine_Power_lim{i}.data(count,1) = 0;
            
        else
            Vel_lim{i}.data(count,1) = Bin_data{i}.data(count,2)/100;
            Power_lim{i}.data(count,1) = (1/2).*A.*ro.*(Vel_lim{i}.data(count,1)).^3;
            Turbine_Power_lim{i}.data(count,1) = Turbine_Efficiency.*Power_lim{i}.data(count,1);
        end
    end
    
    Time_Increment{i}.data(:,1) = Time_Increment{1}.data(:,1);
    
    % Convert date and time from ADCP to usable form for plots
    date_time = [Date{i} Time{i}];
    date_time_new = cell2mat(date_time);
    final_time{i} = datenum(date_time_new,'yyyy-mm-ddHH:MM:SS');

    i = i + 1;
end

% Determines percent difference from bin with "max" vel.
i = 1;
bin = 1;

for bin = 1:number_of_bins
    Vel_Difference{i} = ((abs(Speed{number_of_bins} - Speed{i}))); %./Speed{number_of_bins}).*100;
    Power_Difference{i} = ((abs(Power{number_of_bins} - Power{i})));%./Power_Density{number_of_bins}).*100;
    
    Percent_Vel_Diff{i} = sum(Vel_Difference{i}./(sum(Speed{i})./Data_Length)*100)./Data_Length;
    Percent_Vel_Diff{i}';
    
    Percent_Power_Diff{i} = sum(Power_Difference{i}./(sum(Power{i})./Data_Length)*100)./Data_Length;
    Percent_Power_Diff{i}';
    
    i = i + 1;
end

% Used for plotting
Length_of_Time = length(final_time{i-1});

i = 1;

% Color array for plots
for xxx = .8:-.2:.2
    for yyy = .8:-.2:0
        for zzz = 1:-.2:0
            c1(i) = {[xxx yyy zzz]};
            i = i + 1;
        end
    end
end

% Determine peaks and troughs of power generation to descretize tidal cycles

for i=1:number_of_bins
    [Power_Smoothed{i},L] = wsmooth(Power{i},final_time{i},2);
    [maxtab{i}, mintab{i}]=peakdet(Power_Smoothed{i}, 1000, final_time{i});
    hold on
end

figure(1)
plot(final_time{1},Power_Smoothed{1})
plot(maxtab{1}(:,1), maxtab{1}(:,2), 'r*');
plot(mintab{1}(:,1), mintab{1}(:,2), 'g*');

figure(2)
plot(final_time{1},Speed{1})

for i=1:number_of_bins
    k=1;
    for j=1:length(final_time{i})
        if final_time{i}(j)==mintab{i}(k,1)
            index{i}(k)=j;
            k=k+1;
            if k>length(mintab{i})
                break
            end
        end
    end
    for l=1:length(index{i})-1
        energy{i}(l)=trapz(final_time{i}(index{i}(l):index{i}(l+1)),Power{i}(index{i}(l):index{i}(l+1)));
    end
    ebb_energy{i}=energy{i}(1:2:end);
    flood_energy{i}=energy{i}(2:2:end);
end



% delta=200;
% for i=1:number_of_bins
%     [max_Power{i}, min_Power{i}]=peakdet(Power{i}, delta, final_time{i});
%     k=1;
%     for j=1:length(final_time{i})
%         if final_time{i}(j)==min_Power{i}(k,1)
%             index{i}(k)=j;
%             k=k+1;
%             if j==length(final_time{i})-10
%                 pause(120)
%             end
%         end
%     end
% end



% Reset index for generating plots
i = 1;
bin = 1;



% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% ::::::::::::::::::: Code Used for Generation of Plots :::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

% -------------------------------------------------------------------------
% --- Plot of Magnitude of Velocity [cm/s] and Temperature [F] (scroll) ---
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % figure
% % for bin = 1:number_of_bins
% %     
% %     hold on
% %     plot(final_time{i},Speed{i}.*100,'g')
% %     plot(final_time{i},ADCP_Temp(:,1).*1.8 + 32)
% % 
% %     i = i + 1;
% % end
% % NumTicks = 200;
% % L = get(gca,'XLim');
% % set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% % datetick('x','mmm-dd HH:MM PM','keepticks')
% % title('Magnitude of Velocity and Temperature vs. Time')
% % ylabel('Magnitude of Velocity [cm/s] and Temperature F')
% % legend('Magnitude of Velocity [cm/s]','Temperature F')
% % axis([733264.711805556 733289.370138889 0 200])
% % addScrollbar(gca, 1);
% % hold off


% -------------------------------------------------------------------------
% ------ Plots of Available Power, Operating Power, and Turbine Power -----
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % for bin = 1:number_of_bins
% %     
% %     figure
% %     hold on
% %     plot(final_time{i},Power{i})
% %     plot(final_time{i},Power_lim{i}.data(:,1),'g')
% %     plot(final_time{i},Turbine_Power{i},'m')
% %     plot(final_time{i},Turbine_Power_lim{i}.data(:,1),'k')
% %     datetick('x','mmm-dd','keepticks')
% %     legend('Total Available Power','Available Power for Operating Range','Full Range Turbine','35 Percent Turbine')
% %     hold off
% %     
% %     i = i + 1;
% % end

% -------------------------------------------------------------------------
% --- Plots of Available PWR, Oper. PWR, Turbine PWR w/AREA (and scroll) --
% -------------------------------------------------------------------------

i = 1;
bin = 1;

for bin = 1:number_of_bins

    PD{i} = trapz(Time_Increment{i}.data(:,1),Power{i})./60;
    PL{i} = trapz(Time_Increment{i}.data(:,1),Power_lim{i}.data(:,1))./60;
    PAT{i} = trapz(Time_Increment{i}.data(:,1),Turbine_Power{i})/60;
    PHT{i} = trapz(Time_Increment{i}.data(:,1),Turbine_Power_lim{i}.data(:,1))/60;
    
    PD{i} = cumtrapz(Time_Increment{i}.data(:,1),Power{i})./60;
    PL{i} = cumtrapz(Time_Increment{i}.data(:,1),Power_lim{i}.data(:,1))./60;
    PAT{i} = cumtrapz(Time_Increment{i}.data(:,1),Turbine_Power{i})/60;
    PHT{i} = cumtrapz(Time_Increment{i}.data(:,1),Turbine_Power_lim{i}.data(:,1))/60;
    
    Diff_PDPL{i} = PD{i} - PL{i};
    Diff_PATPHT{i} = PAT{i} - PHT{i};
    
    Diff_PDPL{i}';
    Diff_PATPHT{i}';
    
    figure
    
    subplot(2,1,1);
    hold on
    H1 = area(final_time{i},Power{i});
    H2 = area(final_time{i},Power_lim{i}.data(:,1));
    H3 = area(final_time{i},Turbine_Power{i});
    H4 = area(final_time{i},Turbine_Power_lim{i}.data(:,1));
    set(H2,'FaceColor',[0,0.25,0.25]);
    set(H3,'FaceColor',[0,0.5,0.5]);
    set(H4,'FaceColor',[0,0.75,0.75]);
    datetick('x','mmm-dd','keepticks')
    NumTicks = 300;
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    datetick('x','mmm-dd HH:MM PM','keepticks')
    set(gca,'XMinorTick','on','YMinorTick','on') 
    addScrollbar(gca, 1);
    ylabel('Power [Watts]')
    title('Power')
    grid on
    legend('Total Available Power','Available Power for Operating Range','Full Range Turbine','35 Percent Turbine')
    hold off
    
    subplot(2,1,2);
    hold on
    plot(final_time{i},Speed{i})
    plot(final_time{i},1)
    NumTicks = 300;
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    datetick('x','mmm-dd HH:MM PM','keepticks')
    set(gca,'XMinorTick','on','YMinorTick','on') 
    addScrollbar(gca, 1);
    ylabel('Magnitude of Velocity [m/s]')
    hold off
    
    
    i = i + 1;
end

% -------------------------------------------------------------------------
% ---- Colormap: Depth vs. Time where Velocity is represented by color ----
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % figure
% % for bin = 1:number_of_bins
% %     
% %     Vel_East{i} = Vel_East{i}';
% %     Vel_East{i} = flipud(Vel_East{i});
% %     
% %     % Plot of depth versus time where velocity is shown as color map
% %     hold on
% %     imagesc(final_time{i},Depth{i},Vel_East{i})
% %     set(gca,'Ydir','Normal')
% %     datetick('x','mmm-dd','keepticks')
% %     colorbar
% %     caxis([1 2])
% %     title('Velocity and Depth versus Time','Fontsize',14)
% %     xlabel('Time')
% %     ylabel('Depth')
% %     
% %     i = i + 1;
% % end
% % hold off

% -------------------------------------------------------------------------
% --- Colormap: Depth vs. Time where Mag. of Vel. is shown as color map ---
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % figure
% % for bin = 1:number_of_bins
% %     
% %     Speed{i} = Speed{i}';
% %     Speed{i} = flipud(Speed{i});
% %     
% %     % Plot of depth versus time where velocity is shown as color map
% %     hold on
% %     imagesc(final_time{i},Depth{i},Speed{i})
% %     set(gca,'Ydir','Normal')
% %     datetick('x','mmm-dd','keepticks')
% %     colorbar
% %     caxis([0 2])
% %     title('Velocity and Depth versus Time','Fontsize',14)
% %     xlabel('Time')
% %     ylabel('Depth')
% %     
% %     i = i + 1;
% % end
% % hold off

% -------------------------------------------------------------------------
% ------------ Contour Plot where Speed is represented by color -----------
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % % Construct velocity matrix
% % for bin = 1:number_of_bins
% %     Z(i,:) = Speed{i}';
% %     
% %     i = i + 1;
% % end
% % 
% % i = 1;
% % bin = 1;
% % 
% % figure
% % contour3(final_time{i},1:number_of_bins,Z,30)
% % datetick('x','mmm-dd','keepticks')
% % xlabel('Time')
% % ylabel('Depth [m]')
% % zlabel('Velocity [m/s]')

% -------------------------------------------------------------------------
% ------- Surface Plot of Depth vs. Time vs. Magnitude of Velocity --------
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % % Construct velocity matrix
% % for bin = 1:number_of_bins
% %     Z(i,:) = Speed{i};
% %     
% %     i = i + 1;
% % end
% % 
% % i = 1; % Dont remove this index, it needs to be reset here
% % figure
% % surfc(final_time{i},1:number_of_bins,Z)
% % datetick('x','mmm-dd','keepticks')
% % ylabel('Depth [m]')
% % zlabel('Velocity [m/s]')
% % axis tight

% -------------------------------------------------------------------------
% ---------------------- Plot of Velocity vs. Depth -----------------------
% Add scroll
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % figure
% % set(gcf,'Color',[0,.5,.5])
% % for bin = 1:number_of_bins
% %     hold on
% %     plot(final_time{i},Speed{i},'Color',c1{i})
% %     datetick('x','mmm-dd','keepticks')
% %     legendInfo{i} = ['Bin #: ' num2str(i)]; 
% %     i = i + 1;
% % end
% % legend(legendInfo,'Location','EastOutside')
% % axis tight
% % hold off

% -------------------------------------------------------------------------
% ----------------------- Plot of Power vs. Depth -------------------------
% Add scroll
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % figure
% % set(gcf,'Color',[0,.5,.5])
% % for bin = 1:number_of_bins
% %     hold on
% %     plot(final_time{i},Power{i},'Color',c1{i})
% %     datetick('x','mmm-dd','keepticks')
% %     legendInfo{i} = ['Power Bin #: ' num2str(i)]; 
% %     i = i + 1;
% % end
% % legend(legendInfo,'Location','EastOutside')
% % hold off

% -------------------------------------------------------------------------
% ----------------- Plot of # of Occurrences of Velocity ------------------
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % for bin = 1:number_of_bins
% %     
% %     % Plots the ADCP Experimental Velocity Measurements
% %     figure
% %     hist_bins = 0:1:2;
% %     hist(Vel_lim{i}.data(:,1),hist_bins)
% %     title('Number of Occurrences of Velocity','Fontsize',12);
% %     xlabel('Velocity [m/s]')
% %     ylabel('# of Occurrences')
% %     grid on;
% %     
% %     i = i + 1;
% % end

% -------------------------------------------------------------------------
% -------------- Plot of both Power and Magnitude of Velocity -------------
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % for bin = 1:number_of_bins
% %     
% %     %Plots the ADCP Experimental Velocity Measurements
% %     figure
% %     hold on
% %     subplot(2,1,1);
% %     plot(final_time{i},Power{i});
% %     datetick('x','mmm-dd','keepticks')
% %     ylabel('Power [Watts]')
% %     axis tight
% %     
% %     subplot(2,1,2);
% %     plot(final_time{i},Speed{i});
% %     datetick('x','mmm-dd','keepticks')
% %     ylabel('Magnitude of Velocity [m/s]')
% %     axis tight
% %     hold off
% %     
% %     i = i + 1;
% % end

% -------------------------------------------------------------------------
% ----------- Plot of Magnitude of Velocity [m/s] with scrollbar ----------
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % for bin = 1:number_of_bins
% %     
% %     figure
% %     hold on
% %     plot(final_time{i},Speed{i})
% %     ylabel('Magnitude of Velocity [m/s]')
% %     title('Magnitude of Velocity vs. Time')
% %     NumTicks = 225;
% %     L = get(gca,'XLim');
% %     set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% %     datetick('x','mmm-dd HH:MM PM','keepticks')
% %     set(gca,'XMinorTick','on','YMinorTick','on')
% %     addScrollbar(gca, .6);
% %     
% %     hold off
% %     
% %     i = i + 1;
% % end

% -------------------------------------------------------------------------
% ----------------- Plot of East versus North Velocities ------------------
% -------------------------------------------------------------------------
% Ref: http://depts.washington.edu/nnmrec/docs/20100528_EplerJ_thesis_TidalResourceADCP.pdf

% % i = 1;
% % bin = 1;
% % 
% % for bin = 1:number_of_bins
% %     
% %     figure
% %     hold on
% %     plot(Vel_East{i},Vel_North{i},'*')
% %     ylabel('Velocity North [m/s]')
% %     xlabel('Velocity East [m/s]')
% %     title('North vs. East Velocity')
% %     grid on
% %     hold off
% %     
% %     i = i + 1;
% % end

% -------------------------------------------------------------------------
% ---------------- Plot of East versus Vertical Velocities ----------------
% -------------------------------------------------------------------------
% Ref: http://depts.washington.edu/nnmrec/docs/20100528_EplerJ_thesis_TidalResourceADCP.pdf

% % i = 1;
% % bin = 1;
% % 
% % for bin = 1:number_of_bins
% %     
% %     figure
% %     hold on
% %     plot(Vel_East{i},Vel_Vert{i},'*')
% %     ylabel('Velocity Vertical [m/s]')
% %     xlabel('Velocity East [m/s]')
% %     title('North vs. East Velocity')
% %     grid on
% %     hold off
% %     
% %     i = i + 1;
% % end

% -------------------------------------------------------------------------
% -------- Plot of Vertical Velocity [m/s] vs. Time with scrollbar --------
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % for bin = 1:number_of_bins
% %     
% %     figure
% %     hold on
% %     plot(final_time{i},Vel_Vert{i})
% %     ylabel('Vertical Velocity [m/s]')
% %     title('Vertical Velocity vs. Time')
% %     datetick('x','mmm-dd HH:MM PM','keepticks')
% %     
% %     hold off
% %     
% %     i = i + 1;
% % end

% -------------------------------------------------------------------------
% ---------- Plot of North Velocity [m/s] vs. Time with scrollbar ---------
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % for bin = 1:number_of_bins
% %     
% %     figure
% %     hold on
% %     plot(final_time{i},Vel_North{i})
% %     ylabel('North Velocity [m/s]')
% %     title('North Velocity vs. Time')
% %     datetick('x','mmm-dd HH:MM PM','keepticks')
% %     
% %     hold off
% %     
% %     i = i + 1;
% % end

% -------------------------------------------------------------------------
% -------------- Plot of Magnitude of Velocity with Direction -------------
% -------------------------------------------------------------------------

% % i = 1;
% % bin = 1;
% % 
% % for bin = 1:number_of_bins
% % 
% %         figure
% %         hold on
% %         ax1 = gca;
% %         plot(final_time{i},Speed{i})
% %         ylabel('Magnitude of Velocity [m/s]')
% %         set(gca,'XTickLabel',[])
% %         %plot(final_time{i},Direction{i})
% %         ax2 = axes('Position',get(ax1,'Position'),...
% %                'YAxisLocation','right',...
% %                'Color','none',...
% %                'XColor','k','YColor','r');
% %         linkaxes([ax1 ax2],'x');
% %         hold on
% %         plot(final_time{i},Direction{i},'r+','Parent',ax2);
% %         ylabel('Direction')
% %         title('Velocity and Direction')
% %         NumTicks = 225;
% %         L = get(gca,'XLim');
% %         set(gca,'XTick',linspace(L(1),L(2),NumTicks))
% %         datetick('x','mmm-dd HH:MM PM','keepticks')
% %         addScrollbar(gca, .9);
% %         hold off
% %         
% %     i = i + 1;
% % end

% -------------------------------------------------------------------------
% --------------     Histogram of Energy in Tidal Cycle       -------------
% -------------------------------------------------------------------------


for i=1:number_of_bins
    figure(i+2)
%     Bin_Width(i)=2*iqr(energy{i})/(length(energy{i}.^(1/3)));
%     Total_Plot_Width(i)=max(energy{i})-min(energy{i});
%     Number_of_Hist_Bins(i)=Total_Plot_Width(i)/Bin_Width(i);
%     hist(ebb_energy{i},2*(length(energy{i}))^(1/3));

%     [n2(i),xout2(i)]=hist(ebb_energy{i},2*(length(energy{i}))^(1/3));
%     bar(xout2(i),n2(i),'g');
%     [n3(i),xout3(i)]=hist(flood_energy{i},2*(length(energy{i}))^(1/3));
%     bar(xout3(i),n3(i),'r');
     bar([hist(ebb_energy{i},2*(length(ebb_energy{i}))^(1/3))],'r')
     hold on
     bar([hist(flood_energy{i},2*(length(flood_energy{i}))^(1/3))],'g')
     ylabel('Frequency')
     xlabel('Wh/Cycle?')
     title('Energy per Half Tidal Cycle 12m^2 Turbine')
     legend('Ebb','Flood')
end