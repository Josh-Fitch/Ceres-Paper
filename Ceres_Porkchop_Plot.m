%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% "CERES CREWED MISSION USING DST" PORKCHOP PLOT
%
% Program Description 
% This program utilizes JPL's HORIZONS Ephemeris data files to locate 
% the positions of earth and ceres ffor a user-specified timeframe and 
% uses a lambert solver to find and plot all possible trajectory's in a 
% user-specified number of date steps.
% 
% Function Call:
% Ceres_Porkchop_Plot(date_step, max_timeflight, min_timeflight)
%
% Input Arguments:
% ======================================================================
%    name            :   units    :    description
% ======================================================================
% 1. date_step       :   [days]   :   The size of date intervals to evaluate
% 2. max_timeflight  :   [days]   :   The maximum trajectory duration to be
%                                     evaluated by Lambert solver per launch date
% 3. min_timeflight  :   [days]   :   The minimum trajectory duration to be
%                                     evaluated by Lambert solver per launch date
%
% Output Arguments:
% 1. Pork Chop Contour plot of C3 arrival energy and v_inf departure
%
% Reference Information
%   Author Name:  Joshua Fitch
%   Author Email:  jfitch007@outlook.com
%   Affilliation: JSC USRA Intern Spring 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ceres_Porkchop_Plot(date_step, max_timeflight, min_timeflight)
addpath('C:\Users\jfitc\OneDrive\NASA\JSC Spring 2019\IAC Paper Resources\SPICE Resources\mice');
addpath('C:\Users\jfitc\OneDrive\NASA\JSC Spring 2019\IAC Paper Resources\rodyo-FEX-Lambert-35edc80\rodyo-FEX-Lambert-35edc80');

%% USER INPUTS
% Launch range year and month (data is for 01/01/2035 to 01/01/2050)
launch_min_year = 2035;
launch_min_month = 1;
launch_max_year = 2045;
launch_max_month = 1;
% DEFINES THE GAP BETWEEN AXIS TICKS *IN YEARS*
x_axis_spacing = 1.2; 
y_axis_spacing = 1.6;

%% ESTABLISH DATE RANGES
% Convert dates in MJD
launch_min_date = julian(launch_min_month,1,launch_min_year) - 2400000.5;
launch_max_date = julian(launch_max_month,1,launch_max_year) - 2400000.5;
arrival_min_date = launch_min_date + min_timeflight;
arrival_max_date = launch_max_date + max_timeflight;

return_launch_min_date = launch_min_date;
return_launch_max_date = launch_max_date;
return_arrival_min_date = return_launch_min_date + min_timeflight;
return_arrival_max_date = return_launch_max_date + max_timeflight;

%% HORIZONS Data
% Ceres, every day data from 2035 to 2060
eph_file_Ceres = 'results_Ceres.orb';
orb_eph_Ceres = load(eph_file_Ceres);

% Earth, every day data from 2035 to 2060
eph_file_Earth = 'results_Earth.orb';
orb_eph_Earth = load(eph_file_Earth);

%% Date Range Initialization
% Form list of dates.
launch_dates = (launch_min_date:date_step:launch_max_date);
arrival_dates = (arrival_min_date:date_step:arrival_max_date);
return_launch_dates = (return_launch_min_date:date_step:return_launch_max_date);
return_arrival_dates = (return_arrival_min_date:date_step:return_arrival_max_date);
% Store the number of dates.
n_launch_dates = length(launch_dates);
n_arrival_dates = length(arrival_dates);
n_return_launch = length(return_launch_dates);
n_return_arrival = length(return_arrival_dates);

%% INITIALIZE CERES AND EARTH POSITION AND VELOCITY VECTORS
% Loop over all launch dates and pull relevant position and velocity data
for k=1:n_launch_dates
    % HORIZONS
    [~, ind1] = ismember(launch_dates(k), orb_eph_Earth(:,1));
    r1 = orb_eph_Earth(ind1,2:7);
    % Position Vector Assignments
    earth_x(k) = r1(1,1);  % position of earth in km w.r.t. the Sun
    earth_y(k) = r1(1,2);
    earth_z(k) = r1(1,3);
    % Velocity Vector Assignments
    earth_vx(k) = r1(1,4);  % velocity of earth in km w.r.t. the Sun
    earth_vy(k) = r1(1,5);
    earth_vz(k) = r1(1,6);
end

% Loop over all arrival dates and pull relevant position and velocity data
for k=1:n_arrival_dates
    % HORIZONS
    [~, ind2] = ismember(arrival_dates(k), orb_eph_Ceres(:,1));
    r2 = orb_eph_Ceres(ind2,2:7);
    % Position Vector Assignments
    ceres_x(k) = r2(1,1);  % position of ceres in km w.r.t. the Sun
    ceres_y(k) = r2(1,2);
    ceres_z(k) = r2(1,3);
    % Velocity Vector Assignments
    ceres_vx(k) = r2(1,4);  % velocity of ceres in km w.r.t. the Sun
    ceres_vy(k) = r2(1,5);
    ceres_vz(k) = r2(1,6);
end

% Loop over all return launch dates and pull relevant position and velocity data
for k=1:n_return_launch
    % HORIZONS
    [~, ind3] = ismember(return_launch_dates(k), orb_eph_Ceres(:,1));
    r3 = orb_eph_Ceres(ind3,2:7);
    % Position Vector Assignments
    ceres_x_return(k) = r3(1,1);
    ceres_y_return(k) = r3(1,2);
    ceres_z_return(k) = r3(1,3);
    % Velocity Vector Assignments
    ceres_vx_return(k) = r3(1,4);
    ceres_vy_return(k) = r3(1,5);
    ceres_vz_return(k) = r3(1,6);
end

% Loop over all return arrival dates and pull relevant position and velocity data
for k=1:n_return_arrival
    % HORIZONS
    [~, ind4] = ismember(return_arrival_dates(k), orb_eph_Earth(:,1));
    r4 = orb_eph_Earth(ind4,2:7);
    % Position Vector Assignments
    earth_x_return(k) = r4(1,1);
    earth_y_return(k) = r4(1,2);
    earth_z_return(k) = r4(1,3);
    % Velocity Vector Assignments
    earth_vx_return(k) = r4(1,4);
    earth_vy_return(k) = r4(1,5);
    earth_vz_return(k) = r4(1,6);
end

%% LAMBERT SOLVER
GM_sun = 1.327e11;
for k = 1:n_launch_dates
    r1_vectors(k,:) = [earth_x(k),earth_y(k),earth_z(k)];
end
for k = 1:n_arrival_dates
    r2_vectors(k,:) = [ceres_x(k),ceres_y(k),ceres_z(k)];
end
for k = 1:n_return_launch
    r3_vectors(k,:) = [ceres_x_return(k),ceres_y_return(k),ceres_z_return(k)];
end
for k = 1:n_return_arrival
    r4_vectors(k,:) = [earth_x_return(k),earth_y_return(k),earth_z_return(k)];
end

% POSITIVE TIME OF FLIGHTS
for k = 1:n_launch_dates
    r1 = r1_vectors(k,:);
    for j = k:n_arrival_dates
        time_flight(k,j) = arrival_dates(j) - launch_dates(k);
        if time_flight(k,j) > max_timeflight || time_flight(k,j) < min_timeflight
            v1(k,j,:) = [0,0,0];
            v2(k,j,:) = [0,0,0];
            continue
        end
        r2 = r2_vectors(j,:);
        % Use Lambert Solver
        [v1(k,j,:), v2(k,j,:), ~, ~] = lambert(r1(1,:), r2(1,:), time_flight(k,j), 0, GM_sun);
    end
end

% NEGATIVE TIME OF FLIGHTS
for k = 1:n_launch_dates
    r1 = r1_vectors(k,:);
    for j = k:n_arrival_dates
        time_flight_neg(k,j) = -(arrival_dates(j) - launch_dates(k));
        if -time_flight_neg(k,j) > max_timeflight || -time_flight_neg(k,j) < min_timeflight
            v1_neg(k,j,:) = [0,0,0];
            v2_neg(k,j,:) = [0,0,0];
            continue
        end
        r2 = r2_vectors(j,:);
        % Use Lambert Solver
        [v1_neg(k,j,:), v2_neg(k,j,:), ~, ~] = lambert(r1(1,:), r2(1,:), time_flight_neg(k,j), 0, GM_sun);
    end
end

% RETURN FLIGHT - POSITIVE TIME OF FLIGHTS
for k = 1:n_return_launch
    r3 = r3_vectors(k,:);
    for j = k:n_return_arrival
        time_flight_back(k,j) = return_arrival_dates(j) - return_launch_dates(k);
        if time_flight_back(k,j) > max_timeflight || time_flight_back(k,j) < min_timeflight
            v3(k,j,:) = [0,0,0];
            v4(k,j,:) = [0,0,0];
            continue
        end
        r4 = r4_vectors(j,:);
        % Use Lambert Solver
        [v3(k,j,:), v4(k,j,:), ~, ~] = lambert(r3(1,:), r4(1,:), time_flight_back(k,j), 0, GM_sun);
    end
end

% RETURN FLIGHT - NEGATIVE TIME OF FLIGHTS
for k = 1:n_return_launch
    r3 = r3_vectors(k,:);
    for j = k:n_return_arrival
        time_flight_back_neg(k,j) = -(return_arrival_dates(j) - return_launch_dates(k));
        if -time_flight_back_neg(k,j) > max_timeflight || -time_flight_back_neg(k,j) < min_timeflight
            v3_neg(k,j,:) = [0,0,0];
            v4_neg(k,j,:) = [0,0,0];
            continue
        end
        r4 = r4_vectors(j,:);
        % Use Lambert Solver
        [v3_neg(k,j,:), v4_neg(k,j,:), ~, ~] = lambert(r3(1,:), r4(1,:), time_flight_back_neg(k,j), 0, GM_sun);
    end
end

%% CALCULATE DELTA V

% POSITIVE TIME OF FLIGHTS
for k = 1:n_launch_dates
    for j = k:n_arrival_dates
        if time_flight(k,j) > max_timeflight || time_flight(k,j) < min_timeflight
            dV_earth_vector(k,j,1) = 0;
            dV_earth_vector(k,j,2) = 0;
            dV_earth_vector(k,j,3) = 0;
            dV_ceres_vector(k,j,1) = 0;
            dV_ceres_vector(k,j,2) = 0;
            dV_ceres_vector(k,j,3) = 0;
            dV_earth(k,j) = 0;
            dV_ceres(k,j) = 0;
            total_dV(k,j) = 0;
            c3_E(k,j) = 0;
            v_inf(k,j) = 0;
            continue
        end
        dV_earth_vector(k,j,1) = v1(k,j,1) - earth_vx(k);
        dV_earth_vector(k,j,2) = v1(k,j,2) - earth_vy(k);
        dV_earth_vector(k,j,3) = v1(k,j,3) - earth_vz(k);
        dV_ceres_vector(k,j,1) = v2(k,j,1) - ceres_vx(j);
        dV_ceres_vector(k,j,2) = v2(k,j,2) - ceres_vy(j);
        dV_ceres_vector(k,j,3) = v2(k,j,3) - ceres_vz(j);
        dV_earth(k,j) = sqrt(dV_earth_vector(k,j,1)^2+dV_earth_vector(k,j,2)^2+dV_earth_vector(k,j,3)^2);
        dV_ceres(k,j) = sqrt(dV_ceres_vector(k,j,1)^2+dV_ceres_vector(k,j,2)^2+dV_ceres_vector(k,j,3)^2);
        total_dV(k,j) = dV_earth(k,j) + dV_ceres(k,j);
        v_inf(k,j) = dV_ceres(k,j);
        c3_E(k,j) = (dV_earth(k,j))^2;
    end
end

% NEGATIVE TIME OF FLIGHTS
for k = 1:n_launch_dates
    for j = k:n_arrival_dates
        if -time_flight_neg(k,j) > max_timeflight || -time_flight_neg(k,j) < min_timeflight
            dV_earth_vector_neg(k,j,1) = 0;
            dV_earth_vector_neg(k,j,2) = 0;
            dV_earth_vector_neg(k,j,3) = 0;
            dV_ceres_vector_neg(k,j,1) = 0;
            dV_ceres_vector_neg(k,j,2) = 0;
            dV_ceres_vector_neg(k,j,3) = 0;
            dV_earth_neg(k,j) = 0;
            dV_ceres_neg(k,j) = 0;
            total_dV_neg(k,j) = 0;
            c3_E_neg(k,j) = 0;
            v_inf_neg(k,j) = 0;
            continue
        end
        dV_earth_vector_neg(k,j,1) = v1_neg(k,j,1) - earth_vx(k);
        dV_earth_vector_neg(k,j,2) = v1_neg(k,j,2) - earth_vy(k);
        dV_earth_vector_neg(k,j,3) = v1_neg(k,j,3) - earth_vz(k);
        dV_ceres_vector_neg(k,j,1) = v2_neg(k,j,1) - ceres_vx(j);
        dV_ceres_vector_neg(k,j,2) = v2_neg(k,j,2) - ceres_vy(j);
        dV_ceres_vector_neg(k,j,3) = v2_neg(k,j,3) - ceres_vz(j);
        dV_earth_neg(k,j) = sqrt(dV_earth_vector_neg(k,j,1)^2+dV_earth_vector_neg(k,j,2)^2+dV_earth_vector_neg(k,j,3)^2);
        dV_ceres_neg(k,j) = sqrt(dV_ceres_vector_neg(k,j,1)^2+dV_ceres_vector_neg(k,j,2)^2+dV_ceres_vector_neg(k,j,3)^2);
        total_dV_neg(k,j) = dV_earth_neg(k,j) + dV_ceres_neg(k,j);
        v_inf_neg(k,j) = dV_ceres_neg(k,j);
        c3_E_neg(k,j) = (dV_earth_neg(k,j))^2;
    end
end

% RETURN FLIGHT - POSITIVE TIME OF FLIGHTS
for k = 1:n_return_launch
    for j = k:n_return_arrival
        if time_flight_back(k,j) > max_timeflight || time_flight_back(k,j) < min_timeflight
            dV_earth_vector_return(k,j,1) = 0;
            dV_earth_vector_return(k,j,2) = 0;
            dV_earth_vector_return(k,j,3) = 0;
            dV_ceres_vector_return(k,j,1) = 0;
            dV_ceres_vector_return(k,j,2) = 0;
            dV_ceres_vector_return(k,j,3) = 0;
            dV_earth_return(k,j) = 0;
            dV_ceres_return(k,j) = 0;
            total_dV_return(k,j) = 0;
            c3_E_return(k,j) = 0;
            v_inf_return(k,j) = 0;
            continue
        end
        dV_earth_vector_return(k,j,1) = v4(k,j,1) - earth_vx_return(j);
        dV_earth_vector_return(k,j,2) = v4(k,j,2) - earth_vy_return(j);
        dV_earth_vector_return(k,j,3) = v4(k,j,3) - earth_vz_return(j);
        dV_ceres_vector_return(k,j,1) = v3(k,j,1) - ceres_vx_return(k);
        dV_ceres_vector_return(k,j,2) = v3(k,j,2) - ceres_vy_return(k);
        dV_ceres_vector_return(k,j,3) = v3(k,j,3) - ceres_vz_return(k);
        dV_earth_return(k,j) = sqrt(dV_earth_vector_return(k,j,1)^2+dV_earth_vector_return(k,j,2)^2+dV_earth_vector_return(k,j,3)^2);
        dV_ceres_return(k,j) = sqrt(dV_ceres_vector_return(k,j,1)^2+dV_ceres_vector_return(k,j,2)^2+dV_ceres_vector_return(k,j,3)^2);
        total_dV_return(k,j) = dV_earth_return(k,j) + dV_ceres_return(k,j);
        v_inf_return(k,j) = dV_earth_return(k,j);
        c3_E_return(k,j) = (dV_ceres_return(k,j))^2;
    end
end

% RETURN FLIGHT - NEGATIVE TIME OF FLIGHTS
for k = 1:n_return_launch
    for j = k:n_return_arrival
        if -time_flight_back_neg(k,j) > max_timeflight || -time_flight_back_neg(k,j) < min_timeflight
            dV_earth_vector_return_neg(k,j,1) = 0;
            dV_earth_vector_return_neg(k,j,2) = 0;
            dV_earth_vector_return_neg(k,j,3) = 0;
            dV_ceres_vector_return_neg(k,j,1) = 0;
            dV_ceres_vector_return_neg(k,j,2) = 0;
            dV_ceres_vector_return_neg(k,j,3) = 0;
            dV_earth_return_neg(k,j) = 0;
            dV_ceres_return_neg(k,j) = 0;
            total_dV_return_neg(k,j) = 0;
            c3_E_return_neg(k,j) = 0;
            v_inf_return_neg(k,j) = 0;
            continue
        end
        dV_earth_vector_return_neg(k,j,1) = v4_neg(k,j,1) - earth_vx_return(j);
        dV_earth_vector_return_neg(k,j,2) = v4_neg(k,j,2) - earth_vy_return(j);
        dV_earth_vector_return_neg(k,j,3) = v4_neg(k,j,3) - earth_vz_return(j);
        dV_ceres_vector_return_neg(k,j,1) = v3_neg(k,j,1) - ceres_vx_return(k);
        dV_ceres_vector_return_neg(k,j,2) = v3_neg(k,j,2) - ceres_vy_return(k);
        dV_ceres_vector_return_neg(k,j,3) = v3_neg(k,j,3) - ceres_vz_return(k);
        dV_earth_return_neg(k,j) = sqrt(dV_earth_vector_return_neg(k,j,1)^2+dV_earth_vector_return_neg(k,j,2)^2+dV_earth_vector_return_neg(k,j,3)^2);
        dV_ceres_return_neg(k,j) = sqrt(dV_ceres_vector_return_neg(k,j,1)^2+dV_ceres_vector_return_neg(k,j,2)^2+dV_ceres_vector_return_neg(k,j,3)^2);
        total_dV_return_neg(k,j) = dV_earth_return_neg(k,j) + dV_ceres_return_neg(k,j);
        v_inf_return_neg(k,j) = dV_earth_return_neg(k,j);
        c3_E_return_neg(k,j) = (dV_ceres_return_neg(k,j))^2;
    end
end

%% PLOTTING
figure
%% Figure 1 - TRIP TO CERES
pos_index = 1;
for MJD = launch_dates
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    launch_dates_plotting1(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end
pos_index = 1;
for MJD = arrival_dates
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    arrival_dates_plotting1(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end
hold on;
[C,h1] = contour(launch_dates_plotting1,arrival_dates_plotting1,c3_E.',[50:50:150,200],'r','ShowText','off');
[C,h2] = contour(launch_dates_plotting1,arrival_dates_plotting1,v_inf.',0:2:10,'b','ShowText','off');
%contour(launch_dates,arrival_dates,time_flight.',200:200:max_timeflight-100,'k','on')
[C,h3] = contour(launch_dates_plotting1,arrival_dates_plotting1,c3_E_neg.',[50:50:150,200],'r','ShowText','off');
[C,h4] = contour(launch_dates_plotting1,arrival_dates_plotting1,v_inf_neg.',0:2:10,'b','ShowText','off');
hold off;

% PLOT DETAILS
title("Trip to Ceres for "+launch_min_year+" - "+launch_max_year+" Launch Dates")
xlabel('Launch Date [year/month/day]')
ylabel('Arrival Date [year/month/day]')
legend('C3 Energy (km^2/s^2)','Arrival Vinfinity (km/s)','TOF (days)');
legend('location','northwest');

% FORMATTED LAUNCH DATES
pos_index = 1;
for MJD = launch_min_date:366*x_axis_spacing:launch_max_date
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    formatted_launch_dates(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end

% FORMATTED ARRIVAL DATES
pos_index = 1;
for MJD = arrival_min_date:366*y_axis_spacing:arrival_max_date
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    formatted_arrival_dates(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end

% PLOT FORMATTING
a = gca;
a.FontSize = 20;
a.FontWeight = 'bold';
h1.LineWidth = 2;
h2.LineWidth = 2;
h3.LineWidth = 2;
h4.LineWidth = 2;
xticks(formatted_launch_dates);
yticks(formatted_arrival_dates);
datetick('x',26,'keepticks');
datetick('y',26,'keepticks');
xtickangle(45);
ytickangle(0);

figure
%% Figure 2 - TRIP BACK FROM CERES
pos_index = 1;
for MJD = return_launch_dates
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    return_launch_dates_plotting1(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end
pos_index = 1;
for MJD = return_arrival_dates
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    return_arrival_dates_plotting1(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end
hold on;
contour(return_launch_dates_plotting1,return_arrival_dates_plotting1,c3_E_return.',[25,50,75,100],'r','ShowText','off')
contour(return_launch_dates_plotting1,return_arrival_dates_plotting1,v_inf_return.',0:2:10,'b','ShowText','off')
%contour(return_launch_dates,return_arrival_dates,time_flight_back.',200:200:max_timeflight-100,'k','ShowText','on')
contour(return_launch_dates_plotting1,return_arrival_dates_plotting1,c3_E_return_neg.',[25,50,75,100],'r','ShowText','off')
contour(return_launch_dates_plotting1,return_arrival_dates_plotting1,v_inf_return_neg.',0:2:10,'b','ShowText','off')
hold off;

% PLOT DETAILS
title("Return Trip - Ceres to Earth")
xlabel('Return Launch Date [year/month/day]')
ylabel('Return Arrival Date [year/month/day]')
legend('C3 Energy (Return) (km^2/s^2)','Earth Arrival Vinfinity (km/s)','TOF (days)');
legend('location','northwest');

% FORMATTED LAUNCH DATES
pos_index = 1;
for MJD = return_launch_min_date:366*x_axis_spacing:return_launch_max_date
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    formatted_return_launch_dates(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end

% FORMATTED ARRIVAL DATES
pos_index = 1;
for MJD = return_arrival_min_date:366*y_axis_spacing:return_arrival_max_date
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    formatted_return_arrival_dates(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end

% PLOT FORMATTING
a = gca;
a.FontSize = 20;
a.FontWeight = 'bold';
xticks(formatted_return_launch_dates);
yticks(formatted_return_arrival_dates);
datetick('x',26,'keepticks');
datetick('y',26,'keepticks');
xtickangle(45);
ytickangle(0);

figure
%% FIGURE 3 - TOTAL DELTA V 
hold on;
contourf(launch_dates_plotting1,arrival_dates_plotting1,total_dV.',[5:5:30,40])
contourf(launch_dates_plotting1,arrival_dates_plotting1,total_dV_neg.',[5:5:30,40])
%contour(launch_dates,arrival_dates,time_flight.',200:200:max_timeflight-100,'k','ShowText','on')
hold off;

% PLOT DETAILS
title("Trip to Ceres for "+launch_min_year+" - "+launch_max_year+" Launch Dates")
xlabel('Launch Date [year/month/day]')
ylabel('Arrival Date [year/month/day]')
%legend('C3 Energy (km^2/s^2)','Arrival Vinfinity (km/s)','TOF (days)');
%legend('location','northwest');

% FORMATTED LAUNCH DATES
pos_index = 1;
for MJD = launch_min_date:366*x_axis_spacing:launch_max_date
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    formatted_launch_dates(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end

% FORMATTED ARRIVAL DATES
pos_index = 1;
for MJD = arrival_min_date:366*y_axis_spacing:arrival_max_date
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    formatted_arrival_dates(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end

% PLOT FORMATTING
a = gca;
a.FontSize = 20;
a.FontWeight = 'bold';
xticks(formatted_launch_dates);
yticks(formatted_arrival_dates);
datetick('x',26,'keepticks');
datetick('y',26,'keepticks');
xtickangle(45);
ytickangle(0);

figure
%% FIGURE 4 - Values Plotted Against Launch Date and TOF (Earth to Ceres)
tf_plot_axis = min_timeflight:date_step:max_timeflight;
for i = 1:length(launch_dates)
   for j = i:length(arrival_dates)
       time = arrival_dates(j) - launch_dates(i);
       flight_index = (time + 10 - min_timeflight)/date_step;
       if flight_index < 132
           c3_E_tof(i,flight_index) = c3_E(i,j);
           v_inf_tof(i,flight_index) = v_inf(i,j);
           c3_E_neg_tof(i,flight_index) = c3_E_neg(i,j);
           v_inf_neg_tof(i,flight_index) = v_inf_neg(i,j);
       end
   end    
end

pos_index = 1;
for MJD = launch_dates
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    launch_dates_plotting2(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end

hold on;
[~,h1] = contour(launch_dates_plotting2,tf_plot_axis,c3_E_tof',[50:50:150,200],'r','ShowText','off');
[~,h2] = contour(launch_dates_plotting2,tf_plot_axis,v_inf_tof',0:1:10,'b','ShowText','off');
[~,h3] = contour(launch_dates_plotting2,tf_plot_axis,c3_E_neg_tof',[50:50:150,200],'r','ShowText','off');
[~,h4] = contour(launch_dates_plotting2,tf_plot_axis,v_inf_neg_tof',0:1:10,'b','ShowText','off');
hold off;
%{

%}
% PLOT DETAILS
title("Earth to Ceres: "+launch_min_year+" - "+launch_max_year+" Launch Dates")
xlabel('Launch Date [year/month/day]')
ylabel('Time of Flight (TOF) [days]')
legend('C3 Energy (km^2/s^2)','Arrival Vinfinity (km/s)');
legend('location','northwest');

% FORMATTED LAUNCH DATES
pos_index = 1;
for MJD = launch_min_date:365.5:launch_max_date
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    formatted_launch_dates(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end

% PLOT FORMATTING
a = gca;
a.FontSize = 18;
a.FontWeight = 'bold';
h1.LineWidth = 1.3;
h2.LineWidth = 1.3;
h3.LineWidth = 1.3;
h4.LineWidth = 1.3;
xticks(formatted_launch_dates);
datetick('x',26,'keepticks');
xtickangle(45);
ytickangle(0);
set(a, 'xgrid', 'on')
set(a, 'ygrid', 'on')
set(a,'GridLineStyle','- -')
a.LineWidth = 1.3;

figure
%% FIGURE 5 - Values Plotted Against Launch Date and TOF (Ceres to Earth)
tf_plot_axis = min_timeflight:date_step:max_timeflight;
for i = 1:length(return_launch_dates)
   for j = i:length(return_arrival_dates)
       time = return_arrival_dates(j) - return_launch_dates(i);
       flight_index = (time + 10 - min_timeflight)/date_step;
       if flight_index < 132
           c3_E_return_tof(i,flight_index) = c3_E_return(i,j);
           v_inf_return_tof(i,flight_index) = v_inf_return(i,j);
           c3_E_return_neg_tof(i,flight_index) = c3_E_return_neg(i,j);
           v_inf_return_neg_tof(i,flight_index) = v_inf_return_neg(i,j);
       end
   end    
end

pos_index = 1;
for MJD = return_launch_dates
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    launch_return_dates_plotting2(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end

hold on;
[~,h1] = contour(launch_return_dates_plotting2,tf_plot_axis,c3_E_return_tof',50:25:175,'r','ShowText','off');
[~,h2] = contour(launch_return_dates_plotting2,tf_plot_axis,v_inf_return_tof',0:1:12,'b','ShowText','off');
[~,h3] = contour(launch_return_dates_plotting2,tf_plot_axis,c3_E_return_neg_tof',50:50:150,'r','ShowText','off');
[~,h4] = contour(launch_return_dates_plotting2,tf_plot_axis,v_inf_return_neg_tof',0:1:12,'b','ShowText','off');
hold off;
%{

%}
% PLOT DETAILS
title("Ceres to Earth: "+launch_min_year+" - "+launch_max_year+" Return Launch Dates")
xlabel('Launch Date [year/month/day]')
ylabel('Time of Flight (TOF) [days]')
legend('C3 Energy (km^2/s^2)','Arrival Vinfinity (km/s)');
legend('location','northwest');

% FORMATTED LAUNCH DATES
pos_index = 1;
for MJD = return_launch_min_date:365.5:return_launch_max_date
    [year,month,day,~,~,~,~,~] = julian2greg(MJD + 2400000.5);
    formatted_launch_return_dates(pos_index) = datenum(year,month,day);
    pos_index = pos_index + 1;
end

% PLOT FORMATTING
a = gca;
a.FontSize = 18;
a.FontWeight = 'bold';
h1.LineWidth = 1.3;
h2.LineWidth = 1.3;
h3.LineWidth = 1.3;
h4.LineWidth = 1.3;
xticks(formatted_launch_return_dates);
datetick('x',26,'keepticks');
xtickangle(45);
ytickangle(0);
set(a, 'xgrid', 'on')
set(a, 'ygrid', 'on')
set(a,'GridLineStyle','- -')
a.LineWidth = 1.3;

%% PULL RELEVANT DATA FOR CHECKING
%{
% set all zeros to NaN to avoid zeros being used as minimums

v_inf(v_inf==0) = NaN;
v_inf_return(v_inf_return==0) = NaN; 
v_inf_neg(v_inf_neg==0) = NaN;
v_inf_return_neg(v_inf_return_neg==0) = NaN;
c3_E(c3_E==0) = NaN;
c3_E_return(c3_E_return==0) = NaN;
c3_E_neg(c3_E_neg==0) = NaN;
c3_E_return_neg(c3_E_return_neg==0) = NaN;
total_dV(total_dV==0) = NaN;
total_dV_return(total_dV_return==0) = NaN;
total_dV_neg(total_dV_neg==0) = NaN;
total_dV_return_neg(total_dV_return_neg==0) = NaN;

% Go through and find mins and associated values

[c3_min_c3,I] = min(c3_E(:));
[c3_min_c3_neg,I2] = min(c3_E_neg(:));
if (c3_min_c3 < c3_min_c3_neg)
    [I_row,I_col] = ind2sub(size(c3_E),I);
    c3_min_vinf = v_inf(I_row,I_col);
    c3_min_deltav = total_dV(I_row,I_col); 
else
    c3_min_c3 = c3_min_c3_neg;
    [I_row,I_col] = ind2sub(size(c3_E_neg),I2);
    c3_min_vinf = v_inf_neg(I_row,I_col);
    c3_min_deltav = total_dV_neg(I_row,I_col);
end
c3_min_tof = time_flight(I_row,I_col);
c3_min_launch = launch_dates(I_row) - 693960;
c3_min_arrival = arrival_dates(I_col) - 693960;
min_c3_values = [c3_min_c3,c3_min_vinf,c3_min_deltav,c3_min_tof,c3_min_launch,c3_min_arrival];
xlswrite('C:\Users\jfitc\OneDrive\NASA\JSC Spring 2019\IAC Paper Resources\Ceres Porkchop Data.xlsx',min_c3_values,'Sheet1','B3')

[c3_min_return_c3,I] = min(c3_E_return(:));
[c3_min_c3_return_neg,I2] = min(c3_E_return_neg(:));
if (c3_min_return_c3 < c3_min_c3_return_neg)
    [I_row,I_col] = ind2sub(size(c3_E_return),I);
    c3_min_return_vinf = v_inf_return(I_row,I_col);
    c3_min_return_deltav = total_dV_return(I_row,I_col); 
else
    c3_min_return_c3 = c3_min_c3_return_neg;
    [I_row,I_col] = ind2sub(size(c3_E_return_neg),I2);
    c3_min_return_vinf = v_inf_return_neg(I_row,I_col);
    c3_min_return_deltav = total_dV_return_neg(I_row,I_col); 
end
c3_min_return_tof = time_flight_back(I_row,I_col);
c3_min_return_launch = return_launch_dates(I_row) - 693960;
c3_min_return_arrival = return_arrival_dates(I_col) - 693960;
min_c3_return_values = [c3_min_return_c3,c3_min_return_vinf,c3_min_return_deltav,c3_min_return_tof,c3_min_return_launch,c3_min_return_arrival];
xlswrite('C:\Users\jfitc\OneDrive\NASA\JSC Spring 2019\IAC Paper Resources\Ceres Porkchop Data.xlsx',min_c3_return_values,'Sheet1','B8')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[vinf_min_vinf,I] = min(v_inf(:));
[vinf_min_vinf_neg,I2] = min(v_inf_neg(:));
if (vinf_min_vinf < vinf_min_vinf_neg)
    [I_row,I_col] = ind2sub(size(v_inf),I);
    vinf_min_c3 = c3_E(I_row,I_col);
    vinf_min_deltav = total_dV(I_row,I_col);
else
    vinf_min_vinf = vinf_min_vinf_neg;
    [I_row,I_col] = ind2sub(size(v_inf_neg),I2);
    vinf_min_c3 = c3_E_neg(I_row,I_col);
    vinf_min_deltav = total_dV_neg(I_row,I_col);
end
vinf_min_tof = time_flight(I_row,I_col);
vinf_min_launch = launch_dates(I_row)- 693960;
vinf_min_arrival = arrival_dates(I_col)- 693960;
min_vinf_values = [vinf_min_c3,vinf_min_vinf,vinf_min_deltav,vinf_min_tof,vinf_min_launch,vinf_min_arrival];
xlswrite('C:\Users\jfitc\OneDrive\NASA\JSC Spring 2019\IAC Paper Resources\Ceres Porkchop Data.xlsx',min_vinf_values,'Sheet1','B4')

[vinf_min_return_vinf,I] = min(v_inf_return(:));
[vinf_min_return_vinf_neg,I2] = min(v_inf_return_neg(:));
if (vinf_min_return_vinf < vinf_min_return_vinf_neg)
    [I_row,I_col] = ind2sub(size(v_inf_return),I);
    vinf_min_return_c3 = c3_E_return(I_row,I_col);
    vinf_min_return_deltav = total_dV_return(I_row,I_col);
else
    vinf_min_return_vinf = vinf_min_return_vinf_neg;
    [I_row,I_col] = ind2sub(size(v_inf_return_neg),I2);
    vinf_min_return_c3 = c3_E_return_neg(I_row,I_col);
    vinf_min_return_deltav = total_dV_return_neg(I_row,I_col);
end
vinf_min_return_tof = time_flight_back(I_row,I_col);
vinf_min_return_launch = return_launch_dates(I_row) - 693960;
vinf_min_return_arrival = return_arrival_dates(I_col) - 693960;
min_vinf_return_values = [vinf_min_return_c3,vinf_min_return_vinf,vinf_min_return_deltav,vinf_min_return_tof,vinf_min_return_launch,vinf_min_return_arrival];
xlswrite('C:\Users\jfitc\OneDrive\NASA\JSC Spring 2019\IAC Paper Resources\Ceres Porkchop Data.xlsx',min_vinf_return_values,'Sheet1','B9')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[dv_min_dv,I] = min(total_dV(:));
[dv_min_dv_neg,I2] = min(total_dV_neg(:));
if (dv_min_dv < dv_min_dv_neg)
    [I_row_dv,I_col_dv] = ind2sub(size(total_dV),I);
    dv_min_c3 = c3_E(I_row_dv,I_col_dv);
    dv_min_vinf = v_inf(I_row_dv,I_col_dv); 
else
    dv_min_dv = dv_min_dv_neg;
    [I_row_dv,I_col_dv] = ind2sub(size(total_dV_neg),I2);
    dv_min_c3 = c3_E_neg(I_row_dv,I_col_dv);
    dv_min_vinf = v_inf_neg(I_row_dv,I_col_dv); 
end
dv_min_tof = time_flight(I_row_dv,I_col_dv);
dv_min_launch = launch_dates(I_row_dv) - 693960;
dv_min_arrival = arrival_dates(I_col_dv) - 693960;
min_dv_values = [dv_min_c3,dv_min_vinf,dv_min_dv,dv_min_tof,dv_min_launch,dv_min_arrival];
xlswrite('C:\Users\jfitc\OneDrive\NASA\JSC Spring 2019\IAC Paper Resources\Ceres Porkchop Data.xlsx',min_dv_values,'Sheet1','B5')

[dv_min_return_dv,I] = min(total_dV_return(:));
[dv_min_return_dv_neg,I2] = min(total_dV_return_neg(:));
if (dv_min_return_dv < dv_min_return_dv_neg)
    [I_row_dv_b,I_col_dv_b] = ind2sub(size(total_dV_return),I);
    dv_min_return_c3 = c3_E_return(I_row_dv_b,I_col_dv_b);
    dv_min_return_vinf = v_inf_return(I_row_dv_b,I_col_dv_b); 
else
    dv_min_return_dv = dv_min_return_dv_neg;
    [I_row_dv_b,I_col_dv_b] = ind2sub(size(total_dV_return_neg),I2);
    dv_min_return_c3 = c3_E_return_neg(I_row_dv_b,I_col_dv_b);
    dv_min_return_vinf = v_inf_return_neg(I_row_dv_b,I_col_dv_b); 
end
dv_min_return_tof = time_flight_back(I_row_dv_b,I_col_dv_b);
dv_min_return_launch = launch_dates(I_row_dv_b) - 693960;
dv_min_return_arrival = arrival_dates(I_col_dv_b) - 693960;
min_dv_return_values = [dv_min_return_c3,dv_min_return_vinf,dv_min_return_dv,dv_min_return_tof,dv_min_return_launch,dv_min_return_arrival];
xlswrite('C:\Users\jfitc\OneDrive\NASA\JSC Spring 2019\IAC Paper Resources\Ceres Porkchop Data.xlsx',min_dv_return_values,'Sheet1','B10')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Total Mission Values

max_allow_mission_time = 1200;
min_desired_ceres_stay = 0;
total_dV_return_test1 = total_dV_return;
total_dV_return_neg_test1 = total_dV_return_neg;
total_dV_test2 = total_dV;
total_dV_neg_test2 = total_dV_neg;
total_dV_test3 = total_dV;
total_dV_return_test3 = total_dV_return;
total_dV_neg_test3 = total_dV_neg;
total_dV_return_neg_test3 = total_dV_return_neg;

% Evaluate return trips after min trip there
min_dv_set1_vals = zeros(4,14);
for case_num = 1:4
    min_mission_dv = 10000;
    for launch_index = I_col_dv:n_launch_dates
       for arrival_index = (I_col_dv + min_timeflight/date_step):n_arrival_dates
           test_mission_dv = dv_min_dv + total_dV_return_test1(launch_index,arrival_index);
           mission_time = return_arrival_dates(arrival_index) - (dv_min_launch + 693960);
           ceres_stay = return_launch_dates(launch_index) - (dv_min_arrival + 693960);
           if (test_mission_dv < min_mission_dv) && (mission_time < max_allow_mission_time) && (ceres_stay > min_desired_ceres_stay)
               min_mission_dv = test_mission_dv;
               dVmin_I = [launch_index,arrival_index];
               is_neg = 0;
           end
           test_mission_dv_neg = dv_min_dv + total_dV_return_neg_test1(launch_index,arrival_index);
           if (test_mission_dv_neg < min_mission_dv) && (mission_time < max_allow_mission_time) && (ceres_stay > min_desired_ceres_stay)
               min_mission_dv = test_mission_dv_neg;
               dVmin_I = [launch_index,arrival_index];
               is_neg = 1;
           end
       end
    end
    min_dv_set1_vals(case_num,1) = min_mission_dv;
    min_dv_set1_vals(case_num,2) = dv_min_dv;
    min_dv_set1_vals(case_num,4) = dv_min_c3;
    min_dv_set1_vals(case_num,5) = dv_min_vinf;
    if is_neg
        min_dv_set1_vals(case_num,3) = total_dV_return_neg_test1(dVmin_I(1),dVmin_I(2));
        min_dv_set1_vals(case_num,6) = c3_E_return_neg(dVmin_I(1),dVmin_I(2));
        min_dv_set1_vals(case_num,7) = v_inf_return_neg(dVmin_I(1),dVmin_I(2));
    else 
        min_dv_set1_vals(case_num,3) = total_dV_return_test1(dVmin_I(1),dVmin_I(2));
        min_dv_set1_vals(case_num,6) = c3_E_return(dVmin_I(1),dVmin_I(2));
        min_dv_set1_vals(case_num,7) = v_inf_return(dVmin_I(1),dVmin_I(2));
    end
    min_dv_set1_vals(case_num,8) = dv_min_tof;
    min_dv_set1_vals(case_num,9) = time_flight_back(dVmin_I(1),dVmin_I(2));
    min_dv_set1_vals(case_num,10) = dv_min_launch;
    min_dv_set1_vals(case_num,11) = dv_min_arrival;
    min_dv_set1_vals(case_num,12) = return_launch_dates(dVmin_I(1))-693960;
    min_dv_set1_vals(case_num,13) = return_arrival_dates(dVmin_I(2))-693960;
    min_dv_set1_vals(case_num,14) = min_dv_set1_vals(case_num,12) - dv_min_arrival;
    total_dV_return_test1(dVmin_I(1),dVmin_I(2))=10000;
    total_dV_return_neg_test1(dVmin_I(1),dVmin_I(2))=10000;
end
dVmin_I(1) = NaN;
dVmin_I(2) = NaN;
xlswrite('C:\Users\jfitc\OneDrive\NASA\JSC Spring 2019\IAC Paper Resources\Ceres Porkchop Data.xlsx',min_dv_set1_vals(:,:),'Sheet1','C13')

% Evaluate trips there before return trip
min_dv_set2_vals = zeros(4,14);
for case_num = 1:4
    min_mission_dv = 10000;
    for launch_index = 1:I_col_dv_b
       for arrival_index = 1:(I_col_dv_b - min_timeflight/date_step)
           test_mission_dv = dv_min_return_dv + total_dV_test2(launch_index,arrival_index);
           mission_time = (dv_min_return_arrival + 693960) - launch_dates(launch_index);
           ceres_stay = (dv_min_return_launch + 693960) - arrival_dates(arrival_index);
           if test_mission_dv < min_mission_dv && mission_time < max_allow_mission_time && (ceres_stay > min_desired_ceres_stay)
               min_mission_dv = test_mission_dv;
               dVmin_I = [launch_index,arrival_index];
               is_neg = 0;
           end
           test_mission_dv_neg = dv_min_return_dv + total_dV_neg_test2(launch_index,arrival_index);
           if test_mission_dv_neg < min_mission_dv && mission_time < max_allow_mission_time && (ceres_stay > min_desired_ceres_stay)
               min_mission_dv = test_mission_dv_neg;
               dVmin_I = [launch_index,arrival_index];
               is_neg = 1;
           end
       end
    end
    min_dv_set2_vals(case_num,1) = min_mission_dv;
    min_dv_set2_vals(case_num,3) = dv_min_return_dv;
    min_dv_set2_vals(case_num,6) = dv_min_return_c3;
    min_dv_set2_vals(case_num,7) = dv_min_return_vinf;
    if is_neg
        min_dv_set2_vals(case_num,2) = total_dV_neg_test2(dVmin_I(1),dVmin_I(2));
        min_dv_set2_vals(case_num,4) = c3_E_neg(dVmin_I(1),dVmin_I(2));
        min_dv_set2_vals(case_num,5) = v_inf_neg(dVmin_I(1),dVmin_I(2));
    else 
        min_dv_set2_vals(case_num,2) = total_dV_test2(dVmin_I(1),dVmin_I(2));
        min_dv_set2_vals(case_num,4) = c3_E(dVmin_I(1),dVmin_I(2));
        min_dv_set2_vals(case_num,5) = v_inf(dVmin_I(1),dVmin_I(2));
    end
    min_dv_set2_vals(case_num,8) = time_flight(dVmin_I(1),dVmin_I(2));
    min_dv_set2_vals(case_num,9) = dv_min_return_tof;
    min_dv_set2_vals(case_num,10) = launch_dates(dVmin_I(1))-693960;
    min_dv_set2_vals(case_num,11) = arrival_dates(dVmin_I(2))-693960;
    min_dv_set2_vals(case_num,12) = dv_min_return_launch;
    min_dv_set2_vals(case_num,13) = dv_min_return_arrival;
    min_dv_set2_vals(case_num,14) = dv_min_return_launch - min_dv_set2_vals(case_num,11);
    total_dV_test2(dVmin_I(1),dVmin_I(2))=10000;
    total_dV_neg_test2(dVmin_I(1),dVmin_I(2))=10000;
end
xlswrite('C:\Users\jfitc\OneDrive\NASA\JSC Spring 2019\IAC Paper Resources\Ceres Porkchop Data.xlsx',min_dv_set2_vals(:,:),'Sheet1','C17')
dVmin_I(1) = NaN;
dVmin_I(2) = NaN;

% Evaluate all possible trips
min_dv_set3_vals = zeros(10,14);
previous_min = 0;
for case_num = 1:10
    min_mission_dv = 100000;
    for launch_index = 1:(n_launch_dates - min_timeflight/date_step)
       for arrival_index = (launch_index + min_timeflight/date_step):n_arrival_dates
           for launch_back_index = arrival_index:n_launch_dates
               for arrival_back_index = (launch_back_index + min_timeflight/date_step):n_arrival_dates
                   test_mission_dv = total_dV_test3(launch_index,arrival_index) + total_dV_return_test3(launch_back_index,arrival_back_index);
                   test_mission_dv_neg = total_dV_neg_test3(launch_index,arrival_index) + total_dV_return_neg_test3(launch_back_index,arrival_back_index);
                   test_mission_dv_pos_neg = total_dV_test3(launch_index,arrival_index) + total_dV_return_neg_test3(launch_back_index,arrival_back_index);
                   test_mission_dv_neg_pos = total_dV_neg_test3(launch_index,arrival_index) + total_dV_return_test3(launch_back_index,arrival_back_index);
                   mission_time = return_arrival_dates(arrival_back_index) - launch_dates(launch_index);
                   ceres_stay = return_launch_dates(launch_back_index) - arrival_dates(arrival_index); 
                   if (test_mission_dv<min_mission_dv)&&(mission_time<max_allow_mission_time)&&(ceres_stay>min_desired_ceres_stay)&&(test_mission_dv>previous_min)
                       min_mission_dv = test_mission_dv;
                       dVmin_I = [launch_index,arrival_index];
                       dVmin_I2 = [launch_back_index,arrival_back_index];
                       pos_neg = 1;
                   end
                   if (test_mission_dv_neg<min_mission_dv)&&(mission_time<max_allow_mission_time)&&(ceres_stay>min_desired_ceres_stay)&&(test_mission_dv_neg>previous_min)
                       min_mission_dv = test_mission_dv_neg;
                       dVmin_I = [launch_index,arrival_index];
                       dVmin_I2 = [launch_back_index,arrival_back_index];
                       pos_neg = 2;
                   end
                   if (test_mission_dv_pos_neg<min_mission_dv)&&(mission_time<max_allow_mission_time)&&(ceres_stay>min_desired_ceres_stay)&&(test_mission_dv_pos_neg>previous_min)
                       min_mission_dv = test_mission_dv_pos_neg;
                       dVmin_I = [launch_index,arrival_index];
                       dVmin_I2 = [launch_back_index,arrival_back_index];
                       pos_neg = 3;
                   end
                   if (test_mission_dv_neg_pos<min_mission_dv)&&(mission_time<max_allow_mission_time)&&(ceres_stay>min_desired_ceres_stay)&&(test_mission_dv_neg_pos>previous_min)
                       min_mission_dv = test_mission_dv_neg_pos;
                       dVmin_I = [launch_index,arrival_index];
                       dVmin_I2 = [launch_back_index,arrival_back_index];
                       pos_neg = 4;
                   end
               end
           end
       end
    end
    min_dv_set3_vals(case_num,1) = min_mission_dv;
    switch pos_neg  % 1 = pos both ways; 2 = neg both ways; 3 = pos there neg back; 4 = neg there pos back
        case 1
            min_dv_set3_vals(case_num,2) = total_dV_test3(dVmin_I(1),dVmin_I(2));
            min_dv_set3_vals(case_num,3) = total_dV_return_test3(dVmin_I2(1),dVmin_I2(2));
            min_dv_set3_vals(case_num,4) = c3_E(dVmin_I(1),dVmin_I(2));
            min_dv_set3_vals(case_num,5) = v_inf(dVmin_I(1),dVmin_I(2));
            min_dv_set3_vals(case_num,6) = c3_E_return(dVmin_I2(1),dVmin_I2(2));
            min_dv_set3_vals(case_num,7) = v_inf_return(dVmin_I2(1),dVmin_I2(2));
        case 2
            min_dv_set3_vals(case_num,2) = total_dV_neg_test3(dVmin_I(1),dVmin_I(2));
            min_dv_set3_vals(case_num,3) = total_dV_return_neg_test3(dVmin_I2(1),dVmin_I2(2));
            min_dv_set3_vals(case_num,4) = c3_E_neg(dVmin_I(1),dVmin_I(2));
            min_dv_set3_vals(case_num,5) = v_inf_neg(dVmin_I(1),dVmin_I(2));
            min_dv_set3_vals(case_num,6) = c3_E_return_neg(dVmin_I2(1),dVmin_I2(2));
            min_dv_set3_vals(case_num,7) = v_inf_return_neg(dVmin_I2(1),dVmin_I2(2));
        case 3
            min_dv_set3_vals(case_num,2) = total_dV_test3(dVmin_I(1),dVmin_I(2));
            min_dv_set3_vals(case_num,3) = total_dV_return_neg_test3(dVmin_I2(1),dVmin_I2(2));
            min_dv_set3_vals(case_num,4) = c3_E(dVmin_I(1),dVmin_I(2));
            min_dv_set3_vals(case_num,5) = v_inf(dVmin_I(1),dVmin_I(2));
            min_dv_set3_vals(case_num,6) = c3_E_return_neg(dVmin_I2(1),dVmin_I2(2));
            min_dv_set3_vals(case_num,7) = v_inf_return_neg(dVmin_I2(1),dVmin_I2(2));
        case 4
            min_dv_set3_vals(case_num,2) = total_dV_neg_test3(dVmin_I(1),dVmin_I(2));
            min_dv_set3_vals(case_num,3) = total_dV_return_test3(dVmin_I2(1),dVmin_I2(2));
            min_dv_set3_vals(case_num,4) = c3_E_neg(dVmin_I(1),dVmin_I(2));
            min_dv_set3_vals(case_num,5) = v_inf_neg(dVmin_I(1),dVmin_I(2));
            min_dv_set3_vals(case_num,6) = c3_E_return(dVmin_I2(1),dVmin_I2(2));
            min_dv_set3_vals(case_num,7) = v_inf_return(dVmin_I2(1),dVmin_I2(2));            
    end
    min_dv_set3_vals(case_num,8) = time_flight(dVmin_I(1),dVmin_I(2));
    min_dv_set3_vals(case_num,9) = time_flight_back(dVmin_I2(1),dVmin_I2(2));
    min_dv_set3_vals(case_num,10) = launch_dates(dVmin_I(1))-693960;
    min_dv_set3_vals(case_num,11) = arrival_dates(dVmin_I(2))-693960;
    min_dv_set3_vals(case_num,12) = return_launch_dates(dVmin_I2(1))-693960;
    min_dv_set3_vals(case_num,13) = return_arrival_dates(dVmin_I2(2))-693960;
    min_dv_set3_vals(case_num,14) = min_dv_set3_vals(case_num,12) - min_dv_set3_vals(case_num,11);
    previous_min = min_mission_dv;
end
xlswrite('C:\Users\jfitc\OneDrive\NASA\JSC Spring 2019\IAC Paper Resources\Ceres Porkchop Data.xlsx',min_dv_set3_vals,'Sheet1','C22')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_verification = [c3_min_c3, c3_min_vinf, c3_min_tof;
    c3_min_return_c3, c3_min_return_vinf, c3_min_return_tof;
    vinf_min_vinf, vinf_min_c3, vinf_min_tof;
    vinf_min_return_vinf, vinf_min_return_c3, vinf_min_return_tof];
%disp(data_verification);

%}
