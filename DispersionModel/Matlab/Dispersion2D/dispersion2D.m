close all;clear;clc;
WIDTH = 5;
FONTSIZE = 24;
%% Open & Reading
fid = fopen("../FeildTest/20211109/Vehicle_CO2_1109_JL.csv" );
data = textscan(fid, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','Headerlines',1,'delimiter', ',');
fclose(fid);
WINDSPEED = data{1, 5};
WINDDIREC = data{1, 6};
MEASX = data{1, 19};
MEASY = data{1, 20};
MEASZ = data{1, 4};
MEASCONC = data{1, 13};
MEASFR = data{1, 14};
MEASEMIS = data{1, 22};
%%
timeGap = 0.1; % s
gridSize = round(2*timeGap*10/(2-timeGap*1.1));  % meter
X = 100 * gridSize; % meter
Y = 100 * gridSize; % meter



 
 
% Lane
lane = [0 / gridSize, Y/2 / gridSize;
        X / gridSize, Y/2 / gridSize];
laneWidth = 1;% meter

% Vehicles
vehLength = 2 / gridSize;
vehWidth = 3 / gridSize;
vehicle = [50 * gridSize, 50 * gridSize];
% 10, (X / 2 - vehWidth / 2 + 9) / gridSize - 1;
% 10, (X / 2 - vehWidth / 2 + 3) / gridSize - 1;
% 10, (X / 2 - vehWidth / 2 - 3) / gridSize - 1;
% 10, (X / 2 - vehWidth / 2 - 9) / gridSize - 1;
% 10, (X / 2 - vehWidth / 2 - 15) / gridSize - 1];
%            Y / 2 / gridSize, (X / 2 ) / gridSize ];
       %Y / 2 / gridSize, (X / 2 - vehWidth / 2) / gridSize
vehTurbPadding = 1;
tailpip = [ vehLength / 2,  vehWidth];
%              vehLength / 2 + 1, - vehWidth / 2;];


speed = [0, 0]/gridSize;
% 5, 0;
% 5, 0;
% 5, 0;
% 5, 0;
% 5, 0;];
%         -5, 0;] / gridSize % m/s


maxTimeStep = 1000 / timeGap;
emissionRate = 214; % ppm/s Fuel rate -> 8887 = g/hr -> ppm * 41* Molar mass 
emissionSpeed = 50;
% function
expanding = 0.5;
fading = 0.01;
fadingThreshold = 1.0e-10;
steadyThreshold = 10e-5;
minError = 10^10;
%%
concSameHeight2 = [];
k = 1;
for X = 1:1:3
    for Y = -2:1:2
        j = 1;
        temp = [];
        for Z = [0.3, 0.6, 1.2]
            for i = 1:length(MEASZ)
                if MEASZ(i) == Z && MEASX(i) == X && MEASY(i) == Y && MEASFR(i) < 0.2
                    temp(j) = MEASCONC(i);
                    j = j+1;
                end
            end
        end
         
        if ~isempty(temp) && length(temp) == 3 && temp(1) > temp(2)*4 && temp(2)*4 > temp(3)*3
            temp(2) = temp(2) * 4;
            temp(3) = temp(3) * 3;
            concSameHeight(k, :) = temp / temp(1);
            concSameHeight2 = [concSameHeight2, concSameHeight(k, :)];
            k = k+1;
            [X, Y]
        end
    end
end
Z = [0.3, 0.6, 1.2];
Z4 = [0.3, 0.6, 1.2, 0.3, 0.6, 1.2, 0.3, 0.6, 1.2, 0.3, 0.6, 1.2, 0.3, 0.6, 1.2];

Z2 = [0.3, 0.6, 1.2, 0.3, 0.6, 1.2, 0.3, 0.6, 1.2, 0.3, 0.6, 1.2, 0.3, 0.6, 1.2, 0.3, 0.6, 1.2];
Z3 = [0.3, 0.6, 1.2, 0.3, 0.6, 1.2, 0.3, 0.6, 1.2, 0.3, 0.6, 1.2];
p = 6;
% exp(-(x-0.3)^2/(2*d^2))
% abs(d^(x-0.3))
% cos(d*(x-0.3))
% 2 / (1+exp(d*(x-0.3)))
myExp = fittype(' exp(-(x-0.3)^2/(2*d^2))');
% [cf, gof] = fit(Z', concSameHeight(p,:)', myExp);
[cf, gof] = fit(Z4(1:end)', concSameHeight2(1:end)', myExp);
xi = 0:0.1:10;
cf.d
yi = exp(-(xi-0.3).^2/(2*0.4845^2));%cf.d
%  yi = abs(cf.d).^(xi-0.3);
%  yi = cos(cf.d*(xi-0.3));
% yi = 2 ./ (1+exp(cf.d.*(xi-0.3)));
plot(xi, yi*100)
xlabel("Height (m)", "fontsize", 25)
ylabel("Percentage (%)", "fontsize", 25)
set(gca, "fontsize", 25)
hold on;
plot(Z4(1:end), concSameHeight2(1:end)*100, '*')
% plot(Z, concSameHeight(p,:), '*')
concSameHeight(2,:)
% cd.d = 0.7488, 0.3038
%% Calibration factors
% start from 0.21, 0.002, 0.08
% 0.3best: 0.09, 0.002, 0.15 : 1.0k 
% 0.3best: 0.1, 0.0015, 0.16 : 1.13k 
% 1.2best: 1, 0.009, 0.55 : 383
% 1.2best: 1, 0.003, 0.75 : 419 
for expanding =  3%1.5 %0.75%0.15:0.01:0.23%1.1 %1.6%1.3:0.05:1.7
    for fading = 0.21%0.41%0.003%0.001:0.001:0.003
        for coef = 1%0.05:0.01:0.1%8.5:0.1:9.5
            error = 0;
%% For measured points
height = 0.3;
sigma = 0.3038;
for meas = 1:49 %length(MEASZ)
%     if meas == 5% && meas ~= 12
%         continue;
%     end
    if MEASZ(meas) ~= height
        continue;
%         break;
    end
    
    concCurr = zeros(X / gridSize, Y / gridSize);
    concLast = concCurr;
    
    emissionRate = MEASEMIS(meas)*coef * abs(0.0083)^(height-0.3);%exp(-(height-0.3).^2/(2*sigma^2));
    calibration = [round(vehicle(1, 2) / gridSize + tailpip(2)) + MEASX(meas), round(vehicle(1, 1) / gridSize + tailpip(1)) + MEASY(meas)];
    % Wind
    windSpeed = WINDSPEED(meas);% / timeGap; % m/s
    windDirec = -(WINDDIREC(meas) / 360 * 2 * pi + pi/2) - pi/2; %1/2 * pi;
    % Wind X:1, Y:2
    windEnv(:, :, 1) = ones(size(concCurr, 1), size(concCurr, 2)) * windSpeed * cos(windDirec);
    windEnv(:, :, 2) = ones(size(concCurr, 1), size(concCurr, 2)) * windSpeed * sin(windDirec);
%     windEnv(53, 51, 1) = 0.05;
%     windEnv(53, 51, 2) = -0.05;
    windCurr = windEnv;
    windLast = windCurr;
    
%% For converged dispersion
elapsedTimeTotal = 0;
for i = 1:maxTimeStep
    
    
    tic 
    % Initial
    concCurr = zeros(X / gridSize, Y / gridSize);
    windCurr = zeros(X / gridSize, Y / gridSize, 2);
    windLogic = false(X / gridSize, Y / gridSize);
    % Expanding
    [row,col]  = find(concLast(:, :));
    if ~isempty(row)
        for j = 1:length(row)                                                              
            neighbors = eexpanding2D(concLast(row(j), col(j)), gridSize, expanding, timeGap, windLast(row(j), col(j), :));
            concCurr(row(j), col(j)) = concCurr(row(j), col(j)) + neighbors(5);
            if (row(j) - 1 > 0)
                concCurr(row(j) - 1, col(j)) = concCurr(row(j) - 1, col(j)) + neighbors(2);
                if (col(j) + 1 < Y / gridSize)
                    concCurr(row(j) - 1, col(j) + 1) = concCurr(row(j) - 1, col(j) + 1) + neighbors(3);
                end
                if (col(j) - 1 > 0)
                    concCurr(row(j) - 1, col(j) - 1) = concCurr(row(j) - 1, col(j) - 1) + neighbors(1);
                end
            end 
            if (row(j) + 1 < X / gridSize)
                concCurr(row(j) + 1, col(j)) = concCurr(row(j) + 1, col(j)) + neighbors(8);
                if (col(j) + 1 < Y / gridSize)
                    concCurr(row(j) + 1, col(j) + 1) = concCurr(row(j) + 1, col(j) + 1) + neighbors(9);
                end
                if (col(j) - 1 > 0)
                    concCurr(row(j) + 1, col(j) - 1) = concCurr(row(j) + 1, col(j) - 1) + neighbors(7);
                end
            end
            if (col(j) + 1 < Y / gridSize)
                concCurr(row(j), col(j) + 1) = concCurr(row(j), col(j) + 1) + neighbors(6);
            end
            if (col(j) - 1 > 0)
                concCurr(row(j), col(j) - 1) = concCurr(row(j), col(j) - 1) + neighbors(4);
            end
        end
    end
    
    % New source
    for j = 1:size(vehicle, 1)
        vehicle(j, :) = vehicle(j, :) + speed(j, :) * timeGap;
        
        concCurr(round(vehicle(j, 2) / gridSize + tailpip(2)), round(vehicle(j, 1) / gridSize + tailpip(1)))...
            = concCurr(round(vehicle(j, 2) / gridSize + tailpip(2)), round(vehicle(j, 1) / gridSize + tailpip(1))) + emissionRate * timeGap;
  
       
%         % Wind updates (surrounding vehicles)  
%         for k = 1:2
%             windCurr(round(vehicle(j, 2) - vehWidth / 2 - vehTurbPadding) : round(vehicle(j, 2) + vehWidth / 2 + vehTurbPadding),...
%             round(vehicle(j, 1) - vehLength / 2 - vehTurbPadding) : round(vehicle(j, 1) + vehLength / 2 + vehTurbPadding), k)...
%             = windCurr(round(vehicle(j, 2) - vehWidth / 2 - vehTurbPadding) : round(vehicle(j, 2) + vehWidth / 2 + vehTurbPadding),...
%             round(vehicle(j, 1) - vehLength / 2 - vehTurbPadding) : round(vehicle(j, 1) + vehLength / 2 + vehTurbPadding), k)...
%             + speed(j, k);
%             
%             windLogic(round(vehicle(j, 2) - vehWidth / 2 - vehTurbPadding) : round(vehicle(j, 2) + vehWidth / 2 + vehTurbPadding),...
%             round(vehicle(j, 1) - vehLength / 2 - vehTurbPadding) : round(vehicle(j, 1) + vehLength / 2 + vehTurbPadding), k) = 1;
%         end
    end
    % Concentration Fading
    concCurr(:, :) = fading2D(concCurr(:, :), fading, fadingThreshold);
    % Stop condition
    if sum(sum(concCurr - concLast)) < steadyThreshold
        break;
    end
    concLast = concCurr;
    windLast((windLogic)) = windCurr((windLogic));
    % Wind fading (Wind tends to enviroment)
     windLast = (windEnv - windLast) * timeGap + windLast;
%     elapsedTime = toc
%     elapsedTimeTotal = elapsedTimeTotal + elapsedTime
    
%%
% %     Draw concentration colormap
%     figure(1)
%     set(gcf, 'position',[100,0,2400,1400] )
%     clf
%     hold on
%     surf(concCurr(:, :))
%     colorbar
%     view(2)
%     view(0,-90)
%     axis equal
% %     concCurr(:,:)
% 
%     % vehicle & road display
%     for j = 1:size(vehicle, 1)
%         rectangle('Position',[round(vehicle(j, 1) ), round(vehicle(j, 2) ), vehLength / gridSize, vehWidth / gridSize], 'EdgeColor', 'red', 'LineWidth', 1)
%     end   
% %     laneLine = line(lane(:,1), lane( :,2), 'Color', 'white', 'LineWidth', 20);
%     hold off
%     
%     % Display wind 
% %     tic 
%     figure(2)
%     quiver(windLast(:, :, 1), windLast(:, :, 2))
%     view(2)
%     view(0,-90)
%     % vehicle & road display
%     for j = 1:size(vehicle, 1)
%         rectangle('Position',[round(vehicle(j, 1) ), round(vehicle(j, 2) ), vehLength / gridSize, vehWidth / gridSize], 'EdgeColor', 'red', 'LineWidth', 1)
%     end 
%   
% %     toc
% %     pause(0.01 - elapsedTime)
% %     clc;
    
end
    concCurr = concCurr ;%* 7.9;
    errorInst = abs(MEASCONC(meas) - concCurr(calibration(1), calibration(2)));
%     concCurr(calibration(1), calibration(2))
    error = error + errorInst;
    comparedMEAS2ESTI(meas, :) = [MEASCONC(meas), concCurr(calibration(1), calibration(2))];
    plot(comparedMEAS2ESTI(meas, 1), comparedMEAS2ESTI(meas, 2), "*", 'LineWidth', 10)
    
    hold on;
    
end
    xlabel("Simulation results", "fontsize", 25)
    ylabel("Measurement results", "fontsize", 25)
    set(gca, "fontsize", 25)
   plot([1, 2000], [1, 2000], "--", 'LineWidth', 1)
   figure
%    qqplot(comparedMEAS2ESTI(:, 1), comparedMEAS2ESTI(:, 2))
   if minError > error
       minError = error;
       calibratedExpanding = expanding;
       calibratedFading = fading;
       calibratedCoef = coef;%comparedMEAS2ESTI(:, 2)\ comparedMEAS2ESTI(:, 1);
   end
        end
             end
end