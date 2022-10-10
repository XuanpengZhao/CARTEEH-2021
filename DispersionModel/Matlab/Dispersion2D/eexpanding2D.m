function result = eexpanding2D(conc, gridSize, expanding, timeStep, wind)

 
%coef = ((expanding + 1) ^ timeStep - 1) / (timeStep * expanding);%0.0353;
%(sqrt(timeStep * expanding ^ 3 + timeStep * expanding^2) - expanding) / (timeStep * expanding ^ 2)
% ((expanding + 1) ^ timeStep - 1) / (timeStep * expanding);
L_e = (timeStep * expanding) * gridSize + gridSize;

% if L_e > 3
%     L_e = 3;
% end
L_e_y_plus = ((timeStep * expanding) / 2) + wind(2) * timeStep;
L_e_x_plus = ((timeStep * expanding) / 2) + wind(1) * timeStep;
L_e_y_minus = ((timeStep * expanding) / 2) - wind(2) * timeStep;
L_e_x_minus = ((timeStep * expanding) / 2) - wind(1) * timeStep;
if (gridSize - ((timeStep * expanding) / 2) * gridSize < abs(wind(2) * timeStep)) || ...
       (gridSize - ((timeStep * expanding) / 2) * gridSize < (abs(wind(1) * timeStep)))
   TimeStepORWindSpeedTooLarge__ = 123456
end



% topLeft = max(min(L_e_y_plus, 1), 0) * max(min(L_e_x_plus, 1), 0);
% top = max(min(L_e_y_plus, 1), 0) * max(min(L_e_x_minus + 1, 1), 0);
% topRight =  max(min(L_e_y_plus, 1), 0) * max(min(L_e_x_minus, 1), 0);
% midLeft = max(min(L_e_y_minus + 1, 1), 0) * max(min(L_e_x_plus, 1), 0);
% midRight = max(min(L_e_y_minus + 1, 1), 0) * max(min(L_e_x_minus, 1), 0);
% botLeft =  max(min(L_e_y_minus, 1), 0) * max(min(L_e_x_plus, 1), 0);
% bot = max(min(L_e_y_minus, 1), 0) * max(min(L_e_x_minus + 1, 1), 0);
% botRight = max(min(L_e_y_minus, 1), 0) * max(min(L_e_x_minus, 1), 0);
% 
topLeft = max(min(L_e_y_plus, gridSize), 0) * max(min(L_e_x_minus, gridSize), 0);
top = max(min(L_e_y_plus, gridSize), 0) * max(min(L_e_x_minus + gridSize, gridSize), 0);
topRight =  max(min(L_e_y_plus, gridSize), 0) * max(min(L_e_x_plus, gridSize), 0);
midLeft = max(min(L_e_y_minus + gridSize, gridSize), 0) * max(min(L_e_x_minus, gridSize), 0);
midRight = max(min(L_e_y_minus + gridSize, gridSize), 0) * max(min(L_e_x_plus, gridSize), 0);
botLeft =  max(min(L_e_y_minus, gridSize), 0) * max(min(L_e_x_minus, gridSize), 0);
bot = max(min(L_e_y_minus, gridSize), 0) * max(min(L_e_x_minus + gridSize, gridSize), 0);
botRight = max(min(L_e_y_minus, gridSize), 0) * max(min(L_e_x_plus, gridSize), 0);

mid = 0;
result = [topLeft, top, topRight, midLeft, mid, midRight, botLeft, bot, botRight];
result = result / L_e^2;
result(5) = 1 - sum(result);
result = result * conc;

% reshape(result, 3, 3)'
end