function result = fading2D(conc, fading, threshold)

result = conc(:, :) * (1 - fading);
result((conc(:, :) < threshold & conc(:, :) > 0)) = 0;


end