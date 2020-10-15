function D = calc_dstrain(E, ny)
%CALC_DSTRAIN  Calculated D-strain matrix

D = E/((1+ny)*(1-2*ny)) * [1 - ny, ny, 0;
                           ny, 1 - ny, 0;
                           0, 0, 0.5*(1-2*ny)];

end

