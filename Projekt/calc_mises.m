function mises = calc_mises(s, alpha, E, dT, ny)
% CALC_MISES Calculates the von mises stress for an element
% with stresses = [sigma_xx, sigma_yy, sigma_xy]

s_zz = ny * (s(1) + s(2)) - alpha*E*dT;

mises = sqrt(s(1)^2 + s(2)^2 + s_zz^2 - s(1)*s(2) -s(1)*s_zz - s(2)*s_zz + 3 * s(3)^2);

end

