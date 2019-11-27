function [dJ] = running_cost_minmax(t, J, x_pp, u_pp, v_pp, Q, Ru, Rv)
x = ppval(x_pp, t);
u = ppval(u_pp, t);
v = ppval(v_pp, t);

dJ = x' * Q * x + u' * Ru * u - v' * Rv * v;


end