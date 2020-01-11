function f = x1dot(x1, x2, x3, x4, t)
 global beta rho gamma;
 f = x4 * x2 - x3 * x4 + beta * (x2 - x4) + gamma * (x2^2 - x4 * x1);
end