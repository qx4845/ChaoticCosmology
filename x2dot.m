function f = x2dot(x1, x2, x3, x4, t)
 global beta rho gamma;
 f = rho * x1 * x3 - x4 * x1 + beta * (x3 - x1) + gamma * (x3^2 - x1 * x2);
end