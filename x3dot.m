function f = x3dot(x1, x2, x3, x4, t)
 global beta rho gamma;
 f = x2 * x4 - rho * x1 * x2 + beta * (x4 - x2) + gamma * (x4^2 - x2 * x3);
end