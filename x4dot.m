function f = x4dot(x1, x2, x3, x4, t)
 global beta rho gamma;
 f = x3 * x1 - x2 * x3 + beta * (x1 - x3) + gamma * (x1^2 - x3 * x4);
end