function y = convexsub(l, g, Bk, lk)
% x: the variable
% g: current function gradient (column)
% B: BFGS app. hessian of the Lagrangian
% xk: current value of the variables

y = g'*(l - lk) + (l - lk)'*Bk*(l - lk);
end 
