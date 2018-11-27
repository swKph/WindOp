function Bk=laghessian(lambda, l ,lk, g, ncLk)
B=eye(18);
s=lk-l;

new_g_L = g + (lambda'*ncLk)';
old_g_L = g + (lambda'*ncLk)';


y=new_g_L - old_g_L;

if s'*y >= 0.2*s'*B*s 
    theta = 1; 
else
    theta = (0.8*s'*B*s)/(s'*B*s-s'*y);
end

ytil=theta*y+(1-theta)*B*s;

Bk=B+(1/(ytil'*s))*(ytil*ytil')-(1/(s'*B*s))*(B*s*s'*B);
end