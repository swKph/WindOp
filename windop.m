clear all
thetaw=0; alpha=1/3; R=50; k=2; U=8.9408; A=pi*R^2; rhoa=1.225; Cp=1; D=2*R; % Place Holders
syms l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18
l_sym=[l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18];

for i=1:9
    for j=1:9
theta_ij(i)=sec(dot(l_sym([2*i-1,2*i]),l_sym([2*j-1,2*j])))/(norm(l_sym([2*i-1,2*i])*norm(l_sym([2*j-1,2*j]))));
r_ij(i)=norm(l_sym([2*i-1,2*i])-l_sym([2*j-1,2*j]))*sin(abs(theta_ij(i)-thetaw));
d_ij(i)=norm(l_sym([2*i-1,2*i])-l_sym([2*j-1,2*j]))*cos(abs(theta_ij(i)-thetaw));
du_ij(i)= 2*alpha*(R/(R+k*d_ij(i)))^2*exp(-(r_ij(i)/(R+k*d_ij(i)))^2);
dubar_ij(i)=sqrt(sum(du_ij));
ubar_ij(i)=U*(1-dubar_ij(i));
P=0.5*rhoa*A*Cp*(ubar_ij(i))^3; %Solving for the Power equation in terms of l for every l value
    end
end

% l1=lk(1); l3=lk(3); l5=lk(5); l7=lk(7); l9=lk(9); l11=lk(11); l13=lk(13); l15=lk(15); l17=lk(17);
% l2=lk(2); l4=lk(4); l6=lk(6); l8=lk(8); l10=lk(10); l12=lk(12); l14=lk(14); l16=lk(16); l18=lk(18);

lk=[1;1;6;1;11;1;1;6;6;6;11;6;1;11;6;11;11;11];
     %initial point
     
n = norm((lk-1),2);
while n > 0.001

 % initial guess is lk
l0 = lk;

A = constraint1(l, lk);
b = zeros(1,length(A));
Aeq = []; beq = [];
lb=l0-2.5*D;
g_sym=symbolicgradient(l);
g = cal_g(g_sym, lk);
ub=l0+2.5*D;
[l, fval, exitflag, ~, lambda] = fmincon(@(l)convexsub(l,g,Bk,lk),l0,A,b,...
    Aeq,beq,lb,ub,@(l)noncon(l,rho)); %insert fmin con function

% do trust region
alpha = (P(lk)-P(l))/(convexsub(l,g,Bk,lk)-P(l));
    if alpha >= 0.2
    lk = l;
    rho = 1.1* rho;
    Bk=laghessian(l, cal_g);
    
    else
    rho = 0.5* rho;
    end

end

function y = convexsub(l, g, Bk, lk)
% x: the variable
% g: current function gradient (column)
% B: BFGS app. hessian of the Lagrangian
% xk: current value of the variables

y = g'*(l - lk) + (l - lk)'*Bk*(l - lk);
end 

function g_sym= symbolicgradient(l)
thetaw=0; alpha=1/3; R=50; k=2; U=8.9408; A=pi*R^2; rho=1.225; Cp=1; % Place Holders
syms l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18
l=[l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18];
for i=1:9
    for j=1:9
theta_ij(i)=sec(dot(l([2*i-1,2*i]),l([2*j-1,2*j])))/((norm(l([2*i-1,2*i])))*norm(((l([2*j-1,2*j])))));
r_ij(i)=norm((((l([2*i-1,2*i]))-l([2*j-1,2*j]))))*sin(abs(theta_ij(i)-thetaw));
d_ij(i)=norm(((l([2*i-1,2*i])-l([2*j-1,2*j]))))*cos(abs(theta_ij(i)-thetaw));
du_ij(i)= 2*alpha*(R/(R+k*d_ij(i)))^2*exp(-(r_ij(i)/(R+k*d_ij(i)))^2);
dubar_ij(i)=sqrt(sum(du_ij));
ubar_ij(i)=U*(1-dubar_ij(i));
P=0.5*rho*A*Cp*(ubar_ij(i))^3; %Solving for the Power equation in terms of l for every l value
    end
end
g_sym = gradient(P);

end

function g = cal_g(g_sym, lk)

l1=lk(1); l3=lk(3); l5=lk(5); l7=lk(7); l9=lk(9); l11=lk(11); l13=lk(13); l15=lk(15); l17=lk(17);
l2=lk(2); l4=lk(4); l6=lk(6); l8=lk(8); l10=lk(10); l12=lk(12); l14=lk(14); l16=lk(16); l18=lk(18);
    
g = eval(g_sym);
end


function Bk=laghessian(l ,lk, cal_g) % damped BFGS update
B=eye(18);
s=lk-l;

new_g_L = cal_g(g_sym, lk) + lambda*gradient_constraints(lk, l);
old_g_L = cal_g(g_sym, l) + lambda*gradient_constraints(l, lk);


y=new_g_L - old_g_L;

if s'*y >= 0.2*s'*B*s 
    theta = 1; 
else
    theta = (0.8*s'*B*s)/(s'*B*s-s'*y);
end

ytil=theta*y+(1-theta)*B*s;

Bk=B+(1/(ytil'*s))*(ytil'*ytil)-(1/(s'*B*s))*(B*s'*s*B);
end


function T= noncon(lk,rho) %trust region constraint

T=norm((l-lk),2)-rho;

end


function GC=gradient_constraints(l, lk) %gradient of the constraints symbolically
D=100;
for i=1:9
    for j=1:9
GC=[lk([2*j-1,2*j])-lk([2*i-1,2*i]), -lk([2*j-1,2*j])+lk([2*i-1,2*i]),...
    l([2*j-1,2*j]) - l([2*i-1,2*i]) + (500*sign(lk([2*i-1,2*i]) - lk([2*j-1,2*j]))...
    *abs(lk([2*i-1,2*i]) - lk([2*j-1,2*j])))/(abs(lk([2*i-1,2*i]) - lk([2*j-1,2*j]))^2)^(1/2),...
    l([2*i-1,2*i]) - l([2*j-1,2*j]) - (500*sign(lk([2*i-1,2*i]) - lk([2*j-1,2*j]))...
    *abs(lk([2*i-1,2*i]) - lk([2*j-1,2*j])))/(abs(lk([2*i-1,2*i]) - lk([2*j-1,2*j]))^2)^(1/2);...
    (2*lk(2*i-1) - 2*l(2*i-1))/(2*((lk(2*i-1) - l(2*i-1))^2 + (lk(2*i) - l(2*i))^2)^(1/2)),...
    (2*lk(2*i) - 2*l(2*i))/(2*((lk(2*i-1) - l(2*i-1))^2 + (lk(2*i) - l(2*i))^2)^(1/2)),...
    -(2*lk(2*i-1) - 2*l(2*i-1))/(2*((lk(2*i-1) - l(2*i-1))^2 + (lk(2*i) - l(2*i))^2)^(1/2)),...
    -(2*lk(2*i) - 2*l(2*i))/(2*((lk(2*i-1) - l(2*i-1))^2 + (lk(2*i) - l(2*i))^2)^(1/2))];
    end 
end 
end


function c1=constraint1(l,lk) %constraint function
D=100;
for i=1:9
   for j=1:9
c1(i)=-(lk([2*i-1,2*i])-lk([2*j-1,2*j]))'*(l([2*i-1,2*i])-l([2*j-1,2*j]))+5*D*norm(lk([2*i-1,2*i])-lk([2*j-1,2*j])) ;
c1(i+9)=-(lk([2*i-1,2*i])-lk([2*j-1,2*j]))'*(l([2*i-1,2*i])-l([2*j-1,2*j]))+5*D*norm(lk([2*i-1,2*i])-lk([2*j-1,2*j])) ;

   end
end
end
