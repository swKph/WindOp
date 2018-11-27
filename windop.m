%Design Optimization of Wind Turbine Placement using SQP and BFGS with a
%Trust Region

clear all

thetaw=pi; alpha=1/3; R=50; k=.1; U=12; A=pi*R^2; rhoa=1.225; Cp=1; D=2*R; %Initial conditions for theorized wind farm location

syms l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18 %symbolic placeholders for calculating the power function symbolically
l_sym=[l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18];

%power function loop to calclate power equation
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

lk=[1.1;1.2;6.1;1.4;11.1;1.5;1.1;6.2;6.2;6.3;11.6;6.1;1.3;11.1;6.2;11.7;11.2;11.5;]; %initial point
lk=R.*lk; %original turbine placement location non ideal - 1.5 radii times the grid
l_original = lk; 
l_old = 1e1*ones(size(lk));
rho=1;
Bk=eye(18);

while norm([lk-l_old],2)>0.001

% initial guess is lk
l0 = lk;
A = [];

% constraint1(l, lk);
b = [];

Aeq = []; beq = [];
lb=l0-2.5*D;
ub=l0+2.5*D;

g_sym=symbolicgradient(l_old);
g = cal_g(g_sym,lk);

[l, fval, exitflag, ~, lambda] = fmincon(@(l)convexsub(l,g,Bk,lk),l0,A,b,...
    Aeq,beq,lb,ub,@(l)noncon(l,lk,rho)); %insert fmin con function
nc=sym_nc(l, lk, rho);
ncLk=cal_nc(nc, l ,lk);
lambda = lambda.ineqnonlin;
pLk = cal_p(P,lk);
pL = cal_pL(P,l);

% do trust region
alpha = (pLk-pL)/(convexsub(l_old,g,Bk,lk)-pL);
if alpha <= 0.2
    
    lk = l;
    l_old=lk;
    
    rho = 1.1* rho;
    Bk=laghessian(lambda,l,lk, g, ncLk);
    display(alpha);
    
else
    l_old = l;
    rho = 0.5* rho;
    display(alpha);
end

end

%parse x y data

for i = 1:9
    
    x(i) = lk(2 * i - 1);
    y(i) = lk(2 * i);
    
    xx(i) = l_original(2 * i - 1);
    yy(i) = l_original(2 * i);
    
end

table1 = table(l_original,lk);

figure
hold on
scatter(x,y,'b','marker','x')
scatter(xx,yy,'r','marker','o')
line([x; xx], [y; yy],'Color', 'k', 'LineStyle', '--');
title('Optimized Turbine Locations - Wind Direction - 90* Wind Speed - 12 m/s')
xlabel('x (m)')
ylabel('y (m)')
legend('Optimized Locations','Original Locations','Location','bestoutside')
grid on
grid minor
hold off
