clear all
thetaw=3; alpha=1/3; R=50; k=.005; U=10; A=pi*R^2; rho=1.225; Cp=1.34; % Place Holders
xi=[0 5 10 15 20 0 5 10 15 20];
yi=[0 0 0 0 0 5 5 5 5 5];

xj=[0 5 10 15 20 0 5 10 15 20];
yj=[0 0 0 0 0 5 5 5 5 5];



% for i=1:25 
% l(2*i-1)=x(i);
% l(2*i)=y(i);
% end

% m=1;
% for i=1:25
%     for j=1:25
% 
% theta_ij(i)=sec(dot(l([2*i-1,2*i]),l([2*j-1,2*j])))/(norm(l([2*i-1,2*i])*norm(l([2*j-1,2*j]))));
% r_ij(i)=norm(l([2*i-1,2*i])-l([2*j-1,2*j]))*sin(abs(theta_ij(i)-thetaw));
% d_ij(i)=norm(l([2*i-1,2*i])-l([2*j-1,2*j]))*cos(abs(theta_ij(i)-thetaw));
% du_ij(i)= 2*alpha*(R/(R+k*d_ij(i)))^2*exp(-(r_ij(i)/(R+k*d_ij(i)))^2);
% dubar_ij=sqrt(sum(du_ij));
% ubar_ij(i)=U*(1-dubar_ij);
% P1(i)=0.5*rho*A*Cp*(ubar_ij(i))^3;

P2=0.5*rho*A*Cp*(U*(1-(sqrt(sum(2*alpha*(R/(R+k*(norm(l([2*i-1,2*i])-l([2*j-1,2*j]))...
    *cos(abs((sec(dot(l([2*i-1,2*i]),l([2*j-1,2*j])))/(norm(l([2*i-1,2*i])*norm(l([2*j-1,2*j])))))-thetaw)))))^2....
*exp(-((norm(l([2*i-1,2*i])-l([2*j-1,2*j]))*sin(abs((sec(dot(l([2*i-1,2*i]),l([2*j-1,2*j])))/(norm(l([2*i-1,2*i])...
*norm(l([2*j-1,2*j])))))-thetaw)))/(R+k*(norm(l([2*i-1,2*i])-l([2*j-1,2*j]))*cos(abs((sec(dot(l([2*i-1,2*i]),l([2*j-1,2*j])))/...
(norm(l([2*i-1,2*i])*norm(l([2*j-1,2*j])))))-thetaw)))))^2))))))^3;

%     end
% end
 
for i = 1:10
theta_ij(i) = atand(abs((yj(i)-yi(i)) / abs(xj(i)-xi(i))));
d_ij(i) = sqrt(((xi(i) - xj(i))^2) + (yi(i) - yj(i)))^2 * cos(abs(theta_ij-thetaw));
r_ij(i) = sqrt(((xi(i) - xj(i))^2) + (yi(i) - yj(i)))^2 * sin(abs(theta_ij-thetaw));

end
