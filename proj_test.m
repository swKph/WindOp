li = [0 5 10 15 20 0 5 10 15 20 0 5 10 15 20 0 5 10 15 20 0 5 10 15 20];
li = li';
lj = [0 0 0 0 0 5 5 5 5 5 10 10 10 10 10 15 15 15 15 15 20 20 20 20 20];
lj = lj';

l = [li,lj];

theta_ij = zeros(25,1);

for i = 1:24
theta_ij(i) = atan((l(i,1) - l(i+1,1))/(l(i,2) - l(i+1,2)));
end

% theta_w = pi/2;
% U = 10; %m/s
% r_ij = norm(li - lj) * sin(abs(theta_ij - theta_w));
% d_ij = norm(li - lj) * cos(abs(theta_ij - theta_w));
% 
% alpha = 1/3; %ideal
% R = 100; %meters
% K = (.5 / log(80/.01));
% 
% du = 2*alpha* ( R / (R + K*d_ij))^2 * exp(-(r_ij/(R+K*d_ij))^2);
% 
% du_bar = sum(sqrt(du));
% 
% u_bar = U*(1-du_bar);
% 
% rho = 1.225; %kg/m^3
% A = pi*R^2;
% Cp = 1; %kj/kg
% P = (.5 * rho * A * Cp * u_bar^3);
% 
% f = sum(P);
% 
% li = [0 5 10 15 20 0 5 10 15 20 0 5 10 15 20 0 5 10 15 20 0 5 10 15 20];
% li = li';
% lj = [0 0 0 0 0 5 5 5 5 5 10 10 10 10 10 15 15 15 15 15 20 20 20 20 20];
% lj = lj';
% 
% l = [li,lj];
% 
% 
