function g_sym= symbolicgradient(l_old)

thetaw=90; alpha=1/3; R=50; k=2; U=30; A=pi*R^2; rhoa=1.225; Cp=1; % Place Holders
syms l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18
l_old=[l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18];

for i=1:9
    for j=1:9
        theta_ij(i)=sec(dot(l_old([2*i-1,2*i]),l_old([2*j-1,2*j])))/((norm(l_old([2*i-1,2*i])))*norm(((l_old([2*j-1,2*j])))));
        r_ij(i)=norm((((l_old([2*i-1,2*i]))-l_old([2*j-1,2*j]))))*sin(abs(theta_ij(i)-thetaw));
        d_ij(i)=norm(((l_old([2*i-1,2*i])-l_old([2*j-1,2*j]))))*cos(abs(theta_ij(i)-thetaw));
        du_ij(i)= 2*alpha*(R/(R+k*d_ij(i)))^2*exp(-(r_ij(i)/(R+k*d_ij(i)))^2);
        dubar_ij(i)=sqrt(sum(du_ij));
        ubar_ij(i)=U*(1-dubar_ij(i));
        P=0.5*rhoa*A*Cp*(ubar_ij(i))^3; %Solving for the Power equation in terms of l for every l value
    end
end

g_sym = gradient(P);

end