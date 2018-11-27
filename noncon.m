function [c_n, ceq] = noncon(l,lk,rho)

for i = 1:9
    c(i) = norm(l([2*i-1,2*i])-lk([2*i-1,2*i]),2)-rho;
end

c1 = zeros(36,1);
count = 1;
D=100;
for i=1:8
   for j=i+1:9
        c1(count)=-(lk([2*i-1,2*i])-lk([2*j-1,2*j]))'*(l([2*i-1,2*i])-l([2*j-1,2*j]))...
        +5*D*norm(lk([2*i-1,2*i])-lk([2*j-1,2*j])) ;
        count = count + 1;

   end
end

c_n = [c';c1];
ceq = [];

end
