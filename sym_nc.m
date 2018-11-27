function nc = sym_nc(lk,l,rho)

% rho = 1.1;

syms l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18 lk1 lk2 lk3 lk4 lk5 lk6 lk7 lk8 lk9 lk10 lk11 lk12 lk13 lk14 lk15 lk16 lk17 lk18
l=[l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11 l12 l13 l14 l15 l16 l17 l18];

lk=[lk1 lk2 lk3 lk4 lk5 lk6 lk7 lk8 lk9 lk10 lk11 lk12 lk13 lk14 lk15 lk16 lk17 lk18];
 l=l';
 lk=lk';

nc_sym= cell(9,1);

for i = 1:9
    c(i) = norm(l(2*i-1:2*i,1)-lk(2*i-1:2*i,1),2)-rho;
%     c(i+9)=norm(l([2*i-1,2*i])-lk([2*i-1,2*i]),2)-rho;
    nc_sym{i,1} = gradient(c(i),l');
end

% c = norm([l-lk],2)-rho;
% gcc = gradient(c,l);
D=100;
nc1_sym=cell(9,9);
  for i=1:8
   for j=(i+1):9
        c1(i,j)=-(lk([2*i-1,2*i])-lk([2*j-1,2*j]))'*(l([2*i-1,2*i])-l([2*j-1,2*j]))...
        +5*D*norm(lk([2*i-1,2*i])-lk([2*j-1,2*j])) ;
%         c1(i+9)=-(lk([2*i-1,2*i])-lk([2*j-1,2*j]))'*(l([2*i-1,2*i])-l([2*j-1,2*j]))...
%         +5*D*norm(lk([2*i-1,2*i])-lk([2*j-1,2*j])) ;
nc1_sym{i,j}=gradient(c1(i,j),l');
   end
  end

nc_sym=nc_sym';
nc1_sym=nc1_sym';
nc = [nc_sym ;  nc1_sym];

end

