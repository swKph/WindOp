function c1=constraint1(l,lk)
D=100;
for i=1:9
   for j=1:9
c1(i)=-(lk([2*i-1,2*i])-lk([2*j-1,2*j]))'*(l([2*i-1,2*i])-l([2*j-1,2*j]))+5*D*norm(lk([2*i-1,2*i])-lk([2*j-1,2*j])) ;
c1(i+9)=-(lk([2*i-1,2*i])-lk([2*j-1,2*j]))'*(l([2*i-1,2*i])-l([2*j-1,2*j]))+5*D*norm(lk([2*i-1,2*i])-lk([2*j-1,2*j])) ;

   end
end
end