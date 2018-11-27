function [ncLk] = cal_nc(nc, l, lk)

lk1=lk(1); lk3=lk(3); lk5=lk(5); lk7=lk(7); lk9=lk(9); lk11=lk(11); lk13=lk(13); lk15=lk(15); lk17=lk(17);
lk2=lk(2); lk4=lk(4); lk6=lk(6); lk8=lk(8); lk10=lk(10); lk12=lk(12); lk14=lk(14); lk16=lk(16); lk18=lk(18);

l1=l(1); l3=l(3); l5=l(5); l7=l(7); l9=l(9); l11=l(11); l13=l(13); l15=l(15); l17=l(17);
l2=l(2); l4=l(4); l6=l(6); l8=l(8); l10=l(10); l12=l(12); l14=l(14); l16=l(16); l18=l(18);

c1 = nc(1,:);
c2 = nc(2:end,:)';

c1_val = zeros(9,18);
num_c = 8*9/2;
c2_val = zeros(num_c,18);


for i = 1:9
    c1_val(i,:) = eval(c1{i})';
end
% c1_val_column = reshape(c1_val',9*18,1);
% c1_val = c1_val';

count = 1;
for i=1:8
    for j=(i+1):9
        c2_val(count,:)=eval(c2{i,j})';
        count = count + 1;
    end
end
% c2_val_column = reshape(c2_val',num_c*18,1);
% c2_val = c2_val';

ncLk=[c1_val ; c2_val];


end

