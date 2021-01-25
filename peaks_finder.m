clc; clear;

t_final = 100;
res1 = zeros(t_final, 3);

for t = 1:1:t_final
    mat1 = csvread(strcat(num2str(t),'.csv'),1,0); % x versus h
    u_field = mat1(:,3);
    [umax,ind] = max(u_field);
    res1(t,3) = umax;
    res1(t,2) = mat1(ind,1);
    res1(t,1) = t;
end

csvwrite('max_u.csv',res1)