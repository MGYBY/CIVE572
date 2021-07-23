h_series = [0.0282034071 0.0299590483 0.0312259761 0.0319488122 0.0324026831 0.0326779068 0.0328476701 0.0329167314];
n_series = [1000 2000 4000 8000 16000 32000 64000 128000];
delta_x_series = 160.0./n_series;
num_n = length(n_series);

p_order = zeros(1, (num_n-2));
h_extrap = zeros(1, (num_n-2));
FE_orig = zeros(1, (num_n-2));
FE_mod = zeros(1, num_n);

for l=2:(num_n-1)
    % for each tri-value group
    group_num = l-1;
    r = delta_x_series(l)/delta_x_series(l+1);
    p_order(group_num) = 1.0/log(r)*log((h_series(l)-h_series(l-1))/(h_series(l+1)-h_series(l)));
    h_extrap(group_num) = (r^(p_order(group_num))*h_series(l+1)-h_series(l))/(r^(p_order(group_num))-1);
    FE_orig(group_num) = abs(h_series(l)-h_extrap(group_num))/h_extrap(group_num);
end

h_final_extrap = h_extrap(end);
for l=1:num_n
    FE_mod(l) = abs(h_series(l)-h_final_extrap)/h_final_extrap;
end

dlmwrite('result_orig.txt', [p_order', h_extrap', FE_orig'],'delimiter', ' ');
dlmwrite('result_mod.txt', [h_series', n_series', FE_mod'],'delimiter', ' ');
