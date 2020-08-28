function [fs2] = fscalc_momentum(ind2, u, h, q, dx, tan_theta, cf, un, hn, uhat, hhat)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
g = 9.81;
hhat(ind2) = (h(ind2) + h(ind2-1)) / 2;
u(ind2) = q(ind2) ./ hhat(ind2);
u(1) = un*hn/h(1);
uhat(1) = u(1);
for k = ind2
    if (u(k) > 0)
        uhat(k) = u(k);
    else
        uhat(k) = u(k + 1);
    end
end

qu = uhat.^2.* h;
frictionterm = cf * u .^ 2 ./ 2;
fs2 = -g*(1/dx)*hhat(ind2).*(h(ind2)-h(ind2-1))-(qu(ind2)-qu(ind2-1))/ dx...
    -frictionterm(ind2)+g*tan_theta*hhat(ind2);
end

