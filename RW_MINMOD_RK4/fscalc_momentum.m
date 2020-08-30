function [fs2] = fscalc_momentum(ind2, u, h, q, dx, tan_theta, cf, un, hn, uhat, hhat)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
g = 9.81;
hhat(ind2) = (h(ind2) + h(ind2-1)) / 2;
u(ind2) = q(ind2) ./ hhat(ind2);
u(1) = un*hn/h(1);
uhat(1) = u(1);
% radiating B.C.
q(ind2(end)+1) = q(ind2(end));
q(ind2(end)+2) = q(ind2(end));
u(ind2(end)+1) = u(ind2(end));
u(ind2(end)+2) = u(ind2(end));

for k = ind2(2):ind2(end)
 ua = (u(k)+u(k+1))/2;
    if (ua >= 0)
        phiu = u(k-1);
        phic = u(k);
        phid = u(k+1);
    else
        phiu = u(k+2);
        phic = u(k+1);
        phid = u(k);
    end
    uhat(k) = FluxLim_MINMOD(phiu, phic, phid);
end

qu = uhat.^2.* h;
qu(ind2(1)) = u(ind2(1))^2*h(ind2(1)); % second cell upwind scheme
frictionterm = cf * u .^ 2 ./ 2;
fs2 = -g*(1/dx)*hhat(ind2).*(h(ind2)-h(ind2-1))-(qu(ind2)-qu(ind2-1))/ dx...
    -frictionterm(ind2)+g*tan_theta*hhat(ind2);
end

