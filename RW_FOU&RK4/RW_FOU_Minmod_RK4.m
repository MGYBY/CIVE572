clc;
clear;

mkdir WSP

% %stable case Fr=1.15
% hn = 0.0798;
% un = 0.32176;
% g = 9.810;
% sine_theta = 0.05;
% tan_theta = sine_theta/sqrt(1-sine_theta^2);
% cf = g*tan_theta*2*hn/(un^2);
% epsilon =  0.1; 
% period = 0.933; %large amplitude but still damped
% % period = [0.933, 1.25];
% % period = 0.9:0.025:1.1;
% % period = randn(1,5);


%unstable case1 Fr=3.7
hn = 0.0079756;
un = 1.03774;
g = 9.810;
sine_theta = 0.05011;
tan_theta = sine_theta/sqrt(1-sine_theta^2);
cf = g*sine_theta*2*hn/(un^2);
epsilon =  0.01; 
period = 0.933*8;

% %unstable case2 Fr=5.6
% hn = 0.00533;
% un = 1.2805;
% cf = 9.85e-3;
% g = 9.810;
% sine_theta = 0.119528;
% tan_theta = sine_theta/sqrt(1-sine_theta^2);
% epsilon =  0.01; period = 1.015;

% %Coarsening Effect Imax=40000, xx=400
% hn = 0.0079756;
% un = 1.03774;
% cf = 7.2813e-3;
% g = 9.810;
% sine_theta = 0.05011;
% tan_theta = sine_theta/sqrt(1-sine_theta^2);
% epsilon =  0.01; %period = 0.933;
% period = 0.6:0.3:1.8;

Imax = 4000;
xx = 200; % long channel
dx = xx/Imax;
Ib = 2;
q = zeros(1,Imax+Ib+2); q(1) = hn*un;
h = zeros(1,Imax+Ib)+hn; % h = h+hn;
hhat = zeros(1,Imax+Ib); 
uhat = zeros(1,Imax+Ib);
u = zeros(1,Imax+2+Ib);
qu = zeros(1,Imax+Ib+1);
qmed = q;
hmed = h;
sim_time = 200;
co = 0.18;
dt = co*dx/(1.0*un);
t = 0;
        
%OUTFLOW boundary conditions?????
q(Imax + 1) = q(Imax);
q(Imax + 2) = q(Imax + 1);

x = 0.5*dx:dx:(0.5*dx+(Imax-1)*dx); % computational domain: 1+Ib~Imax+Ib
pic_num=1;
tplot = 1; clf %, drawnow%, set(gcf,'renderer','zbuffer')
plotgap = round(tplot/dt); dt = tplot/plotgap;
nplots = round(sim_time/tplot);

%time loop
for i = 1:nplots
   for n = 1:plotgap
        t = t+dt;
% %         Radiating B.C.
%         q(Imax + 1) = sqrt(h(Imax)*g)*(h(Imax)-hn);

% %         Zero Gradient
%         q(Imax + 1:Imax+1+Ib)=sqrt(g*h(Imax))*(h(Imax)-hn);
%         u(Imax + 1:Imax+1+Ib)=sqrt(g/h(Imax))*(h(Imax)-hn);
%         h(Imax+1) = hn;
        q(Imax+1) = q(Imax);
        q(Imax + 2)= q(Imax+1);
%         h(Imax+1) = h(Imax);
        u(Imax+1) = u(Imax);
        u(Imax+2) = u(Imax);

        %continuity equation
        h(1) = hn*(1+epsilon*sum(sin(2*pi./period*t)));
%         h(Ib) = h(1);
        u(1) = un*hn/h(1);
        q(1) = un*hn;
%         qu(1) = un*hn*u(1);
        qu(1) = q(1)*u(1);
%         q(Ib) = q(1);
        ind1 = Ib:Imax;
        ind2 = Ib:Imax;
        k1A = fscalc_continuity(ind1, q, dx);
        k1B = fscalc_momentum(ind2, u, h, q, dx, tan_theta, cf, un, hn, uhat, hhat);
        qmed(ind2) = q(ind2)+k1B*dt/2;
        hmed(ind1) = h(ind1)+k1A*dt/2;
        k2A = fscalc_continuity(ind1, qmed, dx);
        k2B = fscalc_momentum(ind2, u, hmed, qmed, dx, tan_theta, cf, un, hn, uhat, hhat);
        qmed(ind2) = q(ind2)+k2B*dt/2;
        hmed(ind1) = h(ind1)+k2A*dt/2;
        k3A = fscalc_continuity(ind1, qmed, dx);
        k3B = fscalc_momentum(ind2, u, hmed, qmed, dx, tan_theta, cf, un, hn, uhat, hhat);
        qmed(ind2) = q(ind2)+k3B*dt/2;
        hmed(ind1) = h(ind1)+k3A*dt/2;
        k4A = fscalc_continuity(ind1, qmed, dx);
        k4B = fscalc_momentum(ind2, u, hmed, qmed, dx, tan_theta, cf, un, hn, uhat, hhat);
        h(ind1) = h(ind1) + dt/6*(k1A+2.*k2A+2.*k3A+k4A);
        q(ind2) = q(ind2) + dt/6*(k1B+2.*k2B+2.*k3B+k4B);
       
   end
    
    %%%%%save data and GIF
    %clf
    %hold on
    drawnow
    plot(x(Ib:Imax),h(Ib:Imax),'LineWidth', 1.75);
    aa =  get(gca,'XTickLabel');
   % set(gca,'fontsize',16)
    axis([0 max(x) 0 1.1*max(h(1:end-1))])
%     text(0.5*max(x), 0.25*hn,strcat('$t$=',num2str(t,3),'(s)'),'FontSize',14,'Interpreter','latex')
    title(strcat('$t$=',num2str(t,3),'(s)'),'FontSize',14,'Interpreter','latex');
    xlabel('\fontsize{12}{12}\selectfont {$x$ (m)}','Interpreter','latex','FontWeight','bold')
    ylabel('\fontsize{12}{12}\selectfont {$h$(m)}','Interpreter','latex','FontWeight','bold')
    set(gca,'XTickLabel',aa,'fontsize',12)
    
    grid on
    
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1 %save the initial GIF 
        imwrite(I,map,strcat('Case2','.gif'),'gif', 'Loopcount',inf,'DelayTime',0.2);
    else %append additional frames
        imwrite(I,map,strcat('Case2','.gif'),'gif','WriteMode','append','DelayTime',0.2);
    end
    saveas(gcf,strcat('./WSP/',num2str(t),'.jpg'))
    pic_num = pic_num + 1;
    
    if mod(i,10) == 0
        coord = [x(Ib:Imax)',h(Ib:Imax)'];
        csvwrite(strcat(num2str(t),'.csv'),coord)
    end
end