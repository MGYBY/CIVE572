clc;
clear;

%stable case Fr=1.15
hn = 0.0798;
un = 0.32176;
g = 9.810;
sine_theta = 0.05;
tan_theta = sine_theta/sqrt(1-sine_theta^2);
cf = g*tan_theta*2*hn/(un^2);
epsilon =  0.25; period = 0.933; %large amplitude but still damped
% period = [0.933, 1.25];
% period = 0.9:0.025:1.1;
% period = randn(1,5);

% %unstable case1 Fr=3.7
% hn = 0.0079756;
% un = 1.03774;
% g = 9.810;
% sine_theta = 0.05011;
% tan_theta = sine_theta/sqrt(1-sine_theta^2);
% cf = g*sine_theta*2*hn/(un^2);
% epsilon =  0.01; 
% period = 0.933;

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

Imax = 6000;
xx = 40; % Brock's setup
dx = xx/Imax;
q = zeros(1,Imax+1); q(1) = hn*un;
%q = q+hn*un;
h = zeros(1,Imax); h = h+hn;
hhat = zeros(1,Imax); uhat = zeros(1,Imax);
u = zeros(1,Imax+1);
%u = u+un;
sim_time = 200;
co = 0.10;
dt = co*dx/(1.0*un);
t = 0;

        
%OUTFLOW boundary conditions?????
q(Imax + 1) = sqrt(g*h(Imax))*(h(Imax)-hn);


x = 0.5*dx:dx:(0.5*dx+(Imax-1)*dx);
pic_num=1;
tplot = 0.2; clf %, drawnow%, set(gcf,'renderer','zbuffer')
plotgap = round(tplot/dt); dt = tplot/plotgap;
nplots = round(sim_time/tplot);

%time loop
for i = 1:nplots
   for n = 1:plotgap
        t = t+dt;
% %         Radiating B.C.
%         q(Imax + 1) = sqrt(h(Imax)*g)*(h(Imax)-hn);

%         Zero Gradient
        q(Imax + 1)=q(Imax);
        %continuity equation
         h(1) = hn*(1+epsilon*sum(sin(2*pi./period*t)));
%         rd = epsilon*(rand-0.5); % random dist.
%         h(1) = hn*(1+rd);
        h(2:end) = h(2:end) - dt/dx*(q(3:Imax+1)-q(2:Imax));
        

        hhat(2:Imax) = (h(2:Imax) + h(1:Imax-1)) / 2;
        u(2:Imax) = q(2:Imax) ./ hhat(2:Imax);
         u(1) = un*hn/h(1);
         for k = 1 : Imax
            if (u(k) > 0)
                uhat(k) = u(k);
            else
                uhat(k) = u(k + 1);
            end
         end
         
         qu(1:Imax) = uhat(1:Imax).^2.* h(1:Imax);
         frictionterm = cf * u .^ 2 ./ 2;
         q(2:Imax) = q(2:Imax) - g * dt / dx * hhat(2:Imax) .* (h(2:Imax) - h(1:Imax-1)) + (qu(1:Imax-1) - qu(2:Imax)) * dt / dx - frictionterm(2:Imax) * dt + g*tan_theta*hhat(2:Imax)*dt;
    end
    
    %%%%%save data and GIF
    %clf
    %hold on
    drawnow
    plot(x,h,'LineWidth', 1.75);
    aa =  get(gca,'XTickLabel');
   % set(gca,'fontsize',16)
    axis([0 max(x) 0 1.2*max(h(1:end-1))])
    text(0.5*max(x), 0.25*hn,strcat('$t$=',num2str(t,3),'(s)'),'FontSize',14,'Interpreter','latex')
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
    pic_num = pic_num + 1;
    
    if mod(i,10) == 0
        coord = [x',h'];
        csvwrite(strcat(num2str(t),'.csv'),coord)
    end
end