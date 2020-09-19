
%%%%% in time loop:


clc;
clear;
mkdir WSP2

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

%%%%% non-staggered grid
%unstable case1 Fr=3.7
hn = 0.0079756;
un = 1.03774;
qn = hn*un;
g = 9.810;
sine_theta = 0.05011;
tan_theta = sine_theta/sqrt(1-sine_theta^2);
cf = g*sine_theta*2*hn/(un^2);
epsilon =  0.01; 
period = 0.933;

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

Imax = 8000;
xx = 40; % Brock's setup
dx = xx/Imax;
q = zeros(1,Imax+1)+qn; % q(1) = hn*un;
h = zeros(1,Imax+1)+hn; % h = h+hn;
u = zeros(1,Imax+1)+un;
hstar = zeros(1,Imax);
% left and right flux values
aleft = zeros(1,Imax); % geopotential following flux index
aright = zeros(1,Imax); % geopotential following flux index
UL = zeros(2,Imax);
UR = zeros(2,Imax);
ULstar = zeros(2,Imax); % approximate states for 2 eqns
URstar = zeros(2,Imax);
FL = zeros(2,Imax); %FL and FR follow flux index
FR = zeros(2,Imax);
% qq follow flux index
qqleft = zeros(1,Imax);
qqright = zeros(1,Imax);
% wave speed SL, SR and S* following flux index
SL = zeros(1,Imax);
SR = zeros(1,Imax);
Sstar = zeros(1,Imax);
%HLLC flux (for 2 eqns)
Fhllc = zeros(2,Imax);
sim_time = 200;
co = 0.1;
dt = co*dx/(1.0*un);
t = 0;
        
%OUTFLOW boundary conditions(supersonic)
q(Imax + 1) = q(Imax);
h(Imax+1) = h(Imax);
u(Imax+1) = u(Imax);

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

        %Zero Gradient (supersonic outlet)
        q(Imax+1) = q(Imax);
        h(Imax+1) = h(Imax);
        %inlet with dist.
        h(1) = hn*(1+epsilon*sum(sin(2*pi./period*t)));
        q(1) = qn;
        u = q./h;
        
        % calc. left and right flux val.
        UL = [h(1:Imax);q(1:Imax)];
        UR = [h(2:Imax+1);q(2:Imax+1)];
        FL = [q(1:Imax);(q(1:Imax).^2./h(1:Imax)+1/2*g*h(1:Imax).^2)];
        FR = [q(2:Imax+1);(q(2:Imax+1).^2./h(2:Imax+1)+1/2*g*h(2:Imax+1).^2)];
        aleft = sqrt(g*h(1:Imax));
        aright = sqrt(g*h(2:Imax+1));
        % solve the hyperbolic eqns (without source terms)
        for k = 1:Imax % at flux k
            % calc. estimate for h*
            hstar(k) = 1/2*(h(k+1)+h(k))-1/4*(u(k+1)-u(k))*(h(k+1)+h(k))/(aleft(k)+aright(k));
            %%% calc. the estimated wave speed SL, SR, S*
            % left value
            if hstar(k)>h(k)
                qqleft(k) = sqrt(1/2*((hstar(k)+h(k))*hstar(k))/h(k)^2);
            else
                qqleft(k) = 1.0;
            end
            SL(k) = u(k)-aleft(k)*qqleft(k);
            % right value
            if hstar(k)>h(k+1)
                qqright(k) = sqrt(1/2*((hstar(k)+h(k+1))*hstar(k))/h(k+1)^2);
            else
                qqright(k) = 1.0;
            end
            SR(k) = u(k+1)-aright(k)*qqright(k);
            Sstar(k) = (SL(k)*h(k+1)*(u(k+1)-SR(k))-SR(k)*h(k)*(u(k)-SL(k)))...
                /(h(k+1)*(u(k+1)-SR(k))-h(k)*(u(k)-SL(k)));
            %%% calc. approx states for U*L, U*R
            ULstar(:,k) = h(k)*((SL(k)-u(k))/(SL(k)-Sstar(k)))*[1.0;Sstar(k)];
            URstar(:,k) = h(k+1)*((SR(k)-u(k+1))/(SR(k)-Sstar(k)))*[1.0;Sstar(k)];
            %%% evaluate fluxes Fhllc for both cont. and momentum eqns
            if SL(k)>=0.0
                Fhllc(:,k) = FL(:,k);
            elseif (SL(k)<0.0 && Sstar(k)>0.0)
                Fhllc(:,k) = FL(:,k)+SL(k).*(ULstar(:,k)-UL(:,k));
            elseif (Sstar(k)<0.0 && SR(k)>=0.0)
                Fhllc(:,k) = FR(:,k)+SR(k).*(URstar(:,k)-UR(:,k));
            else
                Fhllc(:,k) = FR(:,k);
            end
        end
        %%% solve the conserved form and obtain h and q (namely U)
        % hyperbolic problem
        h(2:Imax) = h(2:Imax)-dt/dx*(Fhllc(1, 2:Imax)-Fhllc(1, 1:Imax-1));
        q(2:Imax) = q(2:Imax)-dt/dx*(Fhllc(2, 2:Imax)-Fhllc(2, 1:Imax-1));
        %%% solve the ODE for time integration (RK3 or RK4)
        qmed = q(2:Imax)+dt*(source_term(h(2:Imax), q(2:Imax), tan_theta, cf));
        qmed = 3.0/4.0*q(2:Imax)+1.0/4.0*(qmed+dt*source_term(h(2:Imax), qmed, tan_theta, cf));
        q(2:Imax) = 1.0/3.0*q(2:Imax)+2.0/3.0*(qmed+dt*source_term(h(2:Imax), qmed, tan_theta, cf));
   end
    %%%%% save data and GIF
    %clf
    %hold on
    drawnow
    plot(x(1:Imax),h(1:Imax),'LineWidth', 1.75);
    aa =  get(gca,'XTickLabel');
   % set(gca,'fontsize',16)
    axis([0 max(x) 0 1.1*max(h(1:end-1))])
%     text(0.5*max(x), 0.25*hn,strcat('$t$=',num2str(t,3),'(s)'),'FontSize',14,'Interpreter','latex')
    title(strcat('$t$=',num2str(t,3),'(s)'),'FontSize',14,'Interpreter','latex');
    xlabel('\fontsize{12}{12}\selectfont {$x$ (m)}','Interpreter','latex','FontWeight','bold')
    ylabel('\fontsize{12}{12}\selectfont {$h$(m)}','Interpreter','latex','FontWeight','bold')
    set(gca,'XTickLabel',aa,'fontsize',12)
    
    grid on
    saveas(gcf,strcat('./WSP2/',num2str(t),'.jpg'))
    
    
    if mod(i,10) == 0
        coord = [x(1:Imax)',h(1:Imax)'];
        csvwrite(strcat('./WSP2/',num2str(t),'.csv'),coord)
    end
end