clear;
clc;

xmax = 1000;                                         
Ho = 10;                                                  
g = 9.81;
n = 0.014;
t = 0;
dx = 10.0;
c = sqrt(g * Ho);                                     
Courant = 0.40;                                      
dt = Courant * dx / c;
Imax = round(xmax / dx);
sim_time = 1500;

Fr = zeros(1,Imax);
h = zeros(1, Imax);
u = zeros(1,Imax+1);
q = zeros(1,Imax+1);
uhat = zeros(1, Imax);
Fuuh = zeros(1, Imax);

% bottom elevation initial zero
zo(1:Imax) = 0;

% bottom elevation of the dam - sinusoidal from istart to (istart+ilength)
istart = Imax / 5;
ilength = Imax / 5;
zo(istart : istart + ilength) = 8 * sin((0 : (ilength)) / ilength * pi);

% Initial discharge is zero. Initial depth plus bottom elevation is constant.
h(1 : Imax) = Ho - zo(1 : Imax);
q(1 : Imax) = 0;

x = 0.5*dx:dx:(0.5*dx+(Imax-1)*dx);
pic_num=1;
tplot = 1; clf %, drawnow%, set(gcf,'renderer','zbuffer')
plotgap = round(tplot/dt); dt = tplot/plotgap;
nplots = round(sim_time/tplot);

% time loop
for ii = 1:nplots
   for nn = 1:plotgap
        t = t + dt;
        % INFLOW boundary conditions
        q(1) = 10.0;
        h(1) = Ho;
        u(1) = q(1)/Ho;

        % OUTFLOW boundary conditions
        % change the weir height from 0.10*Ho to 0.70*Ho to see
        % how the downstrem depth may alter the flow profile
        weirheight = 0.05 * Ho;
        q(Imax + 1) = sqrt(g * h(Imax)) * (h(Imax) - weirheight);
        u(Imax + 1) = sqrt(g / h(Imax)) * (h(Imax) - weirheight);

        % continuity equation
        h(2:Imax) = h(2:Imax) - dt / dx * (q(3:Imax + 1) - q(2:Imax));
        
        % momentum equation
        hhat(2:Imax) = (h(2:Imax) + h(1:Imax - 1)) / 2;
        % gain computational stability by division of an averaged depth
        u(2:Imax) = q(2:Imax) ./ hhat(2:Imax);
        Fr(2:Imax) = u(2:Imax) ./ sqrt(g * hhat(2:Imax));
        for k = 2 : Imax
        % upwind scheme
            if ((u(k) + u(k + 1)) > 0) 
                uhat(k) = u(k);
            else
                uhat(k) = u(k + 1);
            end
        % central difference
        % uhat(i) = (u(i) + u(i + 1)) / 2
        % momentum flux on the plus side of the finite volume
        Fuuh(k) = uhat(k)^2 * h(k);
        end
        Fuuh(1) = h(1)*u(1)^2;

        frictionterm = g * n ^ 2 * q(2:Imax) .*(q(2:Imax)) ./ (hhat(2:Imax) .^ (7 / 3));
        q(2:Imax) = q(2:Imax) - g * dt / dx * hhat(2:Imax) .* ((h(2:Imax) + zo(2:Imax)) ...
            - (h(1:Imax-1) + zo(1:Imax-1))) + (Fuuh(1:Imax-1) - Fuuh(2:Imax)) * dt / dx - frictionterm * dt;

   end
   
   %%%%%save data and GIF
    %clf
    %hold on
    drawnow
    plot(x,h+zo,'LineWidth', 1.75)
    hold on 
    plot(x,zo,'LineWidth', 2.5)
    hold on
    plot(x,Fr,'LineWidth', 1.75,'Color', 'b')
    hold off
    aa =  get(gca,'XTickLabel');
   % set(gca,'fontsize',16)
    axis([0 max(x) 0 1.2*max(h(1:end-1))])
    title(strcat('$t$=',num2str(t,3),'(s)'),'FontSize',14,'Interpreter','latex')
    xlabel('\fontsize{12}{12}\selectfont {$x$ (m)}','Interpreter','latex','FontWeight','bold')
    ylabel('\fontsize{12}{12}\selectfont {$h$(m)}','Interpreter','latex','FontWeight','bold')
    set(gca,'XTickLabel',aa,'fontsize',12)
    
    grid on
end

        