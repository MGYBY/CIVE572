clc;
clear;

mkdir contoursAndStreamlines

% Open "KEPE.txt" For Output As #1
c1 = 1E-27;
L = 4096;
Hmax = 8;
Hmin = 8;
Have = (Hmax + Hmin) / 2;

massflux = 500;

Qin = 12 * Have;
amplitude = 10;
% ps = L / 50;
g = 9.81;
f = 0 ; %0.01
cf = 0.007;
pi = 3.14159;
period = 200;
EF = 2;
Imax = 100 * EF;
Jmax = 100 * EF;
Ib = 2;
sim_time = 5000; 
selected_data = 20;

Hmax = 10;
dx = L / Imax;
dy = L / Jmax;
t = 0;
n = 0; % time steps
%%%%% Courant Number defined either by the flow velocity Uin
%%%%%%%%%%%%%%%%%%% or by the celerity of the wave c
c = sqrt(g * Hmax);
Co = 0.25;
%dt = Co / c / Sqr((1 / dx ^ 2) + (1 / dy ^ 2))
Uin = Qin/Have;
dt = Co / Uin / sqrt((1 / dx ^ 2) + (1 / dy ^ 2));

%%%%%% Initial Condtions %%%%%%
h = zeros(Imax+Ib+2, Jmax+Ib+2) + Have;
zo = zeros(Imax+Ib, Jmax+Ib);
hu = zeros(Imax+Ib, Jmax+Ib);
hv = zeros(Imax+Ib, Jmax+Ib);
Fuuh = zeros(Imax+Ib, Jmax+Ib);
Fuvh = zeros(Imax+Ib, Jmax+2+Ib);
Fvuh = zeros(Imax+2+Ib, Jmax+Ib);
Fvvh = zeros(Imax+Ib, Jmax+Ib);
qx = zeros(Imax + 2+ Ib, Jmax + Ib +1);
u = zeros(Imax + 2+ Ib, Jmax + Ib +1);
ua = zeros(Imax + 2+ Ib, Jmax + Ib);
qy = zeros(Imax + Ib+1, Jmax + 2 + Ib);
v = zeros(Imax + Ib+1, Jmax + 2 +Ib);
va = zeros(Imax + Ib, Jmax + 2 +Ib);

% %%%%% initial patch f pollutant of mass per unit area hc(i,j) = 1
hc = zeros(Imax+Ib+1, Jmax+Ib+1);
Fcuh = zeros(Imax+Ib+1, Jmax+Ib);
Fcvh = zeros(Imax+Ib, Jmax+Ib+1);

% hc(Imax / 10+Ib : Imax / 10 * 2+Ib, Jmax / 2 - Jmax / 5+Ib : Jmax+Ib) = 1.0;

%%%%% Inflow at the entrance %%%%%
% ind1 = 1 + Ib;
ind2 = 1 + Ib : Jmax / 2 + 2 + Ib;
qx(1:1+Ib, ind2) = Qin;
u = qx / Hmax;

pic_num=1;
tplot = 20.0; clf %, drawnow%, set(gcf,'renderer','zbuffer')
plotgap = round(tplot/dt); dt = tplot/plotgap;
nplots = round(sim_time/tplot);

xx = dx/2:dx:(dx/2+(Imax-1)*dx);
yy = dy/2:dy:(dy/2+(Jmax-1)*dy);
% yy = (dy/2+(Jmax-1)*dy):-dy: dy/2;
[xcoord,ycoord] = meshgrid(xx,yy);
[grad,im] = colorGradient([255,255,255],[255,0,0],20);

%time loop
for ii = 1:nplots
    %%%%%save data and GIF
    clf
    %hold on
    drawnow
    u_plot = u(1+Ib:Imax+Ib,1+Ib:Jmax+Ib);
    v_plot = v(1+Ib:Imax+Ib,1+Ib:Jmax+Ib);
    axis([0 1.05*max(xx) 0 1.05*max(yy)])
%     set(gcf,'Position',[100 100 1200 800])
    rectangle('Position',[0 0 L L],'LineWidth',2)
    hold on 
    p2 = contourf(xcoord, ycoord, hc(1:Imax, 1:Jmax)');
    colormap(grad)
    caxis([0, 1.0])
    aa =  get(gca,'XTickLabel');
   % set(gca,'fontsize',16)
%     axis([0 max(xx) 0 1.2*max(yy)])
    title(strcat('$t$=',num2str(t,3),'(s)'),'FontSize',14,'Interpreter','latex')
%     text(1.1*max(xx), 1.1*max(yy),strcat('$t$=',num2str(t,3),'(s)'),'FontSize',14,'Interpreter','latex')
    xlabel('\fontsize{12}{12}\selectfont {$x$(m)}','Interpreter','latex','FontWeight','bold')
    ylabel('\fontsize{12}{12}\selectfont {$y$(m)}','Interpreter','latex','FontWeight','bold')
    set(gca,'XTickLabel',aa,'fontsize',12)
    hold on
    p1 = streamslice(xcoord,ycoord, u_plot', v_plot', 8);
    hold off
    set(gca, 'Ydir', 'reverse')
    grid on
    saveas(gcf,strcat('./contoursAndStreamlines/',num2str(t),'.jpg'))
    clf;
   for nn = 1:plotgap
        t = t + dt;
        n = n+1;

    %         for i = 1 + Ib : Imax + Ib
    %         for j = 1 + Ib : Jmax + Ib
    %             if (h(i, j) > 0) 
    %                 vel = sqrt(u(i, j) ^ 2 + v(i, j) ^ 2);
    %                 Co = vel * dt * sqrt((1 / dx ^ 2) + (1 / dy ^ 2));
    % %                 if (Co > Comax)
    % %                     Comax = Co;
    % %                 end
    %             end
    %         end
    %         end

    %%%%%%Continuity equation %%%%%
        ind1 = 1 + Ib : Imax + Ib;
        ind2 = 1 + Ib : Jmax + Ib;
        h(ind1, ind2) = h(ind1, ind2) - dt * (qx(ind1+1, ind2) - qx(ind1, ind2)) / dx - dt * (qy(ind1, ind2 + 1) - qy(ind1, ind2)) / dy;

    %%%%% boundary conditions & suppress spurious oscillations
        %%% left wall
        qx(1:Ib+1, 1 + Ib : Jmax + Ib) = 0;
        %%% inflow B.C.
        qx(1:Ib+1, 1 + Ib : Jmax / 2 + 2 + Ib) = Qin;
        u(1:Ib+1, 1 + Ib : Jmax / 2 + 2 + Ib) = qx(1:Ib+1, 1 + Ib : Jmax / 2 + 2 + Ib) / Hmax;
        %%% inner velcity u
        ind2 = 1 + Ib : Jmax + Ib;
        ind1 = 2 + Ib : Imax + Ib;
        hu(ind1, ind2) = (h(ind1, ind2) + h(ind1-1, ind2)) / 2;
        u(ind1, ind2) = qx(ind1, ind2)./hu(ind1, ind2);
        %%% radiating B.C.
        u(Imax +Ib +1:Imax +Ib +2, 1 : Jmax) = sqrt(g / Hmin)*[(h(Imax+Ib, 1 : Jmax)-Hmin);(h(Imax+Ib, 1 : Jmax)-Hmin)]; 
        qx(Imax +Ib +1:Imax +Ib +2, 1 : Jmax) = [u(Imax+ Ib + 1, 1 : Jmax);u(Imax+ Ib + 1, 1 : Jmax)]*Hmin;
%         u(Imax +Ib +1:Imax +Ib +2, 1 : Jmax) = [u(Imax +Ib, 1 : Jmax);u(Imax +Ib, 1 : Jmax)];
%         qx(Imax +Ib +1:Imax +Ib +2, 1 : Jmax) = [qx(Imax +Ib, 1 : Jmax);qx(Imax +Ib, 1 : Jmax)];
        
        %%% top wall and bottom wall
        ind1 = 1 + Ib : Imax + Ib;
        qy(ind1, 1:Ib+1)=0;
        v(ind1, 1:Ib+1)=0;
        qy(ind1, Jmax+Ib+1:Jmax+Ib+2)=0;
        v(ind1, Jmax+Ib+1:Jmax+Ib+2)=0;
        
        qx(ind1, Jmax+Ib+1) = 2*qx(ind1, Jmax+Ib)-qx(ind1, Jmax+Ib-1);
        u(ind1, Jmax+Ib+1) = 2*u(ind1, Jmax+Ib)-u(ind1, Jmax+Ib-1);
%         %%% left and right 'wall'
%         ind1 = 1 + Ib : Imax + Ib;
%         qy(ind1, 1:Ib+1)=0;
%         v(ind1, 1:Ib+1)=0;
%         qy(ind1, Jmax+Ib+1: Jmax+Ib+2)=0;
%         v(ind1, Jmax+Ib+1: Jmax+Ib+2)=0;
        %%% inner velocity v
        ind1 = 1 + Ib : Imax + Ib;
        ind2 = 2 + Ib : Jmax + Ib;
        hv(ind1, ind2) = (h(ind1, ind2) + h(ind1, ind2 - 1)) / 2; 
        v(ind1, ind2) = qy(ind1, ind2) ./ hv(ind1, ind2);
        
        
    %%%%%%%%%%%%%%%% Evaluation of the advection terms for  momentum fluxes
    %%%%%%%%%%%%%%%% Fuuh, Fuvh, Fvuh, Fvvh %%%%%
        for i = 1 + Ib : Imax + Ib
        for j = 1 + Ib : Jmax + Ib
                %upwind scheme
                if (u(i, j) >= 0) 
                    phiu = u(i-1, j) * qx(i-1, j);
                    phic = u(i, j) * qx(i, j);
                    phid = u(i+1, j) * qx(i+1, j);
                else
                    phiu = u(i+2, j) * qx(i+2, j);
                    phic = u(i+1, j) * qx(i+1, j);
                    phid = u(i, j) * qx(i, j);
                end
                Fuuh(i, j) = FluxLim_MINMOD(phiu, phic, phid);
        end
        end

    %Fuvh(2 to Imax, 1 to Jmax+1) at the vorticity nodes-------
    %u(i,1) and u(i,Jmax+1) are specified at the boundary
    %Fuvh(i,1) is a special case because u(i,0) does not exist.
        Fuvh(2 + Ib : Imax + Ib, 1:Ib+1) = 0;         %exchange of momentum with the top wall
        Fuvh(2 + Ib : Imax + Ib, Jmax +Ib + 1) = 0;   %exchange of momentum with the bottom wall
        for i = 2 + Ib : Imax + Ib
        for j = 2 + Ib : Jmax + Ib
                qa = (qy(i, j) + qy(i - 1, j)) / 2;
                %upwind scheme
                if (qa >= 0) 
                    phiu = qa * u(i, j-2);
                    phic = qa * u(i, j-1);
                    phid = qa * u(i, j);
                else
                    phiu = qa * u(i, j+1);
                    phic = qa * u(i, j);
                    phid = qa * u(i, j-1);
                end
                Fuvh(i,j) = FluxLim_MINMOD(phiu, phic, phid);
        end
        end

    %Fvuh(1 to Imax+1,2 to Jmax) at the vorticity nodes-------
    %v(1,j) and v(Imax+1,j) are specified at the boundary
    %Fvuh(1,j) is a special case because v(0,j) does not exist.
        Fvuh(1 : Ib + 1, 2 + Ib : Jmax + Ib) = 0;   % exchange of momentum with the left wall
        Fvuh(Imax + Ib + 1, 2 + Ib : Jmax + Ib) = 0;   % exchange of momentum with the right wall
        for j = 2 + Ib : Jmax + Ib
        for i = 2 + Ib : Imax + Ib
                qa = (qx(i, j) + qx(i, j - 1)) / 2;
                %upwind scheme
                if (qa >= 0) 
                    phiu = qa * v(i-2, j);
                    phic = qa * v(i-1, j);
                    phid = qa * v(i, j);
                else
                    phiu = qa * v(i+1, j);
                    phic = qa * v(i, j);
                    phid = qa * v(i-1, j);
                end
                Fvuh(i,j) = FluxLim_MINMOD(phiu, phic, phid);
        end
        end

    %Fvvh(1 to Imax, 1 to Jmax)  at the h nodes--------------
    %qy(i,1), qy(i,Jmax+1), v(i,1) and v(i,Jmax+1) are specified at the boundary
%         for i = 1 + Ib : Imax + Ib
%                 qy(i, 1) = 0;                  %top
%                 qy(i, Jmax + 1) = 0;      %bottom
%         end
        ind1 = 1 + Ib : Imax + Ib;
        qy(ind1, 1:Ib+1) = 0;
        qy(ind1, Jmax+Ib+1) = 0;

        for j = 1 + Ib : Jmax + Ib
        for i = 1 + Ib : Imax + Ib
                %upwind scheme
                if (v(i, j) >= 0) 
                    phiu = qy(i, j-1) * v(i, j-1);
                    phic = qy(i, j) * v(i, j);
                    phid = qy(i, j+1) * v(i, j+1);
                else
                    phiu = qy(i, j+2) * v(i, j+2);
                    phic = qy(i, j+1) * v(i, j+1);
                    phid = qy(i, j) * v(i, j);
                end
                Fvvh(i,j) = FluxLim_MINMOD(phiu, phic, phid);
        end
        end

    %%%%%%% x-momentum equation %%%%%
        ind1 = 2 + Ib : Imax + Ib;
        ind2 = 1 + Ib : Jmax + Ib;
        adv = -dt * (Fuuh(ind1, ind2) - Fuuh(ind1 - 1, ind2)) / dx - dt * (Fuvh(ind1, ind2 + 1) - Fuvh(ind1, ind2)) / dy;
        Fx = -cf / 2 * u(ind1, ind2) .* sqrt(u(ind1, ind2) .^ 2 + v(ind1, ind2) .^ 2);
        qx(ind1, ind2) = qx(ind1, ind2)- g * dt * hu(ind1, ind2) .* (h(ind1, ind2) + zo(ind1, ind2) - h(ind1 - 1, ind2) - zo(ind1 - 1, ind2)) / dx + adv + Fx;
%         for i = 2 + Ib : Imax + Ib
%         for j = 1 + Ib : Jmax + Ib
%                 adv = -dt * (Fuuh(i, j) - Fuuh(i - 1, j)) / dx - dt * (Fuvh(i, j + 1) - Fuvh(i, j)) / dy;
%                 Fx = -cf / 2 * u(i, j) * sqrt(u(i, j) ^ 2 + v(i, j) ^ 2);
%                 qx(i, j) = qx(i, j) - g * dt * hu(i, j) * (h(i, j) + zo(i, j) - h(i - 1, j) - zo(i - 1, j)) / dx + adv + Fx;
%         end
%         end

    % qx(1,j) and qx(Imax+1,j) are boundary conditions.
    % Fuuh(1 to Imax,1 to Jmax) and Fuvh(2 to Imax, 1 to Jmax+1) are needed.

    %%%%% y-momentum equation
        ind1 = 1 + Ib : Imax + Ib;
        ind2 = 2 + Ib : Jmax + Ib;
        adv = -dt * (Fvuh(ind1 + 1, ind2) - Fvuh(ind1, ind2)) / dx - dt * (Fvvh(ind1, ind2) - Fvvh(ind1, ind2- 1)) / dy;
        Fy = -cf / 2 * v(ind1, ind2) .* sqrt(u(ind1, ind2) .^ 2 + v(ind1, ind2) .^ 2);
        qy(ind1, ind2) = qy(ind1, ind2)-g* dt * hv(ind1, ind2) .* (h(ind1, ind2) + zo(ind1, ind2) - h(ind1, ind2 - 1) - zo(ind1, ind2- 1)) / dy + adv + Fy;

%         for j = 2 + Ib : Jmax + Ib
%         for i = 1 + Ib : Imax + Ib
%                 adv = -dt * (Fvuh(i + 1, j) - Fvuh(i, j)) / dx - dt * (Fvvh(i, j) - Fvvh(i, j - 1)) / dy;
%                 Fy = -cf / 2 * v(i, j) * sqrt(u(i, j) ^ 2 + v(i, j) ^ 2);
%                 qy(i, j) = qy(i, j) - g * dt * hv(i, j) * (h(i, j) + zo(i, j) - h(i, j - 1) - zo(i, j - 1)) / dy + adv + Fy;
%         end
%         end
    % qy(i,1) and q(i,Jmax+1) are boundary conditions.
    % Fvuh(1 to Imax+1,2 to Jmax) and Fvvh(1 to Imax, 1 to Jmax) are needed.

    % %%%%% mass balance for c equation 
    %%%% Mass Flux
    for i = 1 + Ib : Imax + Ib
    for j = 1 + Ib : Jmax + Ib
            if (u(i, j) >= 0) 
                phiu = u(i, j) * hc(i - 2, j);
                phic = u(i, j) * hc(i - 1, j);
                phid = u(i, j) * hc(i, j);
            else
                phiu = u(i, j) * hc(i + 1, j);
                phic = u(i, j) * hc(i, j);
                phid = u(i, j) * hc(i - 1, j);
            end
            Fcuh(i, j) = FluxLim_MINMOD(phiu, phic, phid);

    end
    end
        %%%%% Outflow Boundary Condition 
    Fcuh(Imax + 1 + Ib:Imax + 2 + Ib, 1 + Ib : Jmax + Ib) = [Fcuh(Imax + Ib, 1 + Ib : Jmax + Ib); Fcuh(Imax + Ib, 1 + Ib : Jmax + Ib)];

    for i = 1 + Ib : Imax + Ib
    for j = 2 + Ib : Jmax + Ib
            if (v(i, j) >= 0) 
                phiu = v(i, j) * hc(i, j - 2);
                phic = v(i, j) * hc(i, j - 1);
                phid = v(i, j) * hc(i, j);
            else
                phiu = v(i, j) * hc(i, j + 1);
                phic = v(i, j) * hc(i, j);
                phid = v(i, j) * hc(i, j - 1);
            end
            Fcvh(i, j) = FluxLim_MINMOD(phiu, phic, phid);
    end
    end
    
    %%%%% Mass Balance
    ind1 = 1 + Ib : Imax + Ib;
    ind2 = 1 + Ib : Jmax + Ib;
    hc(ind1, ind2) = hc(ind1, ind2) - dt * (Fcuh(ind1+1, ind2) - Fcuh(ind1, ind2)) / dx - dt * (Fcvh(ind1, ind2+ 1) - Fcvh(ind1, ind2)) / dy;
    %%%%% Outflow Boundary Condition 
    hc(Imax + 1 + Ib:Imax + 2 + Ib, 1 + Ib : Jmax + Ib) = [hc(Imax + Ib, 1 + Ib : Jmax + Ib);hc(Imax + Ib, 1 + Ib : Jmax + Ib)];

    ind1 = 1 + Ib : 1 + EF + Ib;
    ind2 = 1 + Ib : Jmax / 2 + EF + Ib;
    hc(ind1, ind2) = hc(ind1, ind2) + massflux / 4 * dt / dx / EF / dy / EF;
    %         for i = 1 + Ib : Imax + Ib
    %         for j = 1 + Ib : Jmax + Ib
    %                 %upwind scheme
    %                 if (u(i, j) >= 0) 
    %                      Fcuh(i, j) = u(i, j) * hc(i - 1, j);
    %                 else
    %                      Fcuh(i, j) = u(i, j) * hc(i, j);
    %                 end
    %         end
    %         end

    % %%%%% Outflow Boundary Condition 
    %         i = Imax + 1 + Ib;
    %         for j = 1 + Ib : Jmax + Ib
    %                    Fcuh(i, j) = Fcuh(i - 1, j);
    %         end


    % % !!! periodic boundary conditon !!!
    % 'For j = 1 + Ib To Jmax + Ib
    %  '        Fcuh(Imax + 1 + Ib, j) = Fcuh(1 + Ib, j)
    %   '       Fcuh(Ib, j) = Fcuh(Ib + Imax, j)
    % 'Next

    %         for i = 1 + Ib : Imax + Ib
    %         for j = 2 + Ib : Jmax + Ib
    %                 % upwind scheme
    %                 if (qy(i, j) >= 0) 
    %                      Fcvh(i, j) = v(i, j) * hc(i, j - 1);
    %                 else
    %                      Fcvh(i, j) = v(i, j) * hc(i, j);
    %                 end
    %         end
    %         end

    % %%%%% Mass Balance 
    %         for i = 1 + Ib : Imax + Ib
    %             for j = 1 + Ib : Jmax + Ib
    %                 hc(i, j) = hc(i, j) - dt * (Fcuh(i + 1, j) - Fcuh(i, j)) / dx - dt * (Fcvh(i, j + 1) - Fcvh(i, j)) / dy;
    %             end
    %         end

    %%%%% Source Term - mass is added continuously at a rate
    %%%%% equal the mass flux ''''''''''''''''''''''''''''''
     % concentration c(i,j) mass per unit volume
     % hc(i,j) is mass per unit area
     % mass flux is mass per unit time
    %         for i = 1 + Ib : 1 + EF + Ib
    %             for j = Jmax / 2 - EF - 1 + Ib : Jmax / 2 + EF + Ib;
    %                         hc(i, j) = hc(i, j) + massflux / 4 * dt / dx / EF / dy / EF;
    %             end
    %         end




   end
end
    
