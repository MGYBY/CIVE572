clc;
clear;

mkdir U_data % data for xcoord_U_V_H
mkdir contoursAndStreamlines % streamlines graphics

L = 4096;
Hmax = 8;
Hmin = 8;
Have = (Hmax + Hmin) / 2;

massflux = 500;

Qin = 12 * Have;
amplitude = 10;
g = 9.81;
f = 0 ; %0.01
cf = 0.007;
pi = 3.14159;
period = 200;
EF = 2;
Imax = 100 * EF;
Jmax = 100 * EF;
Ib = 0;
sim_time = 5000; 
tsteps = 0;
selected_data = 20;

Hmax = 10;
dx = L / Imax;
dy = L / Jmax;
t = 0;

%%%%% Courant Number defined either by the flow velocity Uin
%%%%%%%%%%%%%%%%%%% or by the celerity of the wave c
% c = sqrt(g * Hmax);
Co = 0.20;
%dt = Co / c / Sqr((1 / dx ^ 2) + (1 / dy ^ 2))
Uin = Qin/Have;
dt = Co / Uin / sqrt((1 / dx ^ 2) + (1 / dy ^ 2))/1.1;
    
%%%%%% Initial Condtions %%%%%%
h = zeros(Imax+Ib, Jmax+Ib) + Have;
zo = zeros(Imax+Ib, Jmax+Ib);
hu = zeros(Imax+Ib, Jmax+Ib);
hv = zeros(Imax+Ib, Jmax+Ib);

Fuuh = zeros(Imax+Ib, Jmax+Ib);
Fuvh = zeros(Imax+Ib, Jmax+1+Ib);
Fvuh = zeros(Imax+1+Ib, Jmax+Ib);
Fvvh = zeros(Imax+Ib, Jmax+Ib);

qx = zeros(Imax + 1+ Ib, Jmax + Ib);
u = zeros(Imax + 1+ Ib, Jmax + Ib);
ua = zeros(Imax + 1+ Ib, Jmax + Ib);


% %%%%% initial patch f pollutant of mass per unit area hc(i,j) = 1
hc = zeros(Imax, Jmax);
Fcuh = zeros(Imax+1, Jmax);
Fcvh = zeros(Imax, Jmax+1);
% ind1 = Imax / 10 : Imax / 10 * 2;
% ind2 = Jmax / 2 - Jmax / 5 : Jmax;
% hc(ind1, ind2) = 1.0;

%%%%% Inflow at the entrance %%%%%
qx(1 + Ib, 1 + Ib : Jmax / 2 + 2 + Ib) = Qin;
u(1 + Ib, 1 + Ib : Jmax / 2 + 2 + Ib) = qx(1 + Ib, 1 + Ib : Jmax / 2 + 2 + Ib)/Hmax;

qy = zeros(Imax + Ib, Jmax +1 + Ib);
v = zeros(Imax + Ib, Jmax + 1 +Ib);
va = zeros(Imax + Ib, Jmax + 1 +Ib);

pic_num=1;
tplot = 20.0; clf %, drawnow%, set(gcf,'renderer','zbuffer')
plotgap = round(tplot/dt); dt = tplot/plotgap;
nplots = round(sim_time/tplot);

xx = dx/2:dx:(dx/2+(Imax-1)*dx);
yy = dy/2:dy:(dy/2+(Jmax-1)*dy);
[xcoord,ycoord] = meshgrid(xx,yy);
[grad,im] = colorGradient([255,255,255],[255,0,0],20);

%time loop
for ii = 1:nplots
   for nn = 1:plotgap
        t = t + dt;
        tsteps = tsteps + 1;

%%%%%%Continuity equation %%%%%
        ind1 = 1 + Ib : Imax + Ib; 
        ind2 = 1 + Ib : Jmax + Ib;
        h(ind1, ind2) = h(ind1, ind2)-dt/dx*(qx(ind1 + 1, ind2)-qx(ind1, ind2))-dt/dy*(qy(ind1, ind2+1)-qy(ind1, ind2));

%%%%% boundary conditions & suppress spurious oscillations
        ind2 = 1 + Ib : Jmax + Ib;
        ind1 = 2 + Ib : Imax + Ib;
        hu(ind1, ind2) = (h(ind1, ind2)+h(ind1-1 ,ind2))/2;
        u(ind1, ind2) = qx(ind1, ind2)./hu(ind1, ind2);

        ind1 = 1 + Ib : Imax + Ib;
        v(ind1, 1)=0; 
        v(ind1, Jmax+1) = 0;
        ind2 = 2 + Ib : Jmax + Ib;
        hv(ind1, ind2) = (h(ind1, ind2)+h(ind1, ind2-1))/2;
        v(ind1, ind2) = qy(ind1, ind2)./hv(ind1, ind2);
        
        qx(1,1 + Ib : Jmax + Ib) = 0; % left wall
        
        qx(1 + Ib, 1 + Ib : Jmax / 2 + 2 + Ib) = Qin; % inflow B.C.
        u(1 + Ib, 1 + Ib : Jmax / 2 + 2 + Ib) = qx(1 + Ib, 1 + Ib : Jmax / 2 + 2 + Ib)/Hmax;
        
        u(Imax + 1, 1 : Jmax) = sqrt(g / Hmin)*(h(Imax, 1 : Jmax)-Hmin); % exit and rad. B.C.
        qx(Imax + 1, 1 : Jmax) = u(Imax + 1, 1 : Jmax)*Hmin;

%%%%%%%%%%%%%%%% Evaluation of the advection terms for  momentum fluxes
%%%%%%%%%%%%%%%% Fuuh, Fuvh, Fvuh, Fvvh %%%%%
        ind1 = 1 + Ib : Imax + Ib;
        ind2 = 1 + Ib : Jmax + Ib;
        ua = (u(ind1, ind2)+u(ind1+1, ind2))/2;
        u_hat = u(ind1, ind2);
        u_ref = ua<0;
        u_move = u(ind1+1, ind2);
        u_hat(u_ref) = u_move(u_ref); % Fuuh(i, j) = qx(i + 1, j) * u(i + 1, j);
        Fuuh(ind1, ind2) = u_hat.^2.*h(ind1, ind2);
        
%Fuvh(2 to Imax, 1 to Jmax+1) at the vorticity nodes-------
%u(i,1) and u(i,Jmax+1) are specified at the boundary
%Fuvh(i,1) is a special case because u(i,0) does not exist.
        ind1 = 2 + Ib : Imax + Ib; 
        Fuvh(ind1, 1) = 0;
        Fuvh(ind1, Jmax+1) = 0;
        for i = 2 + Ib : Imax + Ib
%                 Fuvh(i, 1) = 0;               %exchange of momentum with the top wall
%                 Fuvh(i, Jmax + 1) = 0;   %exchange of momentum with the bottom wall
        for j = 2 + Ib : Jmax + Ib
                qa = (qy(i, j) + qy(i - 1, j)) / 2;
                %upwind scheme
                Fuvh(i,j) = qa*u(i, j+(-1)*(qa>=0));
%                 if (qa >= 0) 
%                        Fuvh(i, j) = qa * u(i, j - 1);
%                 else
%                        Fuvh(i, j) = qa * u(i, j);
%                 end
        end
        end
        
%Fvuh(1 to Imax+1,2 to Jmax) at the vorticity nodes-------
%v(1,j) and v(Imax+1,j) are specified at the boundary
%Fvuh(1,j) is a special case because v(0,j) does not exist.
        ind2 = 2 + Ib : Jmax + Ib;
        Fvuh(1, ind2) = 0;
        Fvuh(Imax+1, ind2) = 0;
        for j = 2 + Ib : Jmax + Ib
%                 Fvuh(1, j) = 0;              % exchange of momentum with the left wall
%                 Fvuh(Imax + 1, j) = 0;   % exchange of momentum with the right wall
        for i = 2 + Ib : Imax + Ib
                qa = (qx(i, j) + qx(i, j - 1)) / 2;
                %upwind scheme
                Fvuh(i,j) = qa*v(i+(-1)*(qa>=0),j);
%                 if (qa >= 0) 
%                        Fvuh(i, j) = qa * v(i - 1, j);
%                 else
%                        Fvuh(i, j) = qa * v(i, j);
%                 end
        end
        end
        
%Fvvh(1 to Imax, 1 to Jmax)  at the h nodes--------------
%qy(i,1), qy(i,Jmax+1), v(i,1) and v(i,Jmax+1) are specified at the boundary
        ind1 = 1 + Ib : Imax + Ib;
        qy(ind1, 1)=0;
        v(ind1, 1)=0;
        qy(ind1, Jmax+1)=0;
        v(ind1, Jmax+1)=0;
%         for i = 1 + Ib : Imax + Ib
%                 qy(i, 1) = 0;                  %top
%                 qy(i, Jmax + 1) = 0;      %bottom
%         end

        ind1 = 1 + Ib : Imax + Ib; 
        ind2 = 1 + Ib : Jmax + Ib;
        va = (v(ind1, ind2)+v(ind1, ind2+1))/2;
        v_hat = v(ind1, ind2);
        v_ref = va<0;
        v_move = v(ind1, ind2+1);
        v_hat(v_ref) = v_move(v_ref); % Fvvh(i, j) = qx(i, j+1) * u(i, j+1);
        Fvvh(ind1, ind2) = v_hat.^2.*h(ind1, ind2);
        
%         for j = 1 + Ib : Jmax + Ib
%         for i = 1 + Ib : Imax + Ib
%                 %upwind scheme
%                 if (v(i, j) >= 0) 
%                        Fvvh(i, j) = qy(i, j) * v(i, j);
%                 else
%                        Fvvh(i, j) = qy(i, j + 1) * v(i, j + 1);
%                 end
%         end
%         end

%%%%%%% x-momentum equation %%%%%

        ind1 = 2 + Ib : Imax + Ib;
        ind2 = 1 + Ib : Jmax + Ib;
        adv = -dt * (Fuuh(ind1, ind2) - Fuuh(ind1 - 1, ind2)) / dx - dt * (Fuvh(ind1, ind2 + 1) - Fuvh(ind1, ind2)) / dy;
        Fx = -cf / 2 * u(ind1, ind2) .* sqrt(u(ind1, ind2) .^ 2 + v(ind1, ind2) .^ 2);
        qx(ind1, ind2) = qx(ind1, ind2)- g * dt * hu(ind1, ind2) .* (h(ind1, ind2) + zo(ind1, ind2) - h(ind1 - 1, ind2) - zo(ind1 - 1, ind2)) / dx + adv + Fx;

    % qx(1,j) and qx(Imax+1,j) are boundary conditions.
    % Fuuh(1 to Imax,1 to Jmax) and Fuvh(2 to Imax, 1 to Jmax+1) are needed.

%%%%% y-momentum equation
        ind1 = 1 + Ib : Imax + Ib;
        ind2 = 2 + Ib : Jmax + Ib;
        adv = -dt * (Fvuh(ind1 + 1, ind2) - Fvuh(ind1, ind2)) / dx - dt * (Fvvh(ind1, ind2) - Fvvh(ind1, ind2- 1)) / dy;
        Fy = -cf / 2 * v(ind1, ind2) .* sqrt(u(ind1, ind2) .^ 2 + v(ind1, ind2) .^ 2);
        qy(ind1, ind2) = qy(ind1, ind2)-g* dt * hv(ind1, ind2) .* (h(ind1, ind2) + zo(ind1, ind2) - h(ind1, ind2 - 1) - zo(ind1, ind2- 1)) / dy + adv + Fy;

    % qy(i,1) and q(i,Jmax+1) are boundary conditions.
    % Fvuh(1 to Imax+1,2 to Jmax) and Fvvh(1 to Imax, 1 to Jmax) are needed.
    
% %%%%% mass balance for c equation 
        ind1 = 2 + Ib : Imax + Ib;
        ind2 = 1 + Ib : Jmax + Ib;
        ua = (u(ind1, ind2)+u(ind1+1, ind2))/2;
        hc_hat = hc(ind1, ind2);
        hc_move = hc(ind1-1, ind2);
        hc_hat(ua>=0) =  hc_move(ua>=0);
        Fcuh(ind1, ind2) = ua.*hc_hat;

%%%%% Outflow Boundary Condition 
        ind2 = 1 + Ib : Jmax + Ib;
        Fcuh(Imax+1+Ib, ind2) = Fcuh(Imax+Ib, ind2);


% % !!! periodic boundary conditon !!!
% 'For j = 1 + Ib To Jmax + Ib
%  '        Fcuh(Imax + 1 + Ib, j) = Fcuh(1 + Ib, j)
%   '       Fcuh(Ib, j) = Fcuh(Ib + Imax, j)
% 'Next
        ind1 = 1 + Ib : Imax + Ib;
        ind2 = 2 + Ib : Jmax + Ib;
        va = (v(ind1, ind2)+v(ind1, ind2+1))/2;
        hc_hat = hc(ind1, ind2);
        hc_move = hc(ind1, ind2-1);
        hc_hat(va>=0) = hc_move(va>=0);
        Fcvh(ind1, ind2) = va.*hc_hat;
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

%%%%% Mass Balance 
        ind1 = 1 + Ib : Imax + Ib;
        ind2 = 1 + Ib : Jmax + Ib;
        hc(ind1, ind2) = hc(ind1, ind2) - dt * (Fcuh(ind1 + 1, ind2) - Fcuh(ind1, ind2)) / dx - dt * (Fcvh(ind1, ind2 + 1) - Fcvh(ind1, ind2)) / dy;


%%%%% Source Term - mass is added continuously at a rate
%%%%% equal the mass flux ''''''''''''''''''''''''''''''
     % concentration c(i,j) mass per unit volume
     % hc(i,j) is mass per unit area
     % mass flux is mass per unit time
        ind1 = 1 + Ib : 1 + EF + Ib;
        ind2 = 1 + Ib : Jmax / 2 + EF + Ib;
        hc(ind1, ind2) = hc(ind1, ind2) + massflux / 4 * dt / dx / EF / dy / EF;
%         for i = 1 + Ib : 1 + EF + Ib
%             for j = Jmax / 2 - EF - 1 + Ib : Jmax / 2 + EF + Ib;
%                         hc(i, j) = hc(i, j) + massflux / 4 * dt / dx / EF / dy / EF;
%             end
%         end
    end
   
    %%%%%save data and GIF
    clf
    %hold on
    drawnow
    u_plot = u(1:Imax,1:Jmax);
    v_plot = v(1:Imax,1:Jmax);
%     p1 = quiver(xcoord(1:selected_data:end), ycoord(1:selected_data:end), u_plot(1:selected_data:end), v_plot(1:selected_data:end));
%     set(p1, 'LineWidth', 1.0)
    % p1 = streamslice(xcoord,ycoord, v_plot, u_plot, 8);
%     figure(1)
%     clf
%     p1 = streamslice(xcoord,ycoord, v_plot, u_plot, 8);
%     aa =  get(gca,'XTickLabel');
%    % set(gca,'fontsize',16)
%     axis([0 max(xx) 0 1.2*max(yy)])
%     title(strcat('$t$=',num2str(t,3),'(s)'),'FontSize',14,'Interpreter','latex')
% %     text(1.1*max(xx), 1.1*max(yy),strcat('$t$=',num2str(t,3),'(s)'),'FontSize',14,'Interpreter','latex')
%     xlabel('\fontsize{12}{12}\selectfont {$x$(m)}','Interpreter','latex','FontWeight','bold')
%     ylabel('\fontsize{12}{12}\selectfont {$y$(m)}','Interpreter','latex','FontWeight','bold')
%     set(gca,'XTickLabel',aa,'fontsize',12)
%     grid on
%     clf;
%     saveas(gcf,strcat('./streamlines/',num2str(t),'.jpg'))
    
    clf
    axis([0 1.05*max(xx) 0 1.05*max(yy)])
    set(gcf,'Position',[100 100 1200 800])
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
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
    
    
%     F=getframe(gcf);
%     I=frame2im(F);
%     [I,map]=rgb2ind(I,256);
%     if pic_num == 1 %save the initial GIF 
%         imwrite(I,map,strcat('Case2','.gif'),'gif', 'Loopcount',inf,'DelayTime',0.2);
%     else %append additional frames
%         imwrite(I,map,strcat('Case2','.gif'),'gif','WriteMode','append','DelayTime',0.2);
%     end
%     pic_num = pic_num + 1;


    
%     if mod(ii,10) == 0
%         coord = [x',h'];
%         csvwrite(strcat(num2str(t),'.csv'),coord)
%     end
end

