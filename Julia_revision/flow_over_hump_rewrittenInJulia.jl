using Plots
using Printf
pyplot()

const xmax = 1000.0
const Ho = 10.0 # initial detph
const g = 9.81
const n = 0.014
const weirheight = 0.05 * Ho
t = 0.0
const dx = 10.0
c = sqrt(g * Ho)
Courant = 0.40
dt = Courant * dx / c
Imax = Int64(round(xmax / dx))
const sim_time = 1500.0

Fr = Array{Float64}(undef, 1, Imax)
h = Array{Float64}(undef, 1, Imax)
u = Array{Float64}(undef, 1, Imax + 1)
q = Array{Float64}(undef, 1, Imax + 1)
uhat = Array{Float64}(undef, 1, Imax)
hhat = Array{Float64}(undef, 1, Imax)
Fuuh = Array{Float64}(undef, 1, Imax)
# bottom elevation initial zero
zo = zeros(1, Imax)
# zo = Array{Float64}(undef,1, Imax)


# bottom elevation of the dam - sinusoidal from istart to (istart+ilength)
istart = Int64(Imax / 5.0)
ilength = Int64(Imax / 5.0)
zo[istart:istart + ilength] = 8 * sin.((0:(ilength)) / ilength * pi)

# Initial discharge is zero. Initial depth plus bottom elevation is constant.
h[1:Imax] = Ho .- zo[1:Imax]
q[1:Imax] .= 0.0

x = 0.5 * dx:dx:(0.5 * dx + (Imax - 1) * dx)
pic_num = 1
tplot = 1 # clf %, drawnow%, set(gcf,'renderer','zbuffer')
plotgap = round(tplot / dt)
dt = tplot / plotgap
nplots = round(sim_time / tplot)

# time loop
for ii = 1:nplots
    for nn = 1:plotgap
        global t += dt
        # INFLOW boundary conditions
        q[1] = 10.0
        h[1] = Ho
        u[1] = q[1] / Ho
        # OUTFLOW boundary conditions (subsonic outlet)
        # change the weir height from 0.10*Ho to 0.70*Ho to see
        # how the downstrem depth may alter the flow profile
        q[Imax + 1] = sqrt(g * h[Imax]) * (h[Imax] - weirheight)
        u[Imax + 1] = sqrt(g / h[Imax]) * (h[Imax] - weirheight)

        # continuity equation
        h[2:Imax] = h[2:Imax] - dt / dx * (q[3:Imax + 1] - q[2:Imax])

        # momentum equation
        hhat[2:Imax] = (h[2:Imax] + h[1:Imax - 1]) / 2.0;
        # gain computational stability by division of an averaged depth
        u[2:Imax] = q[2:Imax] ./ hhat[2:Imax]
        Fr[2:Imax] = u[2:Imax] ./ sqrt.(g * hhat[2:Imax])
        for k = 2:Imax
        # upwind scheme
            if ((u[k] + u[k + 1]) > 0)
                uhat[k] = u[k]
            else
                uhat[k] = u[k + 1]
            end
        # central difference
        # uhat(i) = (u(i) + u(i + 1)) / 2
        # momentum flux on the plus side of the finite volume
            Fuuh[k] = uhat[k]^2 * h[k];
        end
        Fuuh[1] = h[1] * u[1]^2;

        frictionterm = g * n^2 * q[2:Imax].^2.0 ./ (hhat[2:Imax].^(7.0 / 3.0))
        flux = -g * hhat[2:Imax] .* ((h[2:Imax] + zo[2:Imax]) - (h[1:Imax - 1] + zo[1:Imax - 1])) + (Fuuh[1:Imax - 1] - Fuuh[2:Imax])
        q[2:Imax] = q[2:Imax] + flux * dt / dx - frictionterm * dt
    end

    # plot WSP and Fr
    plot(x[1:Imax], zo[1:Imax], lw=3)
    plot!(x[1:Imax], zo[1:Imax].+h[1:Imax], label="H", lw=1.75)
    plot!(x[1:Imax], Fr[1:Imax], label="Fr", lw=1.75)
    title!("t=$t")
    png("./WSP/$t.png")
    # save WSP
    field_final = open("./WSP/$t.txt", "w")
    for k = 1:Imax
        write(field_final, @sprintf("%.16f",x[k]), " ", @sprintf("%.16f", h[k]), " \n")
    end
    close(field_final)
   #= save data and GIF
    # clf
    # hold on
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

    grid on =#
end
