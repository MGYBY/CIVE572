using Plots, CSV, DataFrames, Printf
gr()
mkpath("./WSP1/")

function HLLC_source_term(h, q, tan_theta, cf)
    var_source_term = tan_theta * 9.81 .* h - (cf * q.^2) ./ (2 * h.^2)
    return var_source_term
end

function rhs_hyperbolic(qn, epsilon, period, t, g, Imax, h, q, hc, hn)
    Fhllc = Array{Float64}(undef, 3, Imax)
    hstar = Array{Float64}(undef, 1, Imax)
    qqleft = Array{Float64}(undef, 1, Imax)
    qqright = Array{Float64}(undef, 1, Imax)
    SL = Array{Float64}(undef, 1, Imax)
    SR = Array{Float64}(undef, 1, Imax)
    Sstar = Array{Float64}(undef, 1, Imax)
    ULstar = Array{Float64}(undef, 3, Imax) # approximate states for 2 eqns
    URstar = Array{Float64}(undef, 3, Imax)

    # Zero Gradient (supersonic outlet)
    q[Imax + 1] = q[Imax]
    h[Imax + 1] = h[Imax]
    hc[Imax + 1] = hc[Imax]
    # inlet with dist.
    h[1] = hn * (1 + epsilon * sum(sin(2 * pi ./ period * t)))
    q[1] = qn
    hc[1] = 1.0 * hn

    # left and right states
    u, UL, UR, FL, FR, aleft, aright = left_right_state(g, Imax, h, q, hc)
    # solve the hyperbolic eqns (without source terms)
    for k = 1:Imax # at flux k
        # calc. estimate for h*
        hstar[k] =
            1 / 2 * (h[k + 1] + h[k]) -
            1 / 4 * (u[k + 1] - u[k]) * (h[k + 1] + h[k]) / (aleft[k] + aright[k])
        # calc. the estimated wave speed SL, SR, S*
        # left value
        if hstar[k] > h[k]
            qqleft[k] = sqrt(1 / 2 * ((hstar[k] + h[k]) * hstar[k]) / (h[k]^2))
        else
            qqleft[k] = 1.0
        end
        SL[k] = u[k] - aleft[k] * qqleft[k]
        # right value
        if hstar[k] > h[k + 1]
            qqright[k] = sqrt(1 / 2 * ((hstar[k] + h[k + 1]) * hstar[k]) / (h[k + 1]^2))
        else
            qqright[k] = 1.0
        end
        SR[k] = u[k + 1] - aright[k] * qqright[k]
        Sstar[k] =
            (SL[k] * h[k + 1] * (u[k + 1] - SR[k]) - SR[k] * h[k] * (u[k] - SL[k])) /
            (h[k + 1] * (u[k + 1] - SR[k]) - h[k] * (u[k] - SL[k]))

        # calc. approx states for U*L, U*R
        ULstar[:, k] = UL[1, k] * ((SL[k] - (UL[2, k] / UL[1, k]))) / ((SL[k] - Sstar[k])) .* [1.0; Sstar[k]; (UL[3, k] / UL[1, k])]
        URstar[:, k] = UR[1, k] * ((SR[k] - (UR[2, k] / UR[1, k]))) / ((SR[k] - Sstar[k])) .* [1.0; Sstar[k]; (UR[3, k] / UR[1, k])]
        # evaluate fluxes Fhllc for both cont. and momentum eqns
        # Fhllc(3/2...Imax+1/2)
        if SL[k] >= 0.0
            Fhllc[:, k] = FL[:, k]
        elseif (SL[k] < 0.0 && Sstar[k] > 0.0)
            Fhllc[:, k] = FL[:, k] .+ SL[k] .* (ULstar[:, k] - UL[:, k])
        elseif (Sstar[k] <= 0.0 && SR[k] > 0.0)
            Fhllc[:, k] = FR[:, k] .+ SR[k] .* (URstar[:, k] - UR[:, k])
        else
            Fhllc[:, k] = FR[:, k]
        end
    end
    return Fhllc
end

function left_right_state(g, Imax, h, q, hc)
    # calc. left and right flux val.
    #  UL, UR, FL, FR, aleft, aright
    aleft = Array{Float64}(undef, 1, Imax) # geopotential following flux index
    aright = Array{Float64}(undef, 1, Imax) # geopotential following flux index
    UL = Array{Float64}(undef, 3, Imax)
    UR = Array{Float64}(undef, 3, Imax)
    FL = Array{Float64}(undef, 3, Imax) # FL and FR follow flux index
    FR = Array{Float64}(undef, 3, Imax)

    u = q ./ h
    # UL = [h[1:Imax]'; q[1:Imax]']
    # UR = [h[2:Imax+1]'; q[2:Imax+1]']
    UL = wenoL(Imax, h, q, hc)
    UR = wenoR(Imax, h, q, hc)
    fluxes(Imax, g, FL, UL[2, :], UL[1, :], UL[3, :])
    fluxes(Imax, g, FR, UR[2, :], UR[1, :], UR[3, :])
    aleft = sqrt.((g * UL[1, :]))
    aright = sqrt.((g * UR[1, :]))
    return u, UL, UR, FL, FR, aleft, aright
end

function fluxes(Imax, g, F, q, h, hc)
    for ii = 1:Imax
        F[1, ii] = q[ii]
        F[2, ii] = ((q[ii]^2) / h[ii]) + 1.0 / 2.0 * g * (h[ii]^2)
        F[3, ii] = (q[ii] * hc[ii]) / (h[ii])
    end
end

# left WENO reconstruction f(1...Imax+1)~f(3/2...Imax+3/2)
function wenoL(n, h, q, hc)
    f = Array{Float64}(undef, 3, Imax)

    i = 1
    v1 = 3.0 * h[i] - 2.0 * h[i + 1]
    v2 = 2.0 * h[i] - h[i + 1]
    v3 = h[i]
    v4 = h[i + 1]
    v5 = h[i + 2]
    f[1, i] = wcL(v1, v2, v3, v4, v5)
    v1 = 3.0 * q[i] - 2.0 * q[i + 1]
    v2 = 2.0 * q[i] - q[i + 1]
    v3 = q[i]
    v4 = q[i + 1]
    v5 = q[i + 2]
    f[2, i] = wcL(v1, v2, v3, v4, v5)
    v1 = 3.0 * hc[i] - 2.0 * hc[i + 1]
    v2 = 2.0 * hc[i] - hc[i + 1]
    v3 = hc[i]
    v4 = hc[i + 1]
    v5 = hc[i + 2]
    f[3, i] = wcL(v1, v2, v3, v4, v5)

    i = 2
    v1 = 2.0 * h[i - 1] - h[i]
    v2 = h[i - 1]
    v3 = h[i]
    v4 = h[i + 1]
    v5 = h[i + 2]
    f[1, i] = wcL(v1, v2, v3, v4, v5)
    v1 = 2.0 * q[i - 1] - q[i]
    v2 = q[i - 1]
    v3 = q[i]
    v4 = q[i + 1]
    v5 = q[i + 2]
    f[2, i] = wcL(v1, v2, v3, v4, v5)
    v1 = 2.0 * hc[i - 1] - hc[i]
    v2 = hc[i - 1]
    v3 = hc[i]
    v4 = hc[i + 1]
    v5 = hc[i + 2]
    f[3, i] = wcL(v1, v2, v3, v4, v5)

    for i = 3:n - 1
        v1 = h[i - 2]
        v2 = h[i - 1]
        v3 = h[i]
        v4 = h[i + 1]
        v5 = h[i + 2]
        f[1, i] = wcL(v1, v2, v3, v4, v5)
        v1 = q[i - 2]
        v2 = q[i - 1]
        v3 = q[i]
        v4 = q[i + 1]
        v5 = q[i + 2]
        f[2, i] = wcL(v1, v2, v3, v4, v5)
        v1 = hc[i - 2]
        v2 = hc[i - 1]
        v3 = hc[i]
        v4 = hc[i + 1]
        v5 = hc[i + 2]
        f[3, i] = wcL(v1, v2, v3, v4, v5)
    end

    i = n
    v1 = h[i - 2]
    v2 = h[i - 1]
    v3 = h[i]
    v4 = h[i + 1]
    v5 = 2.0 * h[i + 1] - h[i]
    f[1, i] = wcL(v1, v2, v3, v4, v5)
    v1 = q[i - 2]
    v2 = q[i - 1]
    v3 = q[i]
    v4 = q[i + 1]
    v5 = 2.0 * q[i + 1] - q[i]
    f[2, i] = wcL(v1, v2, v3, v4, v5)
    v1 = hc[i - 2]
    v2 = hc[i - 1]
    v3 = hc[i]
    v4 = hc[i + 1]
    v5 = 2.0 * hc[i + 1] - hc[i]
    f[3, i] = wcL(v1, v2, v3, v4, v5)
    return f
end

function wenoR(n, h, q, hc)
    f = Array{Float64}(undef, 3, Imax)

    i = 1
    v1 = 2.0 * h[i] - h[i + 1]
    v2 = h[i]
    v3 = h[i + 1]
    v4 = h[i + 2]
    v5 = h[i + 3]
    f[1, i] = wcR(v1, v2, v3, v4, v5)
    v1 = 2.0 * q[i] - q[i + 1]
    v2 = q[i]
    v3 = q[i + 1]
    v4 = q[i + 2]
    v5 = q[i + 3]
    f[2, i] = wcR(v1, v2, v3, v4, v5)
    v1 = 2.0 * hc[i] - hc[i + 1]
    v2 = hc[i]
    v3 = hc[i + 1]
    v4 = hc[i + 2]
    v5 = hc[i + 3]
    f[3, i] = wcR(v1, v2, v3, v4, v5)


    for i = 2:n - 2
        v1 = h[i - 1]
        v2 = h[i]
        v3 = h[i + 1]
        v4 = h[i + 2]
        v5 = h[i + 3]
        f[1, i] = wcR(v1, v2, v3, v4, v5)
        v1 = q[i - 1]
        v2 = q[i]
        v3 = q[i + 1]
        v4 = q[i + 2]
        v5 = q[i + 3]
        f[2, i] = wcR(v1, v2, v3, v4, v5)
        v1 = hc[i - 1]
        v2 = hc[i]
        v3 = hc[i + 1]
        v4 = hc[i + 2]
        v5 = hc[i + 3]
        f[3, i] = wcR(v1, v2, v3, v4, v5)
    end

    i = n - 1
    v1 = h[i - 1]
    v2 = h[i]
    v3 = h[i + 1]
    v4 = h[i + 2]
    v5 = 2.0 * h[i + 2] - h[i + 1]
    f[1, i] = wcR(v1, v2, v3, v4, v5)
    v1 = q[i - 1]
    v2 = q[i]
    v3 = q[i + 1]
    v4 = q[i + 2]
    v5 = 2.0 * q[i + 2] - q[i + 1]    
    v1 = hc[i - 1]
    v2 = hc[i]
    v3 = hc[i + 1]
    v4 = hc[i + 2]
    v5 = 2.0 * hc[i + 2] - hc[i + 1]
    f[3, i] = wcR(v1, v2, v3, v4, v5)


    i = n
    v1 = h[i - 1]
    v2 = h[i]
    v3 = h[i + 1]
    v4 = 2.0 * h[i + 1] - h[i]
    v5 = 3.0 * h[i + 1] - 2.0 * h[i]
    f[1, i] = wcR(v1, v2, v3, v4, v5)
    v1 = q[i - 1]
    v2 = q[i]
    v3 = q[i + 1]
    v4 = 2.0 * q[i + 1] - q[i]
    v5 = 3.0 * q[i + 1] - 2.0 * q[i]
    f[2, i] = wcR(v1, v2, v3, v4, v5)
    v1 = hc[i - 1]
    v2 = hc[i]
    v3 = hc[i + 1]
    v4 = 2.0 * hc[i + 1] - hc[i]
    v5 = 3.0 * hc[i + 1] - 2.0 * hc[i]
    f[3, i] = wcR(v1, v2, v3, v4, v5)

    return f
end

function wcL(v1, v2, v3, v4, v5)
    eps = 1.0e-6

    # smoothness indicators
    s1 = (13.0 / 12.0) * (v1 - 2.0 * v2 + v3)^2 + 0.25 * (v1 - 4.0 * v2 + 3.0 * v3)^2
    s2 = (13.0 / 12.0) * (v2 - 2.0 * v3 + v4)^2 + 0.25 * (v2 - v4)^2
    s3 = (13.0 / 12.0) * (v3 - 2.0 * v4 + v5)^2 + 0.25 * (3.0 * v3 - 4.0 * v4 + v5)^2

    # computing nonlinear weights w1,w2,w3
    c1 = 1.0e-1 / ((eps + s1)^2)
    c2 = 6.0e-1 / ((eps + s2)^2)
    c3 = 3.0e-1 / ((eps + s3)^2)

    w1 = c1 / (c1 + c2 + c3)
    w2 = c2 / (c1 + c2 + c3)
    w3 = c3 / (c1 + c2 + c3)

    # candiate stencils
    q1 = v1 / 3.0 - 7.0 / 6.0 * v2 + 11.0 / 6.0 * v3
    q2 = -v2 / 6.0 + 5.0 / 6.0 * v3 + v4 / 3.0
    q3 = v3 / 3.0 + 5.0 / 6.0 * v4 - v5 / 6.0

    # reconstructed value at interface
    f = (w1 * q1 + w2 * q2 + w3 * q3)

    return f
end

function wcR(v1, v2, v3, v4, v5)
    eps = 1.0e-6

    s1 = (13.0 / 12.0) * (v1 - 2.0 * v2 + v3)^2 + 0.25 * (v1 - 4.0 * v2 + 3.0 * v3)^2
    s2 = (13.0 / 12.0) * (v2 - 2.0 * v3 + v4)^2 + 0.25 * (v2 - v4)^2
    s3 = (13.0 / 12.0) * (v3 - 2.0 * v4 + v5)^2 + 0.25 * (3.0 * v3 - 4.0 * v4 + v5)^2

    c1 = 3.0e-1 / (eps + s1)^2
    c2 = 6.0e-1 / (eps + s2)^2
    c3 = 1.0e-1 / (eps + s3)^2

    w1 = c1 / (c1 + c2 + c3)
    w2 = c2 / (c1 + c2 + c3)
    w3 = c3 / (c1 + c2 + c3)

    # candiate stencils
    q1 = -v1 / 6.0 + 5.0 / 6.0 * v2 + v3 / 3.0
    q2 = v2 / 3.0 + 5.0 / 6.0 * v3 - v4 / 6.0
    q3 = 11.0 / 6.0 * v3 - 7.0 / 6.0 * v4 + v5 / 3.0

    # reconstructed value at interface
    f = (w1 * q1 + w2 * q2 + w3 * q3)

    return f
end

# include("./HLLC_source_term.jl")

# % %stable case Fr=1.15
# % hn = 0.0798;
# % un = 0.32176;
# % g = 9.810;
# % sine_theta = 0.05;
# % tan_theta = sine_theta/sqrt(1-sine_theta^2);
# % cf = g*tan_theta*2*hn/(un^2);
# % epsilon =  0.1;
# % period = 0.933; %large amplitude but still damped
# % % period = [0.933, 1.25];
# % % period = 0.9:0.025:1.1;
# % % period = randn(1,5);

# non-staggered grid
# unstable case1 Fr=3.7
const hn = 0.0079756
const un = 1.03774
const qn = hn * un
const g = 9.810
const sine_theta = 0.05011
const tan_theta = sine_theta / sqrt(1 - sine_theta^2)
const cf = g * sine_theta * 2 * hn / (un^2)
const epsilon = 0.05
const period = 0.933

#  # unstable case2 Fr=5.6
#  hn = 0.00533;
#  un = 1.2805;
#  cf = 9.85e-3;
#  g = 9.810;
#  sine_theta = 0.119528;
#  tan_theta = sine_theta/sqrt(1-sine_theta^2);
#  epsilon =  0.01; period = 1.015;

#  # Coarsening Effect Imax=40000, xx=400
#  hn = 0.0079756;
#  un = 1.03774;
#  cf = 7.2813e-3;
#  g = 9.810;
#  sine_theta = 0.05011;
#  tan_theta = sine_theta/sqrt(1-sine_theta^2);
#  epsilon =  0.01; %period = 0.933;
#  period = 0.6:0.3:1.8;

const Imax = 2500
const xx = 40.0  # Brock's setup
const dx = xx / Imax
q = zeros(1, Imax + 1) .+ qn # or I.C. with zero discharge
# q = q+hn*un;
h = zeros(1, Imax + 1) .+ hn
hc = zeros(1, Imax + 1) .+ hn * 1.0 # concentration-color eqn
# c = zeros(1, Imax + 1)
qmed = zeros(1, Imax + 1)
hmed = zeros(1, Imax + 1)
hcmed = zeros(1, Imax + 1)
# HLLC flux (for 2 eqns)
Fhllc = Array{Float64}(undef, 3, Imax)
const sim_time = 50.0
const co = 0.40
dt = co * dx / (1.0 * un)
t = 0.0

# OUTFLOW boundary conditions(supersonic)
q[Imax + 1] = q[Imax]
h[Imax + 1] = h[Imax]
hc[Imax + 1] = hc[Imax]
# u[Imax + 1] = u[Imax]

# peak information arrays
h_global_max = zero(h)
h_local_max = zero(h)
peak_height = zero(h)
peak_loc = zero(h)

x = 0.5 * dx:dx:(0.5 * dx + (Imax - 1) * dx) # computational domain: 1+Ib~Imax+Ib
pic_num = 1
tplot = 2.0 # clf %, drawnow%, set(gcf,'renderer','zbuffer')
plotgap = Int64(round(tplot / dt))
dt = tplot / plotgap
nplots = Int64(round(sim_time / tplot))
nt = Int64(round(sim_time / dt))
const txt_output = 10
const graphic_output = 2
n_txt, n_graphic = 0.0, 0.0
const peak_output = 1.0
peak_count = 0

# time loop
for i = 1:nt
    global t += dt
    #         Radiating B.C.
    #         q(Imax + 1) = sqrt(h(Imax)*g)*(h(Imax)-hn);

    # Zero Gradient (supersonic outlet)
    q[Imax + 1] = q[Imax]
    h[Imax + 1] = h[Imax]
    # inlet with dist.
    h[1] = hn * (1 + epsilon * sum(sin(2 * pi ./ period * t)))
    q[1] = qn

    Fhllc .= rhs_hyperbolic(qn, epsilon, period, t, g, Imax, h, q, hc, hn)
    # solve the conserved form and obtain h and q (namely U)
    # hyperbolic problem
    # RK3
    for k = 2:Imax
        hmed[k] = h[k] - dt / dx * (Fhllc[1, k] - Fhllc[1, k - 1])
        qmed[k] = q[k] - dt / dx * (Fhllc[2, k] - Fhllc[2, k - 1])
        hcmed[k] = hc[k] - dt / dx * (Fhllc[3, k] - Fhllc[3, k - 1])
    end

    Fhllc .= rhs_hyperbolic(qn, epsilon, period, t, g, Imax, hmed, qmed, hcmed, hn)
    for k = 2:Imax
        hmed[k] =
            3.0 / 4.0 * h[k] +
            1.0 / 4.0 * (hmed[k] + (-dt / dx) * (Fhllc[1, k] - Fhllc[1, k - 1]))
        qmed[k] =
            3.0 / 4.0 * q[k] +
            1.0 / 4.0 * (qmed[k] + (-dt / dx) * (Fhllc[2, k] - Fhllc[2, k - 1]))
        hcmed[k] =
            3.0 / 4.0 * hc[k] +
            1.0 / 4.0 * (hcmed[k] + (-dt / dx) * (Fhllc[3, k] - Fhllc[3, k - 1]))
    end

    Fhllc .= rhs_hyperbolic(qn, epsilon, period, t, g, Imax, hmed, qmed, hcmed, hn)
    for k = 2:Imax
        h[k] =
            1.0 / 3.0 * h[k] +
            2.0 / 3.0 * (hmed[k] + (-dt / dx) * (Fhllc[1, k] - Fhllc[1, k - 1]))
        q[k] =
            1.0 / 3.0 * q[k] +
            2.0 / 3.0 * (qmed[k] + (-dt / dx) * (Fhllc[2, k] - Fhllc[2, k - 1]))
        hc[k] =
            1.0 / 3.0 * hc[k] +
            2.0 / 3.0 * (hcmed[k] + (-dt / dx) * (Fhllc[3, k] - Fhllc[3, k - 1]))
    end

    # solve the ODE for time integration (RK3 or RK4)
    # bed-friction and concentration source (if needed)
    for k = 2:Imax
        qmed[k] = q[k] + dt * (HLLC_source_term(h[k], q[k], tan_theta, cf))
    end
    for k = 2:Imax
        qmed[k] =
            3.0 / 4.0 * q[k] +
            1.0 / 4.0 * (qmed[k] + dt * HLLC_source_term(h[k], qmed[k], tan_theta, cf))
    end
    for k = 2:Imax
        q[k] =
            1.0 / 3.0 * q[k] +
            2.0 / 3.0 * (qmed[k] + dt * HLLC_source_term(h[k], qmed[k], tan_theta, cf))
    end

    u = q ./ h
    c = hc ./ h

    if mod(t, graphic_output) < dt
        # save data and GIF
        Plots.plot(
            x[1:Imax],
            h[1:Imax],
            xlabel="x (m)",
            ylabel="h (m)",
            lw=2,
            xlims=(0, maximum(x)),
            ylims=(0, maximum(h) * 1.10),
        )
        png("./WSP1/" * @sprintf("%.0f", t) * ".png")
        Plots.plot(
            x[1:Imax],
            c[1:Imax],
            xlabel="x (m)",
            ylabel="c",
            lw=2,
            color=[:black],
            xlims=(0, maximum(x)),
            ylims=(0, maximum(c) * 1.10),
        )
        png("./WSP1/" * @sprintf("c_%.0f", t) * ".png")
        KE = 1.0 / 2.0 * (h .* u.^2)
        Fr = u ./ sqrt.(g * h)
        CSV.write("./WSP1/" * @sprintf("%.0f", t) * ".csv", DataFrame([x h[1:Imax] u[1:Imax] KE[1:Imax] Fr[1:Imax] c[1:Imax]]))
    end

        # # global and local maximums
    # # local max
    # if t > 33
    #     for ii = 1:Imax
    #         if h[ii] > h_local_max[ii]
    #             h_local_max[ii] = h[ii]
    #         end
    #     end
    #     CSV.write("./WSP1/local_max.csv", DataFrame([x h_local_max[1:Imax]]))
    # end

    if mod(t, peak_output) < dt
        # leading wave tracking
        h_max = maximum(h[1:Imax])
        x_hmax = (findmax(h)[2][2] - 1.0 / 2.0) * dx
        CSV.write("./WSP1/max_height.csv", DataFrame([t x_hmax h_max]), append=true)
        
        # other peaks tracking
        for k = (Imax - 2):-1:3
            if (h[k - 2] < h[k - 1]) && (h[k - 1] < h[k]) && (h[k] > h[k + 1]) && (h[k + 1] > h[k + 2]) && (h[k] > 1.05 * hn) # strong peak identifier & filter
                global peak_count += 1
                peak_height[peak_count] = h[k]
                peak_loc[peak_count] = (k - 1.0 / 2.0) * dx
            end
        end
        if peak_count > 0
            peak_solution = open("./WSP1/peak_info.txt", "a")
            write(peak_solution, string(t), ",")
            for k = 1:peak_count
                write(peak_solution, string(peak_loc[k]), ",", string(peak_height[k]), ",")
            end
            write(peak_solution, "\n", )
            close(peak_solution)
        end
        global peak_count = 0
    end
end
