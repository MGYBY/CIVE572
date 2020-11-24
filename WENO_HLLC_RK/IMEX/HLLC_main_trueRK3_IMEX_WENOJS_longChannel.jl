
using Plots, CSV, DataFrames, Printf
gr()
mkpath("./WSP1/")

function HLLC_source_term(h, q, tan_theta, cf)
    var_source_term = tan_theta * 9.81 .* h - (cf * q.^2) ./ (2 * h.^2)
    return var_source_term
end

function rhs_hyperbolic(qn, epsilon, period, t, g, Imax, h, q, hn)
    Fhllc = Array{Float64}(undef, 2, Imax)
    hstar = Array{Float64}(undef, 1, Imax)
    qqleft = Array{Float64}(undef, 1, Imax)
    qqright = Array{Float64}(undef, 1, Imax)
    SL = Array{Float64}(undef, 1, Imax)
    SR = Array{Float64}(undef, 1, Imax)
    Sstar = Array{Float64}(undef, 1, Imax)
    ULstar = Array{Float64}(undef, 2, Imax) # approximate states for 2 eqns
    URstar = Array{Float64}(undef, 2, Imax)

    # Zero Gradient (supersonic outlet)
    q[Imax + 1] = q[Imax]
    h[Imax + 1] = h[Imax]
    # instant dist.
    if t < period
        h[1] = hn + 0.2 * hn
        q[1] = qn + 0.2 * qn
    else
        h[1] = hn
        q[1] = qn
    end

    # left and right states
    u, UL, UR, FL, FR, aleft, aright = left_right_state(g, Imax, h, q)
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
        ULstar[:, k] = h[k] * ((SL[k] - u[k]) / (SL[k] - Sstar[k])) .* [1.0; Sstar[k]]
        URstar[:, k] = h[k + 1] * ((SR[k] - u[k + 1]) / (SR[k] - Sstar[k])) .* [1.0; Sstar[k]]
        # evaluate fluxes Fhllc for both cont. and momentum eqns
        # Fhllc(3/2...Imax+1/2)
        if SL[k] >= 0.0
            Fhllc[:, k] = FL[:, k]
        elseif (SL[k] < 0.0 && Sstar[k] > 0.0)
            Fhllc[:, k] = FL[:, k] + SL[k] * (ULstar[:, k] - UL[:, k])
        elseif (Sstar[k] < 0.0 && SR[k] >= 0.0)
            Fhllc[:, k] = FR[:, k] + SR[k] * (URstar[:, k] - UR[:, k])
        else
            Fhllc[:, k] = FR[:, k]
        end
    end
    return Fhllc
end

function left_right_state(g, Imax, h, q)
    # calc. left and right flux val.
    #  UL, UR, FL, FR, aleft, aright
    aleft = Array{Float64}(undef, 1, Imax) # geopotential following flux index
    aright = Array{Float64}(undef, 1, Imax) # geopotential following flux index
    UL = Array{Float64}(undef, 2, Imax)
    UR = Array{Float64}(undef, 2, Imax)
    FL = Array{Float64}(undef, 2, Imax) # FL and FR follow flux index
    FR = Array{Float64}(undef, 2, Imax)

    u = q ./ h
    # UL = [h[1:Imax]'; q[1:Imax]']
    # UR = [h[2:Imax+1]'; q[2:Imax+1]']
    UL = wenoL(Imax, h, q)
    UR = wenoR(Imax, h, q)
    fluxes(Imax, g, FL, UL[2, :], UL[1, :])
    fluxes(Imax, g, FR, UR[2, :], UR[1, :])
    aleft = sqrt.(g * UL[1, :])
    aright = sqrt.(g * UR[1, :])
    return u, UL, UR, FL, FR, aleft, aright
end

function fluxes(Imax, g, F, q, h)
    for ii = 1:Imax
        F[1, ii] = q[ii]
        F[2, ii] = ((q[ii]^2) / h[ii]) + 1 / 2 * g * (h[ii]^2)
    end
end

# left WENO reconstruction f(1...Imax+1)~f(3/2...Imax+3/2)
function wenoL(n, h, q)
    f = Array{Float64}(undef, 2, Imax)

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
    return f
end

function wenoR(n, h, q)
    f = Array{Float64}(undef, 2, Imax)

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
    f[2, i] = wcR(v1, v2, v3, v4, v5)

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

    c1 = nw_g(0.1, w1)
    c2 = nw_g(0.6, w2)
    c3 = nw_g(0.3, w3)

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

    c1 = nw_g(0.3, w1)
    c2 = nw_g(0.6, w2)
    c3 = nw_g(0.1, w3)

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

function nw_g(c, w)
    g = w * (c + c^2 - 3 * c * w + w^2) / (c^2 + w * (1 - 2 * c))
    return g
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
const epsilon = 0.01
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

const Imax = 2000 * 4
const xx = 160.0  # Brock's setup
const dx = xx / Imax
q = zeros(1, Imax + 1) .+ qn # or I.C. with zero discharge
# q = q+hn*un;
h = zeros(1, Imax + 1) .+  hn
qmed1 = zeros(1, Imax + 1)
hmed1 = zeros(1, Imax + 1)
qmed2 = zeros(1, Imax + 1)
hmed2 = zeros(1, Imax + 1)
qmed3 = zeros(1, Imax + 1)
hmed3 = zeros(1, Imax + 1)
qmed4 = zeros(1, Imax + 1)
hmed4 = zeros(1, Imax + 1)
# HLLC flux (for 2 eqns)
Fhllcp = Array{Float64}(undef, 2, Imax)
Fhllc1 = Array{Float64}(undef, 2, Imax)
Fhllc2 = Array{Float64}(undef, 2, Imax)
Fhllc3 = Array{Float64}(undef, 2, Imax)
Fhllc4 = Array{Float64}(undef, 2, Imax)
ks = Array{Float64}(undef, 1, Imax)
const sim_time = 200.0
const co = 0.20 # Cou
dt = co * dx / (1.0 * un)
t = 0.0

# OUTFLOW boundary conditions(supersonic)
q[Imax + 1] = q[Imax]
h[Imax + 1] = h[Imax]
# u[Imax + 1] = u[Imax]

x = 0.5 * dx:dx:(0.5 * dx + (Imax - 1) * dx) # computational domain: 1+Ib~Imax+Ib
pic_num = 1
tplot = 2.0 # clf %, drawnow%, set(gcf,'renderer','zbuffer')
plotgap = Int64(round(tplot / dt))
dt = tplot / plotgap
nplots = Int64(round(sim_time / tplot))
nt = Int64(round(sim_time / dt))
const txt_output = 10
const graphic_output = 2.0 
n_txt, n_graphic = 0.0, 0.0
const b_coeff = 1.0
const α, β, η = 0.24169426078821, 0.06042356519705, 0.12915286960590

# time loop
for i = 1:nt
    global t += dt
    #         Radiating B.C.
    #         q(Imax + 1) = sqrt(h(Imax)*g)*(h(Imax)-hn);

    # Zero Gradient (supersonic outlet)
    q[Imax + 1] = q[Imax]
    h[Imax + 1] = h[Imax]
    # instant dist.
    if t < period
        h[1] = hn + 0.2 * hn
        q[1] = qn + 0.2 * qn
    else
        h[1] = hn
        q[1] = qn
    end

    # Fhllcp .= rhs_hyperbolic(qn, epsilon, period, t, g, Imax, h, q, hn) # K1 and K1hat
    # solve the conserved form and obtain h and q (namely U)
    # implicit time integration SSP(4,3,3)
    for k = 2:Imax
        hmed1[k] = h[k]
        a_coeff = α * dt * cf / (hmed1[k]^2) / 2.0
        c_coeff = -(q[k] + 
        α * dt * (g * hmed1[k] * tan_theta))
        qmed1[k] = (-b_coeff + sqrt(b_coeff^2 - 4 * a_coeff * c_coeff)) / (2 * a_coeff)
    end

    Fhllc1 .= rhs_hyperbolic(qn, epsilon, period, t, g, Imax, hmed1, qmed1, hn) # K2 and K2hat
    # solve the conserved form and obtain h and q (namely U)
    # implicit time integration
    for k = 2:Imax
        hmed2[k] = h[k] - 0.0 * dt / dx * (Fhllc1[1, k] - Fhllc1[1, k - 1])
        a_coeff = α * dt * cf / (hmed2[k]^2) / 2.0
        c_coeff = -(q[k] - 0.0 * dt / dx * (Fhllc1[2, k] - Fhllc1[2, k - 1]) + 
        (α) * dt * (g * hmed2[k] * tan_theta) + 
        (-α) * dt * (HLLC_source_term(hmed1[k], qmed1[k], tan_theta, cf)))
        qmed2[k] = (-b_coeff + sqrt(b_coeff^2 - 4 * a_coeff * c_coeff)) / (2 * a_coeff)
    end

    Fhllc2 .= rhs_hyperbolic(qn, epsilon, period, t, g, Imax, hmed2, qmed2, hn) # K3 and K3hat
    # solve the conserved form and obtain h and q (namely U)
    # implicit time integration
    for k = 2:Imax
        hmed3[k] = h[k] - 1.0 * dt / dx * (Fhllc2[1, k] - Fhllc2[1, k - 1])
        a_coeff = α * dt * cf / (hmed3[k]^2) / 2.0
        c_coeff = -(q[k] - 1.0 * dt / dx * (Fhllc2[2, k] - Fhllc2[2, k - 1]) + 
        (1 - α) * dt * (HLLC_source_term(hmed2[k], qmed2[k], tan_theta, cf)) + 
        α * dt * (g * hmed3[k] * tan_theta))
        qmed3[k] = (-b_coeff + sqrt(b_coeff^2 - 4 * a_coeff * c_coeff)) / (2 * a_coeff)
    end

    Fhllc3 .= rhs_hyperbolic(qn, epsilon, period, t + 1.0 * dt, g, Imax, hmed3, qmed3, hn) # K4 and K4hat
    # solve the conserved form and obtain h and q (namely U)
    # implicit time integration
    for k = 2:Imax
        hmed4[k] = h[k] - (1.0 / 4.0) * dt / dx * (Fhllc2[1, k] - Fhllc2[1, k - 1]) - 
        (1.0 / 4.0) * dt / dx * (Fhllc3[1, k] - Fhllc3[1, k - 1])
        a_coeff = α * dt * cf / (hmed4[k]^2) / 2.0
        c_coeff = -(q[k] - (1.0 / 4.0) * dt / dx * (Fhllc2[2, k] - Fhllc2[2, k - 1]) - 
        (1.0 / 4.0) * dt / dx * (Fhllc3[2, k] - Fhllc3[2, k - 1]) +
        α * dt * (g * hmed4[k] * tan_theta) +
        (1.0 / 2.0 - β - η - α) * dt * (HLLC_source_term(hmed3[k], qmed3[k], tan_theta, cf)) +
        η * dt * (HLLC_source_term(hmed2[k], qmed2[k], tan_theta, cf)) +
        β * dt * (HLLC_source_term(hmed1[k], qmed1[k], tan_theta, cf)))
        qmed4[k] = (-b_coeff + sqrt(b_coeff^2 - 4 * a_coeff * c_coeff)) / (2 * a_coeff)
    end

    Fhllc4 .= rhs_hyperbolic(qn, epsilon, period, t + 0.5 * dt, g, Imax, hmed4, qmed4, hn) # K4 and K4hat
    for k = 2:Imax
        h[k] = h[k] - (1.0 / 6.0) * dt / dx * (Fhllc2[1, k] - Fhllc2[1, k - 1]) -
        (1.0 / 6.0) * dt / dx * (Fhllc3[1, k] - Fhllc3[1, k - 1]) -
        (2.0 / 3.0) * dt / dx * (Fhllc4[1, k] - Fhllc4[1, k - 1])

        q[k] = q[k] - (1.0 / 6.0) * dt / dx * (Fhllc2[2, k] - Fhllc2[2, k - 1]) -
        (1.0 / 6.0) * dt / dx * (Fhllc3[2, k] - Fhllc3[2, k - 1]) -
        (2.0 / 3.0) * dt / dx * (Fhllc4[2, k] - Fhllc4[2, k - 1]) +
        (1.0 / 6.0) * dt * (HLLC_source_term(hmed2[k], qmed2[k], tan_theta, cf)) +
        (1.0 / 6.0) * dt * (HLLC_source_term(hmed3[k], qmed3[k], tan_theta, cf)) +
        (2.0 / 3.0) * dt * (HLLC_source_term(hmed4[k], qmed4[k], tan_theta, cf))
    end

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
        png("./WSP1/" * @sprintf("%.2f", t) * ".png")
        if mod(t, txt_output) < dt
            CSV.write("./WSP1/" * @sprintf("%.2f", t) * ".csv", DataFrame([x h[1:Imax]]))
        end
    end
end
