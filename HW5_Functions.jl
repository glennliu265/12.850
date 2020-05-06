using Plots
using LinearAlgebra
using Printf
include("AllFunctions12850.jl")

#immpt
dx      = 1
dt      = 1
u_value = 1

# Time
tsmax = 100

# Make Grid
xmin      = 0
xmax      = 100
ugrid  = [xmin:dx:xmax;]
tgrid  = ocnmod.get_midpoints(ugrid)
xmax   = length(tgrid)


#
q1t = 2
# Prescribe tracer field
if q1t == 1
    c0      = ones(size(tgrid)) * 5
    mask    = 20 .< tgrid .< 50
    c0 = c0 .* mask
elseif q1t == 2
    c0 = 5*sin.(2*pi.* tgrid ./ tgrid[xmax])
end
plot(tgrid,c0)


# Run FEM_UW1
#c50,tgrid50,cour=ocnmod.FEM_UW1(dx,dt,u_value,tsmax,c0,xmin,xmax)


animtest = @animate for t ∈ 1:tsmax
    plot(tgrid75,c75[:,t],
        title ="FEM-UW1, v = "*@sprintf("%.2f",cour)*"; t="*string(t),
        ylims=(0,6))

end
gif(animtest,"HW5_test_uw1_type"*string(q1t)*".gif",fps=5)


## LFM-CD

c50,tgrid50,cour=ocnmod.LFM_CD(dx,dt,u_value,tsmax,c0,xmin,xmax)






## Make Plot For FEM-UW1-Sine
l = @layout[a b]
p1 = plot(tgrid,c0,
            label="I.C.",
            lc=:black,
            lw=2.5,
            xlabel="X",
            ylabel="Amplitude",
            ylims=(0,10),
            title="Tracer @t=50 (FEM-UW1)",
            dpi=600
            )
p1 = plot!(tgrid25,c25[:,50],
            label="v = 0.25",
            lc=:blue,
            lw=2.5
            )
p1 = plot!(tgrid50,c50[:,50],
            label="v = 0.50",
            lc=:orange,
            lw=2.5,
            )

p1 = plot!(tgrid,c100[:,50],
            label="v = 1.00",
            lc=:black,
            ls=:dot,
            lw=2.5,
            )

p1 = hline!([5],
            lc=:black,
            lw=1,
            ls=:dash,
            label="")

# Plot time evolution
p2 = plot(tgrid,c0,
            label="I.C.",
            lc=:black,
            lw=2.5,
            xlabel="X",
            ylabel="Amplitude",
            ylims=(0,10),
            title="Tracer Adv; v = 0.5",
            dpi=600
            )
p2 = plot!(tgrid50,c50[:,1],
            label="t = 1",
            lc=RGB(250/255,217/255,170/255),
            lw=2.5,
            )
p2 = plot!(tgrid50,c50[:,25],
            label="t = 25",
            lc=RGB(248/255,199/255,127/255),
            lw=2.5,
            )
p2 = plot!(tgrid50,c50[:,50],
            label="t = 50",
            lc=RGB(242/255,168/255,59/255),
            lw=2.5,
            )
p2 = hline!([5],
            lc=:black,
            lw=1,
            ls=:dash,
            label="")
pl1 = plot(p2,p1,layout=l)
savefig(pl1,"HW5_FEMUW1_SinWave.png")
JLD.save("/Users/gyl/HW5_FEMUW1_SinWavem.jld","c25",c25,"c50",c50,"c100",c100,"tgrid25",tgrid25,"tgrid50",tgrid50,"tgrid100",tgrid100)

## Make Plot For FEM-UW1-Square
l = @layout[a b]
p1 = plot(tgrid,c0,
            label="I.C.",
            lc=:black,
            lw=2.5,
            xlabel="X",
            ylabel="Amplitude",
            ylims=(0,10),
            title="Tracer @t=50 (FEM-UW1)",
            dpi=600
            )
p1 = plot!(tgrid25,c25[:,50],
            label="v = 0.25",
            lc=:blue,
            lw=2.5
            )
p1 = plot!(tgrid50,c50[:,50],
            label="v = 0.50",
            lc=:orange,
            lw=2.5,
            )
p1 = plot!(tgrid75,c75[:,50],
            label="v = 0.75",
            lc=:green,
            lw=2.5,
            )
p1 = plot!(tgrid,c100[:,50],
            label="v = 1.00",
            lc=:black,
            ls=:dot,
            lw=2.5,
            )

p1 = hline!([5],
            lc=:black,
            lw=1,
            ls=:dash,
            label="")

# Plot time evolution
p2 = plot(tgrid,c0,
            label="I.C.",
            lc=:black,
            lw=2.5,
            xlabel="X",
            ylabel="Amplitude",
            ylims=(0,10),
            title="Tracer Adv; v = 0.5",
            dpi=600
            )
p2 = plot!(tgrid50,c50[:,1],
            label="t = 1",
            lc=RGB(250/255,217/255,170/255),
            lw=2.5,
            )
p2 = plot!(tgrid50,c50[:,25],
            label="t = 25",
            lc=RGB(248/255,199/255,127/255),
            lw=2.5,
            )
p2 = plot!(tgrid50,c50[:,50],
            label="t = 50",
            lc=RGB(242/255,168/255,59/255),
            lw=2.5,
            )
p2 = hline!([5],
            lc=:black,
            lw=1,
            ls=:dash,
            label="")
pl1 = plot(p2,p1,layout=l)
savefig(pl1,"HW5_FEMUW1_SqWave.png")
JLD.save("/Users/gyl/HW5_FEMUW1_SqWavem.jld","c25",c25,"c50",c50,"c75",c75,"c100",c100,"tgrid25",tgrid25,"tgrid50",tgrid50,"tgrid75",tgrid75,"tgrid100",tgrid100)

##
## Make Plot For LFM_CD-Square
l = @layout[a b]
p1 = plot(tgrid,c0,
            label="I.C.",
            lc=:black,
            lw=2.5,
            xlabel="X",
            ylabel="Amplitude",
            ylims=(0,10),
            title="Tracer @t=50 (LFM-CD)",
            dpi=600
            )
p1 = plot!(tgrid25,c25[:,50],
            label="v = 0.25",
            lc=:blue,
            lw=2.5
            )
p1 = plot!(tgrid50,c50[:,50],
            label="v = 0.50",
            lc=:orange,
            lw=2.5,
            )
p1 = plot!(tgrid,c100[:,50],
            label="v = 1.00",
            lc=:black,
            ls=:dot,
            lw=2.5,
            )

p1 = hline!([5],
            lc=:black,
            lw=1,
            ls=:dash,
            label="")

# Plot time evolution
p2 = plot(tgrid,c0,
            label="I.C.",
            lc=:black,
            lw=2.5,
            xlabel="X",
            ylabel="Amplitude",
            ylims=(0,10),
            title="Tracer Adv; v = 0.5",
            dpi=600
            )
p2 = plot!(tgrid50,c50[:,1],
            label="t = 1",
            lc=RGB(250/255,217/255,170/255),
            lw=2.5,
            )
p2 = plot!(tgrid50,c50[:,25],
            label="t = 25",
            lc=RGB(248/255,199/255,127/255),
            lw=2.5,
            )
p2 = plot!(tgrid50,c50[:,50],
            label="t = 50",
            lc=RGB(242/255,168/255,59/255),
            lw=2.5,
            )
p2 = hline!([5],
            lc=:black,
            lw=1,
            ls=:dash,
            label="")
pl1 = plot(p2,p1,layout=l)
savefig(pl1,"HW5_LFMCD_SqWave.png")
JLD.save("/Users/gyl/HW5_LFMCD_SqWavem.jld","c25",c25,"c50",c50,"c100",c100,"tgrid25",tgrid25,"tgrid50",tgrid50,"tgrid75",tgrid75,"tgrid100",tgrid100)



## Make Plot For LFM-CD-Sine
l = @layout[a b]
p1 = plot(tgrid,c0,
            label="I.C.",
            lc=:black,
            lw=2.5,
            xlabel="X",
            ylabel="Amplitude",
            ylims=(0,10),
            title="Tracer @t=50 (LFM-CD)",
            dpi=600
            )
p1 = plot!(tgrid25,c25[:,50],
            label="v = 0.25",
            lc=:blue,
            lw=2.5
            )
p1 = plot!(tgrid50,c50[:,50],
            label="v = 0.50",
            lc=:orange,
            lw=2.5,
            )

p1 = plot!(tgrid,c100[:,50],
            label="v = 1.00",
            lc=:black,
            ls=:dot,
            lw=2.5,
            )

p1 = hline!([5],
            lc=:black,
            lw=1,
            ls=:dash,
            label="")

# Plot time evolution
p2 = plot(tgrid,c0,
            label="I.C.",
            lc=:black,
            lw=2.5,
            xlabel="X",
            ylabel="Amplitude",
            ylims=(0,10),
            title="Tracer Adv; v = 0.5",
            dpi=600
            )
p2 = plot!(tgrid50,c50[:,1],
            label="t = 1",
            lc=RGB(250/255,217/255,170/255),
            lw=2.5,
            )
p2 = plot!(tgrid50,c50[:,25],
            label="t = 25",
            lc=RGB(248/255,199/255,127/255),
            lw=2.5,
            )
p2 = plot!(tgrid50,c50[:,50],
            label="t = 50",
            lc=RGB(242/255,168/255,59/255),
            lw=2.5,
            )
p2 = hline!([5],
            lc=:black,
            lw=1,
            ls=:dash,
            label="")
pl1 = plot(p2,p1,layout=l)
savefig(pl1,"HW5_LFM-CD_SinWave.png")
JLD.save("/Users/gyl/HW5_LFM-CD_SinWavem.jld","c25",c25,"c50",c50,"c100",c100,"tgrid25",tgrid25,"tgrid50",tgrid50,"tgrid100",tgrid100)


## FEM - UW2
# Evaluate stability
#cour = u[1] * dt / dx

# FEM and UW2
c = zeros(Float64,xmax,tsmax)
for t = 1:tsmax
    dt = 1

    if t == 1
        cin = c0
    else
        cin = c[:,t-1]
    end

    ct = zeros(Float64,size(tgrid))
    for i = 1:xmax
        ui = u[i]
        ci = cin

        # Assume Periodic
        if i == 1
            ip1 = i+1
            ip2 = i+2
            im1 = xmax
            im2 = xmax-1

        elseif i == 2

            ip1 = i+1
            ip2 = i+2
            im1 = i-1
            im2 = xmax

        elseif i == xmax

            ip1 = 1
            ip2 = 2
            im1 = i-1
            im2 = i-2

        elseif i == xmax-1

            ip1 = i+1
            ip2 = 1
            im1 = i-1
            im2 = i-2

        else

            ip1 = i+1
            ip2 = i+2
            im1 = i-1
            im2 = i-2

        end

        # Calculate u+ and u-
        up = (ui + abs(ui))/2
        um = (ui - abs(ui))/2



        ct[i] = ci[i] - dt/(dx*2) * (up* (3*ci[i] - 4*ci[im1] + 1*ci[im2])
                + um*(-1*ci[ip2]+4*ci[ip1]-3*ci[i]) )

    end

    c[:,t] = ct
end


animtest = @animate for t ∈ 1:tsmax
    plot(tgrid,c[:,t],
        title ="t="*string(t),
        ylims=(0,6)
        )

end
gif(animtest,"HW5_test_uw1_type"*string(q1t)*".gif",fps=5)





# Amplification Factor Modulus
# FEM-UW1
# Courant Number
nu = u[1]*dt/dx

# Create k vector
kdx = [0:(pi/12):pi;]


a2 = (1 .- nu .+ nu .* cos.(kdx)).^2
b2 = (nu .* sin.(kdx)).^2
lambda = sqrt.(a2.+b2)

Plots.plot(kdx,lambda,proj=:polar,lims=(0,1))
