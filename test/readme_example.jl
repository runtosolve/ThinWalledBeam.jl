using ThinWalledBeam


#Example 2 from Plaut, R.H., Moen, C.D.(2021). "Lateral-Torsional Deformations of Single-Span and Two-Span Thin-Walled Beams with Continuous Bracing".
#Journal of Constructional Steel Research.

#Z-section, no slope, single span, fixed-fixed
#kϕ=1500 N*mm/rad/mm, kx=0.1 N/mm^2, gravity load
#q = 1 kN/m
#Check again Figure 14 in the Plaut and Moen (2020) manuscript which was solved with Mathematica.
#  max ϕ

plaut_solution = 0.0023 #radians

#*********** Julia solution


num_beam_segments = 48
z = 0.0:7620.0/num_beam_segments:7620.0
num_nodes = length(z)

Ix = ones(Float64, num_nodes) * 3.230E6
Iy = ones(Float64, num_nodes) * 449530.0
Ixy = ones(Float64, num_nodes) * -865760.0
J = ones(Float64, num_nodes) * 397.09
Cw = ones(Float64, num_nodes) * 3.4104E9

E = ones(Float64, num_nodes) * 200.0
ν = ones(Float64, num_nodes) * 0.3
G = E./(2 .*(1 .+ ν))

ax = ones(Float64, num_nodes) * 27.826
ay = ones(Float64, num_nodes) * 101.6

kx = ones(Float64, num_nodes) * 0.1/1000  #kN/mm/mm
ay_kx = ones(Float64, num_nodes) * 101.6
kϕ = ones(Float64, num_nodes) * 1500/1000 #kN-mm/rad/mm

qx = ones(Float64, num_nodes) * 0.0
qy = ones(Float64, num_nodes) * 1.0/1000 #kN/mm

#u''=v''=ϕ''=0 (simply-supported), u'=v'=ϕ'=0  (fixed), u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free, e.g., a cantilever)
end_boundary_conditions = ["simply-supported", "simply-supported"]

#supports
#z, u, v, ϕ
supports = [(0.0, "fixed", "fixed", "fixed"),
            (7620.0, "fixed", "fixed", "fixed")]


model = ThinWalledBeam.solve(z, Ix, Iy, Ixy, J, Cw, E, G, kx, kϕ, ay_kx, qx, qy, ax, ay, end_boundary_conditions, supports)


using Plots 
plot(model.inputs.z, model.outputs.ϕ)

ϕmax = maximum(model.outputs.ϕ)

#calculate percent error as compared to Mathematica solution
error= abs.((ϕmax .- plaut_solution))./ plaut_solution

#accept 1% error from numerical solution
@test error <= 0.01
