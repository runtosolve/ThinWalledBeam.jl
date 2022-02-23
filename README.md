ThinWalledBeam.jl
==========

*Second order structural analysis of thin-walled beams*

Calculate deflection and twist of a thin-walled beam under vertical and lateral loads applied at an arbitrary location on the cross-section (top flange, bottom flange, ...).  The beam can be defined with lateral and rotational springs along its length, and it can have multiple support locations.

Install
-----------------------------

```julia
(v1.7) pkg> add ThinWalledBeam
```

(Type `]` to enter package mode.)

Example Usage
-------------

Solve Example 2 from [Plaut, R.H., Moen, C.D.(2021). "Lateral-Torsional Deformations of Single-Span and Two-Span Thin-Walled Beams with Continuous Bracing". Journal of Constructional Steel Research](https://www.sciencedirect.com/science/article/pii/S0143974X21000158?casa_token=52LwhLAs40wAAAAA:LjCOP5af8dSiEek3R40ggBnxfV0Svt89EdiwjRWF1_r-uF4eFvL8-r54udIpN7DeP9BZuQGA), see Figure 14.

Z-section, no slope, simple span
kϕ=1500 N*mm/rad/mm, kx=0.1 N/mm^2, gravity load
qy = 1 kN/m


```julia
using ThinWalledBeam


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

kx = ones(Float64, num_nodes) * 0.1/1000
ay_kx = ones(Float64, num_nodes) * 101.6
kϕ = ones(Float64, num_nodes) * 1500/1000

qx = ones(Float64, num_nodes) * 0.0
qy = ones(Float64, num_nodes) * 5.0/1000 #kN/mm

#u''=v''=ϕ''=0 (simply-supported), u'=v'=ϕ'=0  (fixed), u''=v''=ϕ''=u'''=v'''=ϕ'''=0 (free, e.g., a cantilever)
end_boundary_conditions = ["simply-supported", "simply-supported"]

#supports
#z, u, v, ϕ
supports = [(0.0, "fixed", "fixed", "fixed"),
            (7620.0, "fixed", "fixed", "fixed")]


model = ThinWalledBeam.solve(z, Ix, Iy, Ixy, J, Cw, E, ν, G, kx, kϕ, ay_kx, qx, qy, ax, ay, end_boundary_conditions, supports)



```

References
-----------------------------

[Plaut and Moen (2020)](https://cloud.aisc.org/SSRC/2020/43PlautandMoen2020SSRC.pdf)

[Plaut and Moen (2021)](https://www.sciencedirect.com/science/article/pii/S0143974X21000158?casa_token=52LwhLAs40wAAAAA:LjCOP5af8dSiEek3R40ggBnxfV0Svt89EdiwjRWF1_r-uF4eFvL8-r54udIpN7DeP9BZuQGA)