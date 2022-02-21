module ThinWalledBeam

using DiffEqOperators: CenteredDifference, LinearAlgebra, NLsolve

export define, solve


Base.@kwdef mutable struct Model   #using Parameters.jl macro here to help assign some constants, leave others for later.

   Ix::Union{Array{Float64}, Nothing} = nothing
   Iy::Union{Array{Float64}, Nothing} = nothing
   Ixy::Union{Array{Float64}, Nothing} = nothing
   J::Union{Array{Float64}, Nothing} = nothing
   Cw::Union{Array{Float64}, Nothing} = nothing
   E::Union{Array{Float64}, Nothing} = nothing
   ν::Union{Array{Float64}, Nothing} = nothing
   G::Union{Array{Float64}, Nothing} = nothing
   ax::Union{Array{Float64}, Nothing} = nothing
   ay::Union{Array{Float64}, Nothing} = nothing
   ay_kx::Union{Array{Array{Float64}}, Nothing} = nothing
   kx::Union{Array{Array{Float64}}, Nothing} = nothing
   kϕ::Union{Array{Array{Float64}}, Nothing} = nothing

   qx::Union{Array{Float64}, Nothing} = nothing
   qy::Union{Array{Float64}, Nothing} = nothing

   end_boundary_conditions::Union{Array{Int64}, Nothing} = nothing
   supports::Union{Vector{Tuple{Float64, String, String, String}}, Nothing} = nothing

   z::Union{Array{Float64}, Nothing} = nothing

   Kff::Union{Matrix{Float64}, Nothing} = nothing 
   Ff::Union{Array{Float64}, Nothing} = nothing

   free_dof_u::Union{Array{Int64}, Nothing} = nothing
   free_dof_v::Union{Array{Int64}, Nothing} = nothing
   free_dof_ϕ::Union{Array{Int64}, Nothing} = nothing

   u::Union{Array{Float64}, Nothing} = nothing
   v::Union{Array{Float64}, Nothing} = nothing
   ϕ::Union{Array{Float64}, Nothing} = nothing

end


function calculate_derivative_operators(dz)

   #Define the number of nodes.
   num_nodes = length(dz) + 1

   #Add extra dz on each end for padding nodes used in the CenteredDifference function.
   dz = [dz[1]; dz; dz[end]]

   #Calculate the 4th derivative operator.
   nth_derivative = 4
   derivative_order = 2
   Azzzz = CenteredDifference(nth_derivative, derivative_order, dz, num_nodes)

   #Convert the operator to a matrix.
   Azzzz = Array(Azzzz)   

   #Trim off the ghost nodes.
   Azzzz = Azzzz[:,2:end-1]  

   #Calculate the 2nd derivative operator.
   nth_derivative = 2
   derivative_order = 2
   Azz = CenteredDifference(nth_derivative, derivative_order, dz, num_nodes)
   Azz = Array(Azz)   
   Azz = Azz[:,2:end-1]  

   return Azzzz, Azz

end

function calculate_boundary_stencils(bc_flag, h, nth_derivative)

   #Calculate boundary conditions stencils without ghost nodes using
   #Jorge M. Souza, "Boundary Conditions in the Finite Difference Method"
   #https://cimec.org.ar/ojs/index.php/mc/article/download/2662/2607


   taylor_coeffs =  [1; h; h^2/2; h^3/6; h^4/24]

   RHS = zeros(Float64, 5)
   #row1*[A;B;C;D;E]*u2
   #row2*[A;B;C;D;E]*u2'
   #row3*[A;B;C;D;E]*u2''
   #row4*[A;B;C;D;E]*u2'''
   #row5*[A;B;C;D;E]*u2''''

   if bc_flag == 1 #simply supported end

      LHS = [1 1 1 1 0
         -1 0 1 2 0
         1 0 1 4 2
         -1 0 1 8 -6
         1 0 1 16 12]

      RHS[nth_derivative + 1] = (1/taylor_coeffs[nth_derivative + 1])
      boundary_stencil = LHS \ RHS
      boundary_stencil = ((boundary_stencil[1:4]),(zeros(4)))  #since u''=0

   elseif bc_flag == 2 #fixed end

      LHS = [1 1 1 1 0
         -1 0 1 2 1
         1 0 1 4 -2
         -1 0 1 8 3
         1 0 1 16 -4]

      RHS[nth_derivative+1] = (1/taylor_coeffs[nth_derivative + 1])
      boundary_stencil = LHS\RHS
      boundary_stencil = ((boundary_stencil[1:4]),(zeros(4)))  #since u'=0

   elseif bc_flag == 3 #free end
                #u'' u'''
      LHS = [1 1 1  0   0
           0 1 2  0   0
           0 1 4  0.5 0
           0 1 8  0   6
           0 1 16 0  0]

      RHS[nth_derivative+1] = (1/taylor_coeffs[nth_derivative+1])
      boundary_stencil1 = LHS\RHS  #at free end
      boundary_stencil1 = boundary_stencil1[1:3]   #u''=u'''=0

      # use simply supported BC to find stencil at one node in from free end
      LHS = [1 1 1 1 0
         -1 0 1 2 0
         1 0 1 4 2
         -1 0 1 8 -6
         1 0 1 16 12]

      RHS[nth_derivative+1] = (1/taylor_coeffs[nth_derivative+1])
      boundary_stencil2 = LHS\RHS #at one node over from free end
      boundary_stencil2 = boundary_stencil2[1:4]
      boundary_stencil = ((boundary_stencil1), (boundary_stencil2))  #two stencils are calculated

   end

   return boundary_stencil

end

function apply_end_boundary_conditions(A, end_boundary_conditions, nth_derivative, dz)

   #Consider the left end boundary condition.  Swap out interior stencil for boundary stencil.
   h = dz[1]
   bc_flag = end_boundary_conditions[1]
   boundary_stencil = calculate_boundary_stencils(bc_flag, h, nth_derivative)

   A[1,:] .= 0.0
   A[2,:] .= 0.0

   if (bc_flag == 1) | (bc_flag == 2)   #make this cleaner, combine
      A[2,1:length(boundary_stencil[1])] = boundary_stencil[1]
   else
      A[1,1:length(boundary_stencil[1])] = boundary_stencil[1]
      A[2,1:length(boundary_stencil[2])] = boundary_stencil[2]
   end

   #Consider the right end boundary condition.
   h = dz[end]
   bc_flag = end_boundary_conditions[2]
   boundary_stencil = calculate_boundary_stencils(bc_flag,h,nth_derivative)

   A[end,:] .= 0.0
   A[end-1,:] .= 0.0

   if (bc_flag == 1) | (bc_flag == 2)
      A[end-1,(end-(length(boundary_stencil[1])-1)):end] = reverse(boundary_stencil[1])
   else
      A[end,end-(length(boundary_stencil[1])-1):end] = reverse(boundary_stencil[1])
      A[end-1,end-(length(boundary_stencil[2])-1):end] = reverse(boundary_stencil[2])
   end

   return A

end

function define(z, section_properties, material_properties, kx, kϕ, ay_kx, qx, qy, ax, ay, end_boundary_conditions, supports)

   num_nodes = length(z)

   #Define properties at each node.
   Ix = ones(Float64, num_nodes) * section_properties.Ixx
   Iy = ones(Float64, num_nodes) * section_properties.Iyy
   Ixy = ones(Float64, num_nodes) * section_properties.Ixy
   J = ones(Float64, num_nodes) * section_properties.J
   Cw = ones(Float64, num_nodes) * section_properties.Cw

   E = ones(Float64, num_nodes) * material_properties.E
   ν = ones(Float64, num_nodes) * material_properties.ν
   G = E./(2 .*(1 .+ ν))

   #Store model inputs in data structure.
   model = Model(Ix, Iy, Ixy, J, Cw, E, ν, G, ax, ay, ay_kx, kx, kϕ, qx, qy, end_boundary_conditions, supports, z, nothing, nothing, nothing, nothing, nothing, nothing, nothing,  nothing)
  
   return model

end

function define_fixed_dof(z, support_location, support_fixity)

   if support_fixity == "fixed"

      fixed_dof = findall(x->x≈support_location, z)

   else

      fixed_dof = []

   end

   return fixed_dof

end


function governing_equations(model)

    #The equilbrium equations come from Plaut and Moen(2019) https://cloud.aisc.org/SSRC/2020/43PlautandMoen2020SSRC.pdf.

   #Define the number of nodes.
   num_nodes = length(model.z)

   #Calculate the derivative operators.
   dz = diff(model.z) 
   Azzzz,Azz = calculate_derivative_operators(dz) 

   #Apply left and right end boundary condition stencils to derivative operators.
   nth_derivative = 4
   Azzzz = apply_end_boundary_conditions(Azzzz, model.end_boundary_conditions, nth_derivative, dz)

   nth_derivative = 2
   Azz = apply_end_boundary_conditions(Azz, model.end_boundary_conditions, nth_derivative, dz)

   #Build identity matrix for ODE operations.
   AI = Matrix(1.0I, num_nodes, num_nodes)

   #Build operator matrix that doesn't update with load.

   #LHS first.
   A11 = zeros(Float64, num_nodes,num_nodes)
   A12 = zeros(Float64, num_nodes,num_nodes)
   A13 = zeros(Float64, num_nodes,num_nodes)
   A21 = zeros(Float64, num_nodes,num_nodes)
   A22 = zeros(Float64, num_nodes,num_nodes)
   A23 = zeros(Float64, num_nodes,num_nodes)
   A31 = zeros(Float64, num_nodes,num_nodes)
   A32 = zeros(Float64, num_nodes,num_nodes)
   A33 = zeros(Float64, num_nodes,num_nodes)

   #Sum spring terms if there are more than one on a cross-section.
   sum_kx = sum(model.kx)
   sum_kϕ = sum(model.kϕ)

   num_kx_springs = size(model.kx)[1]
   kx_ay_kx = [model.kx[i][:] .* model.ay_kx[i][:] for i=1:num_kx_springs] #multiply kx*ay_kx terms for each spring first
   sum_kx_ay_kx = sum(kx_ay_kx) #then sum them
   
   kx_ay_kx_ay_kx = [model.kx[i][:] .* model.ay_kx[i][:] .* model.ay_kx[i][:] for i=1:num_kx_springs]
   sum_kx_ay_kx_ay_kx = sum(kx_ay_kx_ay_kx)

   #Calculate operator quantities on LHS  AU=B.
   for i = 1:num_nodes
      A11[i,:] = model.E .* model.Iy .* Azzzz[i,:] .+ sum_kx .* AI[i,:]
      A12[i,:] = model.E .* model.Ixy .* Azzzz[i,:]
      A13[i,:] = sum_kx_ay_kx .*AI[i,:]
      A21[i,:] = model.E .* model.Ixy .* Azzzz[i,:]
      A22[i,:] = model.E .* model.Ix .* Azzzz[i,:]
      A31[i,:] = sum_kx_ay_kx .* AI[i,:]
      A33[i,:] = model.E .* model.Cw .* Azzzz[i,:] .- model.G .* model.J .* Azz[i,:] .+ sum_kx_ay_kx_ay_kx .* AI[i,:] .+ sum_kϕ .* AI[i,:] .+ model.qx .* model.ax .* AI[i,:] .- model.qy .* model.ay .* AI[i,:]
   end

   #Calculate RHS of AU=B.
   B1 = model.qx
   B2 = model.qy
   B3 = model.qx .* model.ay .+ model.qy .* model.ax

   #Reduce problem to free dof.

   fixed_dof_u = Array{Int64}(undef, 1)
   fixed_dof_v = Array{Int64}(undef, 1)
   fixed_dof_ϕ = Array{Int64}(undef, 1)

   for i = 1:length(model.supports)

      if i == 1

         fixed_dof_u = define_fixed_dof(model.z, model.supports[i][1], model.supports[i][2])
         fixed_dof_v = define_fixed_dof(model.z, model.supports[i][1], model.supports[i][3])
         fixed_dof_ϕ = define_fixed_dof(model.z, model.supports[i][1], model.supports[i][4])

      else

         new_fixed_dof_u = define_fixed_dof(model.z, model.supports[i][1], model.supports[i][2])
         new_fixed_dof_v = define_fixed_dof(model.z, model.supports[i][1], model.supports[i][3])
         new_fixed_dof_ϕ = define_fixed_dof(model.z, model.supports[i][1], model.supports[i][4])
         
         fixed_dof_u = [fixed_dof_u; new_fixed_dof_u]
         fixed_dof_v = [fixed_dof_v; new_fixed_dof_v]
         fixed_dof_ϕ = [fixed_dof_ϕ; new_fixed_dof_ϕ]

      end

   end

   #Find free dof. 
   model.free_dof_u = setdiff(1:num_nodes,fixed_dof_u)
   model.free_dof_v = setdiff(1:num_nodes,fixed_dof_v)
   model.free_dof_ϕ = setdiff(1:num_nodes,fixed_dof_ϕ)

   #Assemble global K.
   K = [A11 A12 A13;
   A21 A22 A23;
   A31 A32 A33]

   #Partition stiffness matrix.
   free_dof = [model.free_dof_u; model.free_dof_v .+ num_nodes; model.free_dof_ϕ .+ 2*num_nodes]

   model.Kff = K[free_dof, free_dof]

   #Define external force vector.
   F = [B1; B2; B3]

   model.Ff = F[free_dof]

   return model

end

function residual!(R, U, K, F)

   for i=1:length(F)
    R[i] = transpose(K[i,:])*U-F[i]
   end

   return R

end


function solve(model)

   #Set up solution matrices from governing equations.
   model = governing_equations(model)

   #Define the number of nodes along the beam.
   num_nodes = length(model.z)

   #Define the deformation vectors.
   u = zeros(Float64, num_nodes)
   v = zeros(Float64, num_nodes)
   ϕ = zeros(Float64, num_nodes)

   #Define the deformation initial guess for the nonlinear solver.
   deformation_guess = model.Kff \ model.Ff

   #Solve for the beam deformations.
   solution = nlsolve((R,U) ->residual!(R, U, model.Kff, model.Ff), deformation_guess)

   #Pull the displacements and twist from the solution results.
   u[model.free_dof_u] = solution.zero[1:length(model.free_dof_u)]
   v[model.free_dof_v] = solution.zero[length(model.free_dof_u)+1:(length(model.free_dof_u) + length(model.free_dof_v))]
   ϕ[model.free_dof_ϕ] = solution.zero[(length(model.free_dof_u) + length(model.free_dof_v))+1:(length(model.free_dof_u) + length(model.free_dof_v)+length(model.free_dof_ϕ))]

   #Add deformations to data structure.
   model.u = u
   model.v = v
   model.ϕ = ϕ

   return model

end

end # module
