### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 8733a11b-ab6c-43e3-b4c1-f126a3f5640c
begin
	using PlutoUI, PlutoTeachingTools, LinearAlgebra
	TableOfContents()
end

# ╔═╡ cb298b2f-c4fe-4045-b2fd-64c561ab9da5
begin
	using Plots
	
	Plots.plot([0,0,1,1,0], [0,1,1,0,0], color=:black, linewidth=1, xlim=[-0.2, 1.5], ylim=[-0.1, 1.1], legend=false,size=(700,500), axis = false, grid=false)
	scatter!([0.5, 0.5, 0, 1], [0, 1, 0.5, 0.5], color=:black, markersize=5
	)
	plot!([0.5, 0.65], [0,0], arrow=true, linewidth=3, color=:black)
	plot!([0.5, 0.65], [1,1], arrow=true, linewidth=3, color=:black)
	plot!([0,0], [0.5, 0.65], arrow=true, linewidth=3, color=:black)
	plot!([1,1], [0.5, 0.65], arrow=true, linewidth=3, color=:black)
	scatter!([0.5], [0.5], markershape=:circle, markercolor=:white, markersize=8)
	scatter!([0.5], [0.5], markershape=:circle, markercolor=:black, markersize=3)
	annotate!(0.5, 0.58, "H_z(i+1/2, j+1/2)", annotationfontsize=10)
	annotate!(0.5, -0.08, "E_x(i+1/2, j)", annotationfontsize=10)
	annotate!(0.5, 1.08, "E_x(i+1/2, j+1)", annotationfontsize=10)
	annotate!(-0.15, 0.5, "E_y(i, j+1/2)", annotationfontsize=10)
	annotate!(1.18, 0.5, "E_y(i+1, j+1/2)", annotationfontsize=10)
end

# ╔═╡ befb8334-b8cc-46d5-bede-246612c9feae
using PlotlyJS

# ╔═╡ c3ebeb58-24c4-11ef-2fd7-a55d234fb215
md"""
# Finite-Difference Time-Domain Simulations of Electromagnetic Waves in 2D

The behavior of electromagnetic waves at interfaces between different media is a fundamental problem in electromagnetics. The finite-difference time-domain method is a powerful and straightforward approach to simulate such problems. 
In the following, we will implement a finite-difference time-domain (FDTD) method and simulate the behavior of electromagnetic waves at interfaces between dielectric media and scattering on perfectly electrically conducting (PEC) surfaces.
"""

# ╔═╡ f2290f80-3ff5-42cf-8d3b-a152d8bbf578
md"""
## Propaedeutics and Theory
Electromagnetic waves in a linear and non-dispersive medium are typically described by the wave equation
```math
\begin{equation}
\nabla \times \nabla \times \mathbfit E = \Delta E = - \mu \epsilon \frac{\partial^2 \mathbfit E}{\partial t^2}\,.
\end{equation}
```
As an alternative to the wave equation, which is a second-order partial differential equation, we can describe electromagnetic waves by a system of coupled first-order partial differential equations, the Maxwell equations 
```math
\begin{aligned}
& \varepsilon \frac{\partial E_x}{\partial t} = \frac{\partial H_z}{\partial y} - \frac{\partial H_y}{\partial z}\,,\\ 
& \varepsilon \frac{\partial E_y}{\partial t} = \frac{\partial H_x}{\partial z} - \frac{\partial H_z}{\partial x}\,,\\ 
& \varepsilon \frac{\partial E_z}{\partial t} = \frac{\partial H_y}{\partial x} - \frac{\partial H_x}{\partial y}\,,\\ 
\\
&\mu \frac{\partial H_x}{\partial t} = \frac{\partial E_y}{\partial z} - \frac{\partial E_z}{\partial y}\,,\\
&\mu \frac{\partial H_y}{\partial t} = \frac{\partial E_z}{\partial x} - \frac{\partial E_x}{\partial z}\,,\\
&\mu \frac{\partial H_z}{\partial t} = \frac{\partial E_x}{\partial y} - \frac{\partial E_y}{\partial x}\,.\\
\end{aligned}
```
"""

# ╔═╡ 341beaad-0527-4f37-8666-71cfc80a4bae
md"""
### TE-Waves in 2D
Implementing this problem in three dimensions can be quite tedious. To reduce the implementation effort, we will limit ourselves to simulating transverse electric waves (TE waves) in two dimensions. 
Therefore, we only need to consider the field components $H_z$, $E_x$, and $E_y$, leaving only the equations
```math
\begin{aligned}
& \varepsilon \frac{\partial E_x}{\partial t} = \frac{\partial H_z}{\partial y}\,,\\ 
& \varepsilon \frac{\partial E_y}{\partial t} = - \frac{\partial H_z}{\partial x}\,,\\ 
&\mu \frac{\partial H_z}{\partial t} = \frac{\partial E_x}{\partial y} - \frac{\partial E_y}{\partial x}\,.\\
\end{aligned}
```

### Yee Cell
To achieve good results in the simulation, we use the so-called Yee cell, which defines the positions of the discretized electric and magnetic fields. 
The basic principle is that we place the electric field on an integer grid and the magnetic field on a half grid.
"""

# ╔═╡ 61a4d03e-cbf1-4ff2-a915-b68dddfacba5
md"""
### Discretization
!!! note \"Problem 1: Discretized Equations\"
	Discretize the following equations using the the Yee Cell and a staggered grid finite difference approach.
	```math
		\begin{aligned}
			& \varepsilon \frac{\partial E_x}{\partial t} = \frac{\partial H_z}{\partial y}\,,\\ 
			& \varepsilon \frac{\partial E_y}{\partial t} = - \frac{\partial H_z}{\partial x}\,,\\ 
			&\mu \frac{\partial H_z}{\partial t} = \frac{\partial E_x}{\partial y} - \frac{\partial E_y}{\partial x}\,.\\
		\end{aligned}
	```
!!! tip
	You can use the nomenclature $E_x(i, j, n)$, where $i$ is the dicretization in $x$-direction, $j$ in $y$-direction, and $n$ the time step.
!!! tip
	An example of how to code equations in Markdown is shown below.
"""

# ╔═╡ c9433bc7-274e-439b-bfa7-4da9d9b32480
md"""
#### Example:
```math
\begin{multline}
\frac{A(i, j, n)+B(i, j, n)+C(i, j, n)}{\Delta abc} = \\
\frac{A(i, j, n)+B(i, j, n)+C(i, j, n)}{\Delta abc} + \frac{A(i, j, n)+B(i, j, n)+C(i, j, n)}{\Delta abc}
\end{multline}
```
"""

# ╔═╡ fd5de002-9e52-4483-81c8-edffa078a83a
#Write down your solution here


md"""
#### Discretized Equations:
```math
\begin{aligned}

& \varepsilon \frac{E_x(i+\frac{1}{2},j,n+1)-E_x(i+\frac{1}{2},j,n)}{\Delta t} \\
& = \frac{H_z(i+\frac{1}{2},j+\frac{1}{2},n+\frac{1}{2})-H_z(i+\frac{1}{2},j-\frac{1}{2},n+\frac{1}{2})}{\Delta y}

\\ \\ \\

& \varepsilon \frac{E_y(i,j+\frac{1}{2},n+1)-E_y(i,j+\frac{1}{2},n)}{\Delta t} \\ & = -\frac{H_z(i+\frac{1}{2},j+\frac{1}{2},n+\frac{1}{2})-H_z(i-\frac{1}{2},j+\frac{1}{2},n+\frac{1}{2})}{\Delta x}

\\ \\ \\

& \mu \frac{H_z(i+\frac{1}{2},j+\frac{1}{2},n+\frac{1}{2})-H_z(i+\frac{1}{2},j+\frac{1}{2},n-\frac{1}{2})}{\Delta t} \\ 
& = \frac{E_x(i+\frac{1}{2},j+1,n)-E_x(i+\frac{1}{2},j,n)}{\Delta y} - \frac{E_y(i+1,j+\frac{1}{2},n)-E_y(i,j+\frac{1}{2},n)}{\Delta x}


\end{aligned}
```
"""

# ╔═╡ 77501480-6a29-4e3d-9d71-2c2d91a0f6bd
md"""
## Implementation
### Update equations
For the implementation, we cannot directly use the halfgrid as in the theory. The halfgrid is enforced by using different arrays to store the electric and magnetic field components and updating them in the correct order. 
To get an additional spatial step for the electric field, we need to add an entry to the array of $E_x$ in $y$-direction and one entry to the array of $E_y$ in $x$-direction.
To implement the equations of Problem 1, we solve each equation for the $n+1$ time step and "neglect" the half steps. 
!!! note \"Problem 2: Discretized Equations\"
	Write down the equations from Problem 1, solved for the $n+1$ time step.
"""

# ╔═╡ 0f4b74d2-23f5-4061-beaa-e6abc44d7afb
#Write down your solution here
md"""
#### Discretized Equations:
```math
\begin{aligned}

& E_x(i,j,n+1) = E_x(i,j,n) + \frac{\Delta t}{\varepsilon \Delta y} [H_z(i,j,n)-H_z(i,j-1,n)]

\\ \\ \\

& E_y(i,j,n+1) = E_y(i,j,n) - \frac{\Delta t}{\varepsilon \Delta x} [H_z(i,j,n)-H_z(i-1,j,n)]

\\ \\ \\

& H_z(i,j,n)=H_z(i,j,n-1) \\ &
+ \frac{\Delta t}{\mu \Delta y}[E_x(i,j+1,n)-E_x(i,j,n)] -\frac{\Delta t}{\mu \Delta x}[E_y(i+1,j,n)-E_y(i,j,n)]

\end{aligned}
```
"""

# ╔═╡ 46dc4504-69c7-49c7-a217-804b3a651c84
md"""
### Time stepping
!!! note \"Problem 3: Time Stepping\"
	The choice of $\Delta t$ is in the FDTD rather important to guarantee a numerically stable simulation. How do you choose $\Delta t$ for this simulation?
"""

# ╔═╡ 57673298-47f9-402b-a63b-0263adbb238a
#Write down your solution here

md"""

We choose $\Delta t$ according to the Courant-Friedrichs-Lewy (CFL) condition. 

CFL condition in 2D is:
```math


\Delta t \leq \frac{1}{c} \left( \frac{1}{\frac{1}{\Delta x^2} + \frac{1}{\Delta y^2}} \right)^{1/2}
```

where:
- The speed of light in the medium is $c$, where $c = \frac{1}{\sqrt{\mu \varepsilon}}$.
- The spatial step sizes in the $x$ and $y$ directions are $\Delta x$ and $\Delta z$.

For $\Delta x = \Delta y = h$: 
```math
\Delta t \leq \frac{h}{c \sqrt{2}} 
```

"""

# ╔═╡ 4dd8f7a1-f5cd-49ae-89ab-ce89da43ff65
md"""
### Constants
The constants we need for the simulation are given by the following struct and are passed to our main function.
"""

# ╔═╡ f5c03c4a-067d-432f-9205-c0693394faa4
struct PhysicalConstants{F}
	c₀::F
	ϵ₀::F
	μ₀::F
end

# ╔═╡ 239595e8-89e2-4f1e-a503-ec30dbdf1a85
function constants(F::Type)
    return PhysicalConstants(F(2.997925e8), F(8.854187817e-12), F(4 * π * 1e-7))
end

# ╔═╡ 3ae46ec8-1755-4dee-8c51-8406d0cb5e3c
md"""
### Simulation Parameter
All parameters needed for the simulation are defined in the struct `SimRegion`. The parameter `M` is of the type `Material{F}`, which is itself a struct containing the material properties at each spatial point. 
The parameters in the struct `Material{F}` are defined by
```math
\begin{align}
	m_{E} &= \frac{\Delta t}{\varepsilon(x, y)\Delta s}\,,\\
	m_{H_z} &= \frac{\Delta t}{\mu(x, y)\Delta s}\,.
\end{align}
```
"""

# ╔═╡ c3522930-99d2-4e61-a591-2fb2ae27b1e8
struct Material{F}
	mE_x::Matrix{F}
	mE_y::Matrix{F}
	mH_z::Matrix{F}
end

# ╔═╡ ecd3a994-3c72-432d-bac2-b4fa884631cb
struct SimRegion{I, F}
	N_x::I 	#number of spatial steps in x
	N_y::I 	#number of spatial steps in y
	T::I 	#number of time steps
	M::Material{F} 	#material parameters on each spatial point
	Δt::F 	#time step
	Δs::F 	#spatial step
end

# ╔═╡ 4cf99d73-fc9b-468e-8b47-862a63a3f5a7
md"""
### Excitations
Later we will simulate a point source and a plane wave, the signal shape for both will be a Gaussian pulse.
"""

# ╔═╡ 583cef9d-16a4-424d-a42e-e115540b0ae8
function gaussian_pulse(
	t::F, 	#time
	μ::F, 	#offset
	σ::F 	#width of the gaussian pulse
) where F
   	return exp(-(t-μ)^2/σ^2)
end

# ╔═╡ d35215d1-f8b3-42e7-b222-610763e0476a
md""" 
For dispaching, we define two structs, one for the point source and one for the plane wave:
"""

# ╔═╡ b726aaf2-e02b-4959-9058-b86ff6c98dd7
abstract type Excitation{F} end

# ╔═╡ 218d9e15-d242-4625-ae9d-37f4408be3e9
struct PointSource{F} <: Excitation{F}
	μ::F
	σ::F
end

# ╔═╡ 060ec34f-9653-46ef-9a6e-12f0cac61153
struct PlaneWave{F} <: Excitation{F}
	μ::F
	σ::F
end

# ╔═╡ 7cc1227e-199d-4497-ab14-cf61c9d24834
md"""
The main function will receive one of the two structs and will call an overloaded function `excite!()`, which will be able to dispatch, depending on the function it receives, waether to excite with a plane wave or a point source.
"""

# ╔═╡ 7644de2a-dbbd-475a-8926-29ee439579f2
function excite!(
	H_z::Array{F}, 				#array containing the magnetic field
	n::I,  						#time step
	sr::SimRegion{I, F}, 		#struct defining the simulation region 
	inc_wave::PointSource{F} 	#excitation either a PointSource{F} or a PlaneWave
) where {I, F}

	H_z[50, Int(size(H_z)[1]/2), n] = gaussian_pulse(
		n*sr.Δt,
		inc_wave.μ,
		inc_wave.σ
	)
end

# ╔═╡ a6dc0326-67ea-4a2c-8bcd-bff2b97b3847
md"""
### Simulation
If implemented correctly, we only need one main function to simulate different scenarios. The main function has two parameters: `sr::SimRegion{I, F}`, which contains all parameters describing the simulation region including the material, and `excitation::Excitation{F}`, which is either of the type `PointSource{F}` or of the type `PlaneWave{F}`, defining the excitation of the simulation. 

The last part is how to handle the boundaries of our simulation domain. For the sake of simplicity, we will assume a PEC at the boundary, which means that we have to set the components of the electric field parallel to the boundary to zero after each time step update.

!!! warning \" Problem 4: FDTD Main Function\"
	Complete the implementation of the `FDTD_2D()` function following the given steps.
!!! tip
	I would recommend to implement the complete function in the following steps:

    1. Start by allocating the arrays to store the electric and magnetic fields. Note that we need to store every time step for the magnetic field. Therefore, we will use only a two-dimensional array for the electric field components and a three-dimensional array for the magnetic field.

    2. Implement the update equations in the following order: a) Update the magnetic field. b) Call the excitation function. c) Update the electric field. d) Enforce the PEC boundary condition.
"""

# ╔═╡ 8b4e721b-3a99-4070-a7b6-7999956afba2
function excite!(
	H_z::Array{F}, 			#array containing the magnetic field
	n::I, 					#time step
	sr::SimRegion{I, F}, 	#struct defining the simulation region
	inc_wave::PlaneWave{F} 	#excitation either a PointSource{F} or a PlaneWave
) where {I, F}
	H_z[1, :, n] .= 1/2 * gaussian_pulse(
		n*sr.Δt,
		inc_wave.μ,
		inc_wave.σ,
	)
end

# ╔═╡ 8c0155fe-493b-4221-8658-64e68491d8e7
begin
	F= Float64
	C = constants(F) # contains struct PhysicalConstants of type float64
end

# ╔═╡ 5b1c7301-2e94-460f-b798-ecb8611e6017
md"""
### Visualization
To visualize our simulation we can use the following function `visualize()`.
The only parameter the function takes is a three-dimensional array, where the first and second dimensions are the spatial dimensions and the third is the time step.
"""

# ╔═╡ e33ca02f-6a79-480c-849d-32a32e3b9739
function visualize(H_z)
	trace = PlotlyJS.surface(z = H_z[:,:,1],cmin=-1,cmax=1)
	n_frames = Int(length(H_z[1,1,:])/5)
	frames = Vector{PlotlyFrame}(undef, n_frames)
	
	for k in 1:n_frames
	    frames[k] = PlotlyJS.frame(
			data=[attr(z = H_z[:,:,k*5])],
	        name="$k", #frame name; it is passed to slider
	        traces=[0] # this means that the above data update the first trace (here the unique one)
	    )
	end
	
	updatemenus = [attr(
	    type="buttons",
	    active=0,
	    y=1, #(x,y) button position
	    x=1.1,
	    buttons=[attr(
	        label="Play",
	        method="animate",
	        args=[
	            nothing,
	            attr(
	                frame=attr(duration=0.1, redraw=true),
	                transition=attr(duration=0.1),
	                fromcurrent=true,
	                mode="immediate"
	            )]
	    )]
	)]
	
	sliders = [attr(
	    active=1,
	    minorticklen=0.1,
	    steps=[attr(
	        label="$k",
	        method="animate",
	        args=[["$k"], attr(
	            mode="immediate",
	            transition=attr(duration=0.1),
	            frame=attr(duration=0.1,
	            redraw=true)
	        )]) for k in 1:n_frames]
	)]


	layout = Layout(
	    title_text="",
	    title_x=0.5,
	    width=500,
	    height=500,
	    scene = attr(zaxis=attr(nticks=4, range=[-1,1])),
	    updatemenus=updatemenus,
	    sliders=sliders
	)
	
	pl = Plot(trace, layout, frames)
end

# ╔═╡ 90e035ef-eded-4d79-b652-5a4643df57a6
function FDTD_2D(
    sr::SimRegion{I, F},
    excitation::Excitation{F}
 ) where {I, F}
    # Allocate zero arrays for storing the fields
    E_x = zeros(F, sr.N_x+1, sr.N_y) # Stores the x-component of the electric field
    E_y = zeros(F, sr.N_x, sr.N_y+1) # Stores the y-component of the electric field
    H_z = zeros(F, sr.N_x, sr.N_y, sr.T) # Stores the z-component of the magnetic field
    
    # Main FDTD loop
    for n in 2:sr.T
    # Update magnetic field at each grid point
    for ii in 1:sr.N_x-1
        for jj in 1:sr.N_y-1
			H_z[ii, jj, n] = H_z[ii, jj, n-1] .+ sr.M.mH_z[ii, jj] .* (E_x[ii, jj+1] .- E_x[ii, jj] .- E_y[ii+1, jj] .+ E_y[ii, jj])
            end
        end
        
        # Apply excitation
		excite!(H_z, n, sr, excitation)
        
        # Update electric field E_x components at each grid point
        for ii in 1:sr.N_x-1
            for jj in 2:sr.N_y-1
				E_x[ii, jj] = E_x[ii, jj] .+ sr.M.mE_x[ii, jj] .* (H_z[ii, jj, n] .- H_z[ii, jj-1, n])
			end
		end

		# Update electric field E_y components at each grid point
		for ii in 2:sr.N_x-1
            for jj in 1:sr.N_y-1				
				E_y[ii, jj] = E_y[ii, jj] .- sr.M.mE_x[ii, jj] .* (H_z[ii, jj, n] .- H_z[ii-1, jj, n])
            end
        end
		
		# Implement the Perfect Electric Conductor (PEC) boundary conditions on the electric field components at the edges of the simulation region
		E_x[:, 1] .= 0  
		E_x[:, end] .= 0  
		E_y[1, :] .= 0  
		E_y[end, :] .= 0  
    end
	
 return H_z
 end

# ╔═╡ b613c16e-66d2-49f7-8527-112335f97d58
md"""
## Results
If the `FDTD_2D`function is implemented correctly, the following code will produce a simulation of a point source.
"""

# ╔═╡ 2dd13c57-59a3-4c72-9aad-eb7c603ae5e1
begin 
	N_x= 200 # number of grid points in the x-direction 
	N_y = 200 # number of grid points in the y-direction 
	T = 300 # total number of time steps
	Δs = 1/100 # spatial step
	Δt = Δs/(C.c₀*sqrt(2)) # time step that satisfy CFL condition
	M = Material(
		Δt / (C.ϵ₀ * Δs) .* ones(F, N_x, N_y+1), 
		Δt / (C.ϵ₀ * Δs) .* ones(F, N_x+1, N_y),
		Δt / (C.μ₀ * Δs) .* ones(F, N_x, N_y)
	)
	sr1 = SimRegion(N_x, N_y, T, M, Δt, Δs) # create a simulation region
	excitation1 = PointSource(40*Δt, 10*Δt) # excitation with point source
	H_z1 = FDTD_2D(sr1, excitation1) # Calls the FDTD_2D function to perform the FDTD simulation with the excitation of type point source
	visualize(H_z1)
end

# ╔═╡ ed7ecc33-9774-4c6d-86a6-7eb78befe270
md"""
!!! warning \"Problem 5: Plan Wave\"
	Simulate and visualize a plane wave with the same parameters.
"""

# ╔═╡ bfbc714e-ca3a-4a4e-a704-1342f477a5c7
begin 	
	sr2 = SimRegion(N_x, N_y, T, M, Δt, Δs) # create a simulation region
	excitation2 = PlaneWave(40*Δt, 10*Δt) # excitation of type Plane wave
	H_z2 = FDTD_2D(sr2, excitation2) # Calls the FDTD_2D function to perform the FDTD simulation with the excitation of type plane wave
	visualize(H_z2) # visualize the plane wave
end

# ╔═╡ cbf17821-1b31-4558-8855-c5451e8e2cf5
md"""
!!! warning \"Problem 5: PEC Scatterer\"
	Write a function 
	```
	pec_cylinder!(
		pos::Vector{I}, r::I, m::Material{F}
	) where {I, F}
	```
	which modifies the material properties according to a PEC cylinder placed at the position `pos::Vector{I}` with the radius `r::I`. Both the position and the radius are given in spatial steps and must be integers.
"""

# ╔═╡ 9da8113b-43fb-4439-9e5f-6e4d3bd69028
# The PEC imposes boundary conditions where the electric field inside and at the surface of the cylinder must be zero

function pec_cylinder!(pos::Vector{I}, r::I, M::Material{F}) where {I, F}
    # Define position of the cylinder
    x_c, y_c = pos 

    # Modify material properties in the simulation region
    for ii in 1:size(M.mE_x, 1)
        for jj in 1:size(M.mE_x, 2)
            # Calculate distance from the cylinder center
            distance = sqrt((ii - x_c)^2 + (jj - y_c)^2)
            if distance <= r
				# Modify the following components to be 0 inside and at the cylinder
                M.mE_x[ii, jj] = 0.0  # Set E-field components to zero
                M.mE_y[ii, jj] = 0.0  # Set E_y to zero
                M.mH_z[ii, jj] = 0.0  # Set H-field components to zero
            end
        end
    end
end

# ╔═╡ 5940de23-eca9-411b-b43c-30bc8f8f23d3
md"""
!!! warning \"Problem 6: PEC Scatterer\"
	Simulate and visualize the scattering of a plane wave on a cylinder using the same parameters as before.
"""

# ╔═╡ a3a612b9-a454-4692-b85b-74d6017f7818
begin
	
	# Define the simulation region
	sr3 = SimRegion(N_x, N_y, T, M, Δt, Δs)
	
	# Define the PEC cylinder parameters
	pos = [Int(N_x/2), Int(N_y/2)]  # Position origin of the cylinder at the center of the grid
	r = 10            # Radius of the cylinder 
	
	# Apply the PEC scatterer (this modifies the material properties)
	pec_cylinder!(pos, r, M)
	
	# Define the plane wave excitation 
	excitation3 = PlaneWave(40 * Δt, 10 * Δt)
	
	# Run the FDTD simulation
	H_z3 = FDTD_2D(sr3, excitation3)
	
	# Visualize the results (scattered field)
	visualize(H_z3)
end

# ╔═╡ 1153a472-58f2-4417-a679-b9d77fb15937
md"""
!!! warning \" Problem 7: Interface\"
	Implement a function 
	```
	interface!(
		pos::I,
		μ₁::F,
		μ₂::F,
		ϵ₁::F,
		ϵ₂::F,
		M::Material{F}
	) where {I, F}
	```
	which modifies the material properties according to an interface between two different materials. Note that $\mu_1, \mu_2, \varepsilon_1,$ and $\varepsilon_2$ are the ralative permitivities and permeabilities of the two materials.
"""

# ╔═╡ 86808aa3-68ac-49a4-b7d4-df4fa95873b3
struct PhysicalConstants2{F} # define a new struct for the relative permittivity and permeability of water
	ϵ₁::F
	μ₁::F
end

# ╔═╡ 8334b234-e6d6-486d-a7bb-8d4b4be1a5ff
begin
	function constants2(F::Type)
	    return PhysicalConstants2(F(1.5671912*10^-11), 
			# permittivity of water ϵ₁ = ϵ₀ϵr = (8.854187817e-12)*1.77 = 1.5671912e-11
			F(1.2566*10^-6)) # permeability of water μ₁ = μ₀μr = (4 * π * 1e-7)*(1−9.1·10^{−6}) = 1.2566×10^-6
	end
	C2 = constants2(F)
end

# ╔═╡ 84a992bc-0166-4485-bdef-d6f6aa1b21fe
 function interface!(
    pos2::I, # position of the interface
    μ₀::F, # permeability of air
    μ₁::F, # permeability of water
    ϵ₀::F, # permittivity of air
    ϵ₁::F, # # permittivity of water
    M::Material{F}
 ) where {I, F}
    
    # Iterate through the grid to modify material properties at the interface
    for ii in 1:N_x
        for jj in 1:N_y
		if ii >= pos2  # Change properties to water after it hits the PEC interface
                M.mH_z[ii, jj] = Δt / (C2.μ₁ * Δs)
                M.mE_x[ii, jj] = Δt / (C2.ϵ₁ * Δs)
                M.mE_y[ii, jj] = Δt / (C2.ϵ₁ * Δs)
            else  # Keep properties as air (vacuum) before the PEC interface
                M.mH_z[ii, jj] = Δt / (C.μ₀ * Δs)
                M.mE_x[ii, jj] = Δt / (C.ϵ₀ * Δs)
                M.mE_y[ii, jj] = Δt / (C.ϵ₀ * Δs)
            end
        end
    end
 end

# ╔═╡ 8858623a-c619-4f4c-8e80-41b471ee5fa2
md"""
!!! warning \" Problem 8: Interface \"
	Simulate the interface between air (vacuum) and water using the same parameters as before ($\epsilon_{r}=1.77$, $\mu_{r}=1−9.1·10^{−6}$) 
"""

# ╔═╡ fbea9fa5-2ce0-4db1-bc28-b5e7a85bdd57
begin

	# Define the simulation region
	sr4 = SimRegion(N_x, N_y, T, M, Δt, Δs)
	
	# Define the position of the interface 
	pos2 = Int(N_x/2)  
	
	# Apply the interface between air and water
	interface!(pos2, C.μ₀, C2.μ₁, C.ϵ₀, C2.ϵ₁, M)
	
	# Define the plane wave excitation 
	excitation4 = PlaneWave(40 * Δt, 10 * Δt)
	
	# Run the FDTD simulation
	H_z4 = FDTD_2D(sr4, excitation4)
	
	# Visualize the results 
	visualize(H_z4)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlotlyJS = "f0f68f2c-4968-5e81-91da-67840de0976a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlotlyJS = "~0.18.13"
Plots = "~1.40.4"
PlutoTeachingTools = "~0.2.15"
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0"
manifest_format = "2.0"
project_hash = "2c6d51a6cca8b2aa4c082a4e2c8512d939bb3a7a"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AssetRegistry]]
deps = ["Distributed", "JSON", "Pidfile", "SHA", "Test"]
git-tree-sha1 = "b25e88db7944f98789130d7b503276bc34bc098e"
uuid = "bf4720bc-e11a-5d0c-854e-bdca1663c893"
version = "0.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BinDeps]]
deps = ["Libdl", "Pkg", "SHA", "URIParser", "Unicode"]
git-tree-sha1 = "1289b57e8cf019aede076edab0587eb9644175bd"
uuid = "9e28174c-4ba2-5203-b857-d8d62c4213ee"
version = "1.0.2"

[[deps.Blink]]
deps = ["Base64", "BinDeps", "Distributed", "JSExpr", "JSON", "Lazy", "Logging", "MacroTools", "Mustache", "Mux", "Reexport", "Sockets", "WebIO", "WebSockets"]
git-tree-sha1 = "08d0b679fd7caa49e2bca9214b131289e19808c0"
uuid = "ad839575-38b3-5650-b840-f874b8c74a25"
version = "0.12.5"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "71acdbf594aab5bbb2cec89b208c41b4c411e49f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.24.0"

[[deps.ChangesOfVariables]]
deps = ["InverseFunctions", "LinearAlgebra", "Test"]
git-tree-sha1 = "799b25ca3a8a24936ae7b5c52ad194685fc3e6ef"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.9"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "7eee164f122511d3e4e1ebadb7956939ea7e1c77"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.6"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "a33b7ced222c6165f624a3f2b55945fac5a598d9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.7"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.FunctionalCollections]]
deps = ["Test"]
git-tree-sha1 = "04cb9cfaa6ba5311973994fe3496ddec19b6292a"
uuid = "de31a74c-ac4f-5751-b3fd-e18cd04993ca"
version = "0.5.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "629693584cef594c3f6f99e76e7a7ad17e60e8d5"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a8863b69c2a0859f2c2c87ebdc4c6712e88bdf0d"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.7+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "401e4f3f30f43af2c8478fc008da50096ea5240f"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.3.1+0"

[[deps.Hiccup]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "6187bb2d5fcbb2007c39e7ac53308b0d371124bd"
uuid = "9fb69e20-1954-56bb-a84f-559cc56a8ff7"
version = "0.2.2"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Dates", "Test"]
git-tree-sha1 = "2787db24f4e03daf859c6509ff87764e4182f7d1"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.16"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "39d64b09147620f5ffbf6b2d3255be3c901bec63"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.8"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "f389674c99bfcde17dc57454011aa44d5a260a40"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.6.0"

[[deps.JSExpr]]
deps = ["JSON", "MacroTools", "Observables", "WebIO"]
git-tree-sha1 = "b413a73785b98474d8af24fd4c8a975e31df3658"
uuid = "97c1335a-c9c5-57fe-bc5d-ec35cebe8660"
version = "0.5.4"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "4b415b6cccb9ab61fec78a621572c82ac7fa5776"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.35"

[[deps.Kaleido_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "43032da5832754f58d14a91ffbe86d5f176acda9"
uuid = "f7e6163d-2fa5-5f23-b69c-1db539e41963"
version = "0.2.1+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e16271d212accd09d52ee0ae98956b8a05c4b626"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "17.0.6+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "ce5f5621cac23a86011836badfedf664a612cee4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.5"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "1ce1834f9644a8f7c011eb0592b7fd6c42c90653"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "3.0.1"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.Mustache]]
deps = ["Printf", "Tables"]
git-tree-sha1 = "3b2db451a872b20519ebb0cec759d3d81a1c6bcb"
uuid = "ffc61752-8dc7-55ee-8c37-f3e9cdd09e70"
version = "1.0.20"

[[deps.Mux]]
deps = ["AssetRegistry", "Base64", "HTTP", "Hiccup", "Pkg", "Sockets", "WebSockets"]
git-tree-sha1 = "82dfb2cead9895e10ee1b0ca37a01088456c4364"
uuid = "a975b10e-0019-58db-a62f-e48ff68538c9"
version = "0.7.6"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e127b609fb9ecba6f201ba7ab753d5a605d53801"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.54.1+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.PlotlyBase]]
deps = ["ColorSchemes", "Dates", "DelimitedFiles", "DocStringExtensions", "JSON", "LaTeXStrings", "Logging", "Parameters", "Pkg", "REPL", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "56baf69781fc5e61607c3e46227ab17f7040ffa2"
uuid = "a03496cd-edff-5a9b-9e67-9cda94a718b5"
version = "0.8.19"

[[deps.PlotlyJS]]
deps = ["Base64", "Blink", "DelimitedFiles", "JSExpr", "JSON", "Kaleido_jll", "Markdown", "Pkg", "PlotlyBase", "PlotlyKaleido", "REPL", "Reexport", "Requires", "WebIO"]
git-tree-sha1 = "e62d886d33b81c371c9d4e2f70663c0637f19459"
uuid = "f0f68f2c-4968-5e81-91da-67840de0976a"
version = "0.18.13"

[[deps.PlotlyKaleido]]
deps = ["Base64", "JSON", "Kaleido_jll"]
git-tree-sha1 = "2650cd8fb83f73394996d507b3411a7316f6f184"
uuid = "f2990250-8cf9-495f-b13a-cce12b45703c"
version = "2.2.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "45470145863035bb124ca51b320ed35d071cc6c2"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.8"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "5d9ab1a4faf25a62bb9d07ef0003396ac258ef1c"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.15"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "7b7850bb94f75762d567834d7e9802fc22d62f9c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.18"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIParser]]
deps = ["Unicode"]
git-tree-sha1 = "53a9f49546b8d2dd2e688d216421d050c9a31d0d"
uuid = "30578b45-9adc-5946-b283-645ec420af67"
version = "0.4.1"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.WebIO]]
deps = ["AssetRegistry", "Base64", "Distributed", "FunctionalCollections", "JSON", "Logging", "Observables", "Pkg", "Random", "Requires", "Sockets", "UUIDs", "WebSockets", "Widgets"]
git-tree-sha1 = "0eef0765186f7452e52236fa42ca8c9b3c11c6e3"
uuid = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"
version = "0.8.21"

[[deps.WebSockets]]
deps = ["Base64", "Dates", "HTTP", "Logging", "Sockets"]
git-tree-sha1 = "f91a602e25fe6b89afc93cf02a4ae18ee9384ce3"
uuid = "104b5d7c-a370-577a-8038-80a2059c5097"
version = "1.5.9"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "1165b0443d0eca63ac1e32b8c0eb69ed2f4f8127"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "936081b536ae4aa65415d869287d43ef3cb576b2"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.53.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─c3ebeb58-24c4-11ef-2fd7-a55d234fb215
# ╟─8733a11b-ab6c-43e3-b4c1-f126a3f5640c
# ╟─f2290f80-3ff5-42cf-8d3b-a152d8bbf578
# ╟─341beaad-0527-4f37-8666-71cfc80a4bae
# ╟─cb298b2f-c4fe-4045-b2fd-64c561ab9da5
# ╟─61a4d03e-cbf1-4ff2-a915-b68dddfacba5
# ╟─c9433bc7-274e-439b-bfa7-4da9d9b32480
# ╟─fd5de002-9e52-4483-81c8-edffa078a83a
# ╟─77501480-6a29-4e3d-9d71-2c2d91a0f6bd
# ╟─0f4b74d2-23f5-4061-beaa-e6abc44d7afb
# ╟─46dc4504-69c7-49c7-a217-804b3a651c84
# ╟─57673298-47f9-402b-a63b-0263adbb238a
# ╟─4dd8f7a1-f5cd-49ae-89ab-ce89da43ff65
# ╠═f5c03c4a-067d-432f-9205-c0693394faa4
# ╠═239595e8-89e2-4f1e-a503-ec30dbdf1a85
# ╟─3ae46ec8-1755-4dee-8c51-8406d0cb5e3c
# ╠═c3522930-99d2-4e61-a591-2fb2ae27b1e8
# ╠═ecd3a994-3c72-432d-bac2-b4fa884631cb
# ╟─4cf99d73-fc9b-468e-8b47-862a63a3f5a7
# ╠═583cef9d-16a4-424d-a42e-e115540b0ae8
# ╟─d35215d1-f8b3-42e7-b222-610763e0476a
# ╠═b726aaf2-e02b-4959-9058-b86ff6c98dd7
# ╠═218d9e15-d242-4625-ae9d-37f4408be3e9
# ╠═060ec34f-9653-46ef-9a6e-12f0cac61153
# ╟─7cc1227e-199d-4497-ab14-cf61c9d24834
# ╠═7644de2a-dbbd-475a-8926-29ee439579f2
# ╟─a6dc0326-67ea-4a2c-8bcd-bff2b97b3847
# ╠═8b4e721b-3a99-4070-a7b6-7999956afba2
# ╠═8c0155fe-493b-4221-8658-64e68491d8e7
# ╟─5b1c7301-2e94-460f-b798-ecb8611e6017
# ╟─befb8334-b8cc-46d5-bede-246612c9feae
# ╟─e33ca02f-6a79-480c-849d-32a32e3b9739
# ╠═90e035ef-eded-4d79-b652-5a4643df57a6
# ╟─b613c16e-66d2-49f7-8527-112335f97d58
# ╠═2dd13c57-59a3-4c72-9aad-eb7c603ae5e1
# ╟─ed7ecc33-9774-4c6d-86a6-7eb78befe270
# ╠═bfbc714e-ca3a-4a4e-a704-1342f477a5c7
# ╟─cbf17821-1b31-4558-8855-c5451e8e2cf5
# ╠═9da8113b-43fb-4439-9e5f-6e4d3bd69028
# ╟─5940de23-eca9-411b-b43c-30bc8f8f23d3
# ╠═a3a612b9-a454-4692-b85b-74d6017f7818
# ╟─1153a472-58f2-4417-a679-b9d77fb15937
# ╠═86808aa3-68ac-49a4-b7d4-df4fa95873b3
# ╠═8334b234-e6d6-486d-a7bb-8d4b4be1a5ff
# ╠═84a992bc-0166-4485-bdef-d6f6aa1b21fe
# ╟─8858623a-c619-4f4c-8e80-41b471ee5fa2
# ╠═fbea9fa5-2ce0-4db1-bc28-b5e7a85bdd57
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
