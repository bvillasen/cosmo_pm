using HDF5
current_dir = pwd()
# modules_dir = current_dir * "/modules"
push!(LOAD_PATH, current_dir)
using functions_CIC
using gravity_functions
using particles_dynamics
using io_functions

dataDir = "/home/bruno/Desktop/data/"
cosmoDir = dataDir * "cosmo_sims/cosmo_pm/cosmo_256_dm/"
outputDir = dataDir * "cosmo_sims/cosmo_pm/cosmo_256_dm/kdk/set_1/"


#Grid Parameters
const xMin = 0
const yMin = 0
const zMin = 0
const xMax = 115000
const yMax = 115000
const zMax = 115000
const Lx = xMax - xMin
const Ly = yMax - yMin
const Lz = zMax - zMin
const nx = 256
const ny = 256
const nz = 256
const dx = ( xMax - xMin ) / nx
const dy = ( yMax - yMin ) / ny
const dz = ( zMax - zMin ) / nz


# Parameters for cosmology
H = 67.74                        #[km/s / Mpc]
const cosmo_h = H/100
const H0 = H/1000                #[km/s / kpc]
const Omega_M = 0.3089
const Omega_L = 0.6911
const Omega_K = 0.0
# const cosmo_G = 4.29448e-06      #kpc km^2 s^-2 Msun^-1
const cosmo_G = 4.30087630e-06

const max_delta_a = 2e-4
delta_a = max_delta_a

const r_0_dm   = dx
const t_0_dm   = 1. / H0
const v_0_dm   = r_0_dm /t_0_dm / cosmo_h
const rho_0_dm = 3*H0*H0 / ( 8*pi*cosmo_G ) * Omega_M /cosmo_h/cosmo_h


icsFileName = cosmoDir * "ics/0_particles.h5"
icsFile = h5open( icsFileName, "r")

#Load Initial Conditions
println( "Loading Ics: $(icsFileName)")
const n_local = read( attrs(icsFile), "n_particles_local")
current_a = read( attrs(icsFile), "current_a")
current_z = read( attrs(icsFile), "current_z")
mass = read( icsFile, "mass")
pos_x = read( icsFile, "pos_x")
pos_y = read( icsFile, "pos_y")
pos_z = read( icsFile, "pos_z")
vel_x = read( icsFile, "vel_x") * current_a
vel_y = read( icsFile, "vel_y") * current_a
vel_z = read( icsFile, "vel_z") * current_a
close(icsFile)
grav_x = zeros( n_local )
grav_y = zeros( n_local )
grav_z = zeros( n_local )

#Load Initial Conditions
outputsFileName = "./scale_output_files/outputs_enzo_dm.txt"
println( "Loading Scale Outpus: $(outputsFileName)")
outputs = readdlm( outputsFileName )
n_outpus = length( outputs )
println( " Outputs Loaded: $(n_outpus)")


#Density CIC parameters
nGHST_CIC = 1
const nCells_x = nx + 2*nGHST_CIC
const nCells_y = ny + 2*nGHST_CIC
const nCells_z = nz + 2*nGHST_CIC

#Initialize memory for gravity
nGHST_POT = 2
potential = zeros( nz+2*nGHST_POT, ny+2*nGHST_POT, nx+2*nGHST_POT)
#Gravity Memory
nx_g = nx+2*nGHST_CIC
ny_g = ny+2*nGHST_CIC
nz_g = nz+2*nGHST_CIC
gravity_x = zeros( nz+2*nGHST_CIC, ny+2*nGHST_CIC, nx+2*nGHST_CIC )
gravity_y = zeros( nz+2*nGHST_CIC, ny+2*nGHST_CIC, nx+2*nGHST_CIC )
gravity_z = zeros( nz+2*nGHST_CIC, ny+2*nGHST_CIC, nx+2*nGHST_CIC )

# part_method = "LeapFrog"
part_method = "KDK"
println( "Particle Method: $(part_method)")


#FFT parameters
println( "Initializing FFT")
G, fft_plan_fwd, fft_plan_bkwd = initialize_FFT( nx, ny, nz, Lx, Ly, Lz, dx, dy, dz )

#Get density
print( "Computing Density CIC")
time_density = @elapsed density =  get_density_CIC( n_local, mass, pos_x, pos_y, pos_z, nCells_x, nCells_y, nCells_z, dx, dy, dz, xMin, yMin, zMin )
dens_avrg = mean( density )
print( "   $(time_density)    $(dens_avrg)\n")

if part_method == "KDK"
  #Get Potential
  print( "Computing Potential")
  time_potential = @elapsed  get_potential( nx, ny, nz, nGHST_POT, cosmo_G, density, G, fft_plan_fwd, fft_plan_bkwd, potential, dens_avrg, current_a )
  print( "   $(time_potential)\n")

  #Get Gravitational Field
  print( "Computing Gravity Field")
  time_grav = @elapsed  get_gravity( nx, ny, nz, nGHST_CIC, nGHST_POT, dx, dy, dz, potential, gravity_x, gravity_y, gravity_z )
  print( "   $(time_grav)\n")

  #Get gravity CIC accelerations
  print( "Computing Gravity CIC")
  time_grav_CIC = @elapsed get_gravity_CIC( n_local, pos_x, pos_y, pos_z, xMin, yMin, zMin, dx, dy, dz, gravity_x, gravity_y, gravity_z, grav_x, grav_y, grav_z  )
  print( "   $(time_grav_CIC)\n")
end


#Save Snapshot
n_file = 0
next_output, n_file = write_snapshot(  n_file,  outputDir, current_a, n_local, density,  pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, outputs )

println( "\nStarting Simulation")
tic()
n_step = 0
while current_a < outputs[n_outpus]
  n_step += 1
  delta_a = max_delta_a

  if ( (current_a + delta_a) >= next_output )
    delta_a = (next_output - current_a)
  end

  println( "\nStep: $(n_step)   current_a: $(current_a)   delta_a: $(delta_a)"   )


  if part_method == "KDK"
    print( " Updating Particles step 1")
    time_particles_update = @elapsed update_particles_KDK_step_1( n_local, current_a, delta_a, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, grav_x, grav_y, grav_z, Omega_M, Omega_L, Omega_K, H0, cosmo_h , Lx, Ly, Lz, xMin, yMin, zMin, xMax, yMax, zMax )
    print( "   $(time_particles_update)\n")
  end


  #Get density
  print( " Computing Density CIC")
  time_density = @elapsed density =  get_density_CIC( n_local, mass, pos_x, pos_y, pos_z, nCells_x, nCells_y, nCells_z, dx, dy, dz, xMin, yMin, zMin )
  dens_avrg = mean( density )
  print( "   $(time_density)    $(dens_avrg)\n")

  #Get Potential
  print( " Computing Potential")
  time_potential = @elapsed  get_potential( nx, ny, nz, nGHST_POT, cosmo_G, density, G, fft_plan_fwd, fft_plan_bkwd, potential, dens_avrg, current_a )
  print( "   $(time_potential)\n")

  #Get Gravitational Field
  print( " Computing Gravity Field")
  time_grav = @elapsed  get_gravity( nx, ny, nz, nGHST_CIC, nGHST_POT, dx, dy, dz, potential, gravity_x, gravity_y, gravity_z )
  print( "   $(time_grav)\n")

  #Get gravity CIC accelerations
  print( " Computing Gravity CIC")
  time_grav_CIC = @elapsed get_gravity_CIC( n_local, pos_x, pos_y, pos_z, xMin, yMin, zMin, dx, dy, dz, gravity_x, gravity_y, gravity_z, grav_x, grav_y, grav_z  )
  print( "   $(time_grav_CIC)\n")

  if part_method == "LeapFrog"
    print( " Updating Particles")
    time_particles_update = @elapsed update_particles_LeapFrog( n_local, current_a, delta_a, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, grav_x, grav_y, grav_z, Omega_M, Omega_L, Omega_K, H0, cosmo_h , Lx, Ly, Lz, xMin, yMin, zMin, xMax, yMax, zMax )
    print( "   $(time_particles_update)\n")
  end

  if part_method == "KDK"
    print( " Updating Particles step 2")
    time_particles_update = @elapsed update_particles_KDK_step_2( n_local, current_a, delta_a, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, grav_x, grav_y, grav_z, Omega_M, Omega_L, Omega_K, H0, cosmo_h , Lx, Ly, Lz, xMin, yMin, zMin, xMax, yMax, zMax )
    print( "   $(time_particles_update)\n")
  end

  current_a += delta_a

  if abs(current_a - next_output) < 1e-7
    next_output, n_file = write_snapshot(  n_file,  outputDir, current_a, n_local, density,  pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, outputs )
  end

end

println( "\nSimulation Finished")
toc()














# for k in 0:nz-1
#   for j in 0:ny-1
#     for i in 0:nx-1
#       id = i + j*nx + k*nx*ny
#       p_x = xMin + i*dx + 0.5*dx + 0.1*dx
#       p_y = yMin + j*dy + 0.5*dy - 0.3*dy
#       p_z = zMin + k*dz + 0.5*dz + 0.2*dz
#       pos_x[id+1] = p_x
#       pos_y[id+1] = p_y
#       pos_z[id+1] = p_z
#       mass[id+1] = dx*dy*dz
#     end
#   end
# end
