module gravity_functions

importall Base

export get_potential, initialize_FFT, get_gravity


function get_potential( nx, ny, nz, nGHST, Gconst, density, G, fft_plan_f, fft_plan_b, potential, dens_avrg, current_a )
  # Apply FFT to rho
  # 4 * M_PI * G.Cosmo.cosmo_G * ( G.Particles.G.density[id_CIC] - G.Cosmo.rho_0_dm ) / G.Cosmo.current_a ;
  rho_trans = fft_plan_f * (4*pi*Gconst * (density -dens_avrg) ) /current_a
  rho_trans = G .* rho_trans
  rho_trans[1] = 0
  # Apply inverse FFTW
  phi = real( fft_plan_b * rho_trans )
  # potential = zeros( nz+2*nGHST, ny+2*nGHST, nx+2*nGHST )
  potential[nGHST+1:end-nGHST, nGHST+1:end-nGHST, nGHST+1:end-nGHST] = phi
  potential[1, :, : ] = potential[end-2*nGHST+1,:,:]
  potential[2, :, : ] = potential[end-2*nGHST+2,:,:]
  potential[:, 1, : ] = potential[:,end-2*nGHST+1,:]
  potential[:, 2, : ] = potential[:,end-2*nGHST+2,:]
  potential[:, :, 1 ] = potential[:,:,end-2*nGHST+1]
  potential[:, :, 2 ] = potential[:,:,end-2*nGHST+2]
  potential[end,:,:] = potential[nGHST+2,:,:]
  potential[end-1,:,:] = potential[nGHST+1,:,:]
  potential[:,end,:] = potential[:,nGHST+2,:]
  potential[:,end-1,:] = potential[:,nGHST+1,:]
  potential[:,:,end] = potential[:,:,nGHST+2]
  potential[:,:,end-1] = potential[:,:,nGHST+1]
  return potential
end

function initialize_FFT( nCells_x, nCells_y, nCells_z, Lx, Ly, Lz, dx, dy, dz )
  dens = zeros( nCells_z, nCells_y, nCells_x )
  # Initialize FFTW
  FFTW.set_num_threads(8)
  fft_plan_fwd  = plan_fft( dens, flags=FFTW.ESTIMATE, timelimit=20 )
  fft_plan_bkwd = plan_ifft( dens, flags=FFTW.ESTIMATE, timelimit=20 )

  fft_kx = sin.( pi/Lx*dx * linspace(0, nCells_x-1, nCells_x) ).^2
  fft_ky = sin.( pi/Ly*dy * linspace(0, nCells_y-1, nCells_y) ).^2
  fft_kz = sin.( pi/Lz*dz * linspace(0, nCells_z-1, nCells_z) ).^2
  G = zeros( nCells_z, nCells_y, nCells_x )
  for i in 1:nCells_x
    sin_kx = fft_kx[i]
    for j in 1:nCells_y
      sin_ky = fft_ky[j]
      for k in 1:nCells_z
        sin_kz =  fft_kz[k]
        G[k,j,i] = -1/ ( sin_kx + sin_ky + sin_kz )  *dx^2/4  #NOTE: DX^2 only for dx=dy=dz
      end
    end
  end
  return [ G, fft_plan_fwd, fft_plan_bkwd ]
end


function get_gravity( nx, ny, nz, nGHST_CIC, nGHST_POT, dx, dy, dz, potential, grav_x, grav_y, grav_z )
  nx_g = nx+2*nGHST_CIC
  ny_g = ny+2*nGHST_CIC
  nz_g = nz+2*nGHST_CIC
  # nx_p = nx+2*nGHST_POT
  # ny_p = ny+2*nGHST_POT
  # nz_p = nz+2*nGHST_POT
  for i in 1:nx_g
    for j in 1:ny_g
      for k in 1:nz_g
        i_p, j_p, k_p = i+1, j+1, k+1
        pot_l = potential[ k_p, j_p, i_p-1 ]
        pot_r = potential[ k_p, j_p, i_p+1 ]
        grav_x[k, j, i] = -0.5 * ( pot_r - pot_l ) / dx
      end
    end
  end

  for i in 1:nx_g
    for j in 1:ny_g
      for k in 1:nz_g
        i_p, j_p, k_p = i+1, j+1, k+1
        pot_l = potential[ k_p, j_p-1, i_p ]
        pot_r = potential[ k_p, j_p+1, i_p ]
        grav_y[k, j, i] = -0.5 * ( pot_r - pot_l ) / dy
      end
    end
  end

  for i in 1:nx_g
    for j in 1:ny_g
      for k in 1:nz_g
        i_p, j_p, k_p = i+1, j+1, k+1
        pot_l = potential[ k_p-1, j_p, i_p ]
        pot_r = potential[ k_p+1, j_p, i_p ]
        grav_z[k, j, i] = -0.5 * ( pot_r - pot_l ) / dz
      end
    end
  end

end



end
