module particles_dynamics

importall Base

export update_particles_LeapFrog, update_particles_KDK_step_1, update_particles_KDK_step_2


function Scale_Function(  a, Omega_M, Omega_L,  Omega_K )
  a3 = a * a * a;
  factor = ( Omega_M + a*Omega_K + a3*Omega_L ) / a;
  return 1. / sqrt(factor);
end

function update_particles_LeapFrog( n_local, current_a, delta_a, p_pos_x, p_pos_y, p_pos_z, p_vel_x, p_vel_y, p_vel_z, p_grav_x,  p_grav_y, p_grav_z, omega_M, omega_L, omega_K, H0, cosmo_h, Lx, Ly, Lz, xMin, yMin, zMin, xMax, yMax, zMax  )

  scale_factor = Scale_Function( current_a, omega_M, omega_L, omega_K ) / H0 * cosmo_h;
  scale_factor_1 = Scale_Function( current_a + 0.5*delta_a, omega_M, omega_L, omega_K  ) / H0 * cosmo_h;
  a2_inv = 1./( ( current_a + 0.5*delta_a )*( current_a + 0.5*delta_a ));

  for i in 1:n_local
    p_vel_x[i] += scale_factor * delta_a * p_grav_x[i]
    p_vel_y[i] += scale_factor * delta_a * p_grav_y[i]
    p_vel_z[i] += scale_factor * delta_a * p_grav_z[i]

    p_pos_x[i] += a2_inv * scale_factor_1 * delta_a * p_vel_x[i]
    p_pos_y[i] += a2_inv * scale_factor_1 * delta_a * p_vel_y[i]
    p_pos_z[i] += a2_inv * scale_factor_1 * delta_a * p_vel_z[i]

    if p_pos_x[i] >= xMax
      p_pos_x[i] -= Lx
    end
    if p_pos_x[i] < xMin
      p_pos_x[i] += Lx
    end

    if p_pos_y[i] >= yMax
      p_pos_y[i] -= Ly
    end
    if p_pos_y[i] < yMin
      p_pos_y[i] += Ly
    end

    if p_pos_z[i] >= zMax
      p_pos_z[i] -= Lz
    end
    if p_pos_z[i] < zMin
      p_pos_z[i] += Lz
    end

  end

end



function update_particles_KDK_step_1( n_local, current_a, delta_a, p_pos_x, p_pos_y, p_pos_z, p_vel_x, p_vel_y, p_vel_z, p_grav_x,  p_grav_y, p_grav_z, omega_M, omega_L, omega_K, H0, cosmo_h, Lx, Ly, Lz, xMin, yMin, zMin, xMax, yMax, zMax  )

  scale_factor = Scale_Function( current_a, omega_M, omega_L, omega_K ) / H0 * cosmo_h;
  # scale_factor_1 = Scale_Function( current_a + 0.5*delta_a, omega_M, omega_L, omega_K  ) / H0 * cosmo_h;
  # a2_inv = 1./( ( current_a + 0.5*delta_a )*( current_a + 0.5*delta_a ));
  scale_factor_1 = Scale_Function( current_a , omega_M, omega_L, omega_K  ) / H0 * cosmo_h;
  a2_inv = 1./( ( current_a  )*( current_a  ));

  for i in 1:n_local
    p_vel_x[i] += scale_factor * 0.5 * delta_a * p_grav_x[i]
    p_vel_y[i] += scale_factor * 0.5 * delta_a * p_grav_y[i]
    p_vel_z[i] += scale_factor * 0.5 * delta_a * p_grav_z[i]

    p_pos_x[i] += a2_inv * scale_factor_1 * delta_a * p_vel_x[i]
    p_pos_y[i] += a2_inv * scale_factor_1 * delta_a * p_vel_y[i]
    p_pos_z[i] += a2_inv * scale_factor_1 * delta_a * p_vel_z[i]

    if p_pos_x[i] >= xMax
      p_pos_x[i] -= Lx
    end
    if p_pos_x[i] < xMin
      p_pos_x[i] += Lx
    end

    if p_pos_y[i] >= yMax
      p_pos_y[i] -= Ly
    end
    if p_pos_y[i] < yMin
      p_pos_y[i] += Ly
    end

    if p_pos_z[i] >= zMax
      p_pos_z[i] -= Lz
    end
    if p_pos_z[i] < zMin
      p_pos_z[i] += Lz
    end

  end

end

function update_particles_KDK_step_2( n_local, current_a, delta_a, p_pos_x, p_pos_y, p_pos_z, p_vel_x, p_vel_y, p_vel_z, p_grav_x,  p_grav_y, p_grav_z, omega_M, omega_L, omega_K, H0, cosmo_h, Lx, Ly, Lz, xMin, yMin, zMin, xMax, yMax, zMax  )

  scale_factor = Scale_Function( current_a + delta_a, omega_M, omega_L, omega_K ) / H0 * cosmo_h;
  # scale_factor_1 = Scale_Function( current_a + 0.5*delta_a, omega_M, omega_L, omega_K  ) / H0 * cosmo_h;
  # a2_inv = 1./( ( current_a + delta_a )*( current_a + 0.5*delta_a ));

  for i in 1:n_local
    p_vel_x[i] += scale_factor * 0.5 * delta_a * p_grav_x[i]
    p_vel_y[i] += scale_factor * 0.5 * delta_a * p_grav_y[i]
    p_vel_z[i] += scale_factor * 0.5 * delta_a * p_grav_z[i]

  end

end




end
