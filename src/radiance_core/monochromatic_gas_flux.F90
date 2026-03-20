! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate fluxes including only gaseous absorption.
!
! Method:
!   Transmission coefficients for each layer are calculated
!   from the gaseous absorption alone. fluxes are propagated
!   upward or downward through the column using these
!   coefficients and source terms.
!
!- ---------------------------------------------------------------------
MODULE monochromatic_gas_flux_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'MONOCHROMATIC_GAS_FLUX_MOD'
CONTAINS
SUBROUTINE monochromatic_gas_flux(n_profile, n_layer                    &
    , tau_gas                                                           &
    , isolir, sec_0, flux_inc_direct, flux_inc_down                     &
    , diff_planck, d_planck_flux_surface                                &
    , diffuse_albedo, direct_albedo                                     &
    , diffusivity_factor                                                &
    , flux_direct, flux_diffuse                                         &
    , nd_profile, nd_layer                                              &
    , source, trans &
    )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_infra_red, ip_solar
  USE vectlib_mod, ONLY : exp_v
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer
!       Maximum number of layers


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , isolir
!       Spectral region
  REAL (RealK), INTENT(IN) ::                                           &
      tau_gas(nd_profile, nd_layer)                                     &
!       Gaseous optical depths
    , sec_0(nd_profile)                                                 &
!       Secant of zenith angle
    , flux_inc_direct(nd_profile)                                       &
!       Incident direct flux
    , flux_inc_down(nd_profile)                                         &
!       Incident diffuse flux
    , d_planck_flux_surface(nd_profile)                                 &
!       Difference in Planckian fluxes between the surface
!       and the overlying air
    , diff_planck(nd_profile, nd_layer)                                 &
!       Difference in Planckian function
    , diffuse_albedo(nd_profile)                                        &
!       Diffuse surface albedo
    , direct_albedo(nd_profile)                                         &
!       Direct surface albedo
    , diffusivity_factor
!       Diffusivity factor
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_profile, 0: nd_layer)                              &
!       Direct flux
    , flux_diffuse(nd_profile, 2*nd_layer+2)
!       Diffuse flux

! Work arrays
  REAL (RealK) :: &
      trans(:, :) &
!       Transmissivities
    , source(:, :)
!       Source function


! Local variables.
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , l
!       Loop variable
! Variables related to the treatment of ill-conditioning
  REAL (RealK) ::                                                       &
      eps_r                                                             &
!       The smallest real number such that 1.0-EPS_R is not 1
!       to the computer's precision
    , sq_eps_r
!       The square root of the above

  REAL (RealK) :: temp

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='MONOCHROMATIC_GAS_FLUX'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the tolerances used in avoiding ill-conditioning, testing
! on any variable.
  eps_r=EPSILON(tau_gas(1, 1))
  sq_eps_r=SQRT(eps_r)

  IF (isolir == ip_infra_red) THEN

    !$omp target teams distribute parallel do simd
    DO l=1, n_profile

      DO i=1, n_layer
        trans(l, i)=exp(-diffusivity_factor*tau_gas(l, i))
        source(l, i)=(1.0e+00_RealK-trans(l, i)+sq_eps_r)            &
          *diff_planck(l, i)                                            &
          /(diffusivity_factor*tau_gas(l, i)+sq_eps_r)
      END DO

      ! Downward fluxes.
      flux_diffuse(l, 2)=flux_inc_down(l)
      DO i=1, n_layer
        flux_diffuse(l, 2*i+2)=trans(l, i)*flux_diffuse(l, 2*i)           &
          -source(l, i)
      END DO

      ! Upward fluxes.
      flux_diffuse(l, 2*n_layer+1)=d_planck_flux_surface(l)             &
        +diffuse_albedo(l)*flux_diffuse(l, 2*n_layer+2)
      DO i=n_layer, 1, -1
        flux_diffuse(l, 2*i-1)=trans(l, i)*flux_diffuse(l, 2*i+1)         &
          +source(l, i)
      END DO

    END DO

  ELSE IF (isolir == ip_solar) THEN
    STOP __LINE__

    DO l=1, n_profile

      DO i=1, n_layer
        trans(l, i)= exp(-diffusivity_factor*tau_gas(l, i))
      END DO

      ! The direct flux.
      flux_direct(l, 0)=flux_inc_direct(l)
      DO i=1, n_layer
        temp=exp(-tau_gas(l, i)*sec_0(l))
        flux_direct(l, i)=flux_direct(l, i-1)*temp
      END DO

      ! Downward fluxes.
      flux_diffuse(l, 2)=flux_inc_down(l)
      DO i=1, n_layer
        flux_diffuse(l, 2*i+2)=trans(l, i)*flux_diffuse(l, 2*i) 
      END DO

      ! Upward fluxes.
      flux_diffuse(l, 2*n_layer+1)=                                     &
        +diffuse_albedo(l)*flux_diffuse(l, 2*n_layer+2)                 &
        +direct_albedo(l)*flux_direct(l, n_layer)
      DO i=n_layer, 1, -1
        flux_diffuse(l, 2*i-1)=trans(l, i)*flux_diffuse(l, 2*i+1) 
      END DO
    END DO

  ELSE
    STOP __LINE__
  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE monochromatic_gas_flux
END MODULE monochromatic_gas_flux_mod
