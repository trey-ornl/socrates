! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate differences in source functions.
!
! Method:
!   The Planck fluxes for the band at the boundaries of layers
!   are found and differences across layers are determined. If
!   the Planckian is being taken to vary quadratically across
!   the layer second differences are found.
!
!------------------------------------------------------------------------------
MODULE diff_planck_source_mod

IMPLICIT NONE
PUBLIC :: diff_planck_source
PRIVATE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'DIFF_PLANCK_SOURCE_MOD'
CONTAINS

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Evaluate Planck flux for a given spectral band
!
!------------------------------------------------------------------------------

REAL(RealK) FUNCTION planck_flux_band(spectrum, i_band, temperature)

  USE def_spectrum, ONLY: StrSpecData
  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE
  !$omp declare target

  TYPE (StrSpecData), INTENT(IN) :: spectrum
  INTEGER, INTENT(IN) :: i_band
  REAL(RealK), INTENT(IN) :: temperature

  INTEGER :: j
  REAL(RealK) :: t_ratio

  
  t_ratio=temperature/spectrum%planck%t_ref_planck

  !STOP __LINE__
  ! Use the polynomial fit to the Planck function
  planck_flux_band &
    = spectrum%planck%thermal_coeff(spectrum%planck%n_deg_fit, i_band)
  DO j=spectrum%planck%n_deg_fit-1, 0, -1
    planck_flux_band = planck_flux_band * t_ratio &
      + spectrum%planck%thermal_coeff(j, i_band)
  END DO

END FUNCTION planck_flux_band

SUBROUTINE diff_planck_source(control, dimen, spectrum, atm, bound, &
  i_band, i_rad_layer, frac_rad_layer, planck)

  USE def_control,  ONLY: StrCtrl
  USE def_dimen,    ONLY: StrDim
  USE def_spectrum, ONLY: StrSpecData
  USE def_atm,      ONLY: StrAtm
  USE def_bound,    ONLY: StrBound
  USE def_planck,   ONLY: StrPlanck
  USE realtype_rd,  ONLY: RealK
  USE rad_pcf,      ONLY: ip_spherical_harmonic, ip_surf_alb_diff
  USE rad_ccf,      ONLY: pi
  USE yomhook,      ONLY: lhook, dr_hook
  USE parkind1,     ONLY: jprb, jpim

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),     INTENT(IN)  :: control

! Dimensions:
  TYPE(StrDim),      INTENT(IN)  :: dimen

! Spectral data:
  TYPE(StrSpecData), INTENT(IN)  :: spectrum

! Atmospheric properties:
  TYPE(StrAtm),      INTENT(IN)  :: atm

! Boundary conditions:
  TYPE(StrBound),    INTENT(IN)  :: bound

  INTEGER, INTENT(IN) :: &
    i_band, &
!     Spectral band
    i_rad_layer(dimen%nd_viewing_level)
!     Layers in which to intercept radiances
  REAL (RealK), INTENT(IN) :: &
    frac_rad_layer(dimen%nd_viewing_level)
!     Fractions below the tops of the layers

! Planckian emission fields
  TYPE(StrPlanck), INTENT(INOUT) :: planck


! Local variables.
  INTEGER :: i, k, l, n_layer, n_profile
!     Loop variables
  REAL (RealK) :: flux_i, flux_im1
  REAL (RealK), PARAMETER :: frac_tile_tol = TINY(1.0_RealK)
!     Tolerance for tile fraction when setting Planckian to zero

  LOGICAL :: l_flux_tile = .FALSE.
!     Local flag to use tile emission rather than temperature

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='DIFF_PLANCK_SOURCE'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF (control%i_angular_integration == ip_spherical_harmonic) THEN
    STOP __LINE__
!   Calculate the Planckian radiance on viewing levels.
    DO i=1, atm%n_viewing_level
      DO l=1, atm%n_profile
        planck%radiance(l, i) = planck_flux_band( spectrum, i_band, &
          atm%t_level(l, i_rad_layer(i)-1) &
          +(atm%t_level(l, i_rad_layer(i))-atm%t_level(l, i_rad_layer(i)-1)) &
          *frac_rad_layer(i) )
        planck%radiance(l, i)=planck%radiance(l, i)/pi
      END DO
    END DO
  END IF

  IF (spectrum%planck%l_planck_tbl) THEN
    STOP __LINE__
  END IF

  n_layer = atm%n_layer
  n_profile = atm%n_profile

! Calculate the change in the Planckian flux across each layer.
  !$omp target teams distribute parallel do simd collapse(2) &
  !$omp& map(to: spectrum%planck%n_deg_fit, spectrum%planck%t_ref_planck) &
  !$omp& map(to: atm%t_level, spectrum%planck%thermal_coeff) &
  !$omp& map(from: planck%diff, planck%flux) &
  !$omp& private(flux_i, flux_im1)
  DO i=0, n_layer
    DO l=1, n_profile
      flux_i = planck_flux_band(spectrum, i_band, atm%t_level(l,i))
      planck%flux(l,i) = flux_i
      IF (i > 0) THEN
        flux_im1 = planck_flux_band(spectrum, i_band, atm%t_level(l,i-1))
        planck%diff(l,i) = flux_i - flux_im1
      END IF
    END DO
  END DO

  IF (control%l_ir_source_quad) THEN
    STOP __LINE__
    ! Calculate the second difference if required.
    DO i=1, atm%n_layer
      DO l=1, atm%n_profile
        planck%diff_2(l, i) = planck_flux_band(spectrum, i_band, atm%t(l, i))
        planck%diff_2(l, i) = 2.0e+00_RealK*(planck%flux(l, i) &
          + planck%flux(l, i-1) - 2.0e+00_RealK * planck%diff_2(l, i))
      END DO
    END DO
  END IF

! Planckian flux at the surface.
! Calculated from surface emission or skin temperature.
  IF (control%l_flux_ground) THEN
    STOP __LINE__
    DO l=1, atm%n_profile
      planck%flux_ground(l) &
        = bound%flux_ground(l, i_band) &
        /(1.0_RealK - bound%rho_alb(l, ip_surf_alb_diff, i_band))
    END DO
  ELSE
    !STOP __LINE__
    !$omp target teams distribute parallel do simd &
    !$omp& map(to: spectrum%planck%n_deg_fit, spectrum%planck%t_ref_planck) &
    !$omp& map(to: bound%t_ground, spectrum%planck%thermal_coeff) &
    !$omp& map(from: planck%flux_ground)
    DO l=1, n_profile
      planck%flux_ground(l) &
        = planck_flux_band(spectrum, i_band, bound%t_ground(l))
    END DO
  END IF

! Local Planckian fluxes will be required on tiled surfaces.
! Furthermore, the overall Planckian will be calculated as a
! weighted sum of the individual components: this allows for
! variations in the Planckian between spectral bands more
! satisfactorily than the use of an equivalent temperature can.
  IF (control%l_tile) THEN

    STOP __LINE__
    ! Fluxes for tiles calculated from surface emission or skin temperature.
    DO k=1, bound%n_tile
      IF (allocated(control%l_flux_tile)) l_flux_tile = control%l_flux_tile(k)
      DO l=1, bound%n_point_tile
        IF (l_flux_tile .AND. bound%frac_tile(l, k) > frac_tile_tol) THEN
          planck%flux_tile(l, k) &
            = bound%flux_tile(l, k, i_band) &
            /(1.0_RealK - bound%rho_alb_tile(l, ip_surf_alb_diff, k, i_band))
        ELSE IF (bound%frac_tile(l, k) > frac_tile_tol) THEN
          planck%flux_tile(l, k) &
            = planck_flux_band(spectrum, i_band, bound%t_tile(l, k))
        ELSE
          planck%flux_tile(l, k) = 0.0_RealK
        END IF
      END DO
    END DO

    IF (control%l_tile_emissivity) THEN
      DO l=1, bound%n_point_tile
        planck%flux_ground(bound%list_tile(l)) &
          =bound%frac_tile(l, 1)*planck%flux_tile(l, 1) &
          *(1.0_RealK - bound%rho_alb_tile(l, ip_surf_alb_diff, 1, i_band))
      END DO
      DO k=2, bound%n_tile
        DO l=1, bound%n_point_tile
          planck%flux_ground(bound%list_tile(l)) &
            =planck%flux_ground(bound%list_tile(l)) &
            +bound%frac_tile(l, k)*planck%flux_tile(l, k) &
            *(1.0_RealK - bound%rho_alb_tile(l, ip_surf_alb_diff, k, i_band))
        END DO
      END DO
      DO l=1, bound%n_point_tile
        planck%flux_ground(bound%list_tile(l)) &
          =planck%flux_ground(bound%list_tile(l)) / ( 1.0_RealK &
          -bound%rho_alb(bound%list_tile(l), ip_surf_alb_diff, i_band) )
      END DO
    ELSE
      DO l=1, bound%n_point_tile
        planck%flux_ground(bound%list_tile(l)) &
          =bound%frac_tile(l, 1)*planck%flux_tile(l, 1)
      END DO
      DO k=2, bound%n_tile
        DO l=1, bound%n_point_tile
          planck%flux_ground(bound%list_tile(l)) &
            =planck%flux_ground(bound%list_tile(l)) &
            +bound%frac_tile(l, k)*planck%flux_tile(l, k)
        END DO
      END DO
    END IF

  END IF


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE diff_planck_source
END MODULE diff_planck_source_mod
