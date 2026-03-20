! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to calculate cloudy two-stream coefficients.
!
! Method:
!   The coeffients for each type of cloud are determined and
!   averaged.
!
!- ---------------------------------------------------------------------
MODULE two_coeff_cloud_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'TWO_COEFF_CLOUD_MOD'
PRIVATE
PUBLIC :: two_coeff_cloud
CONTAINS
SUBROUTINE trans_source_coeff_scalar(l_ip_ir, l_ip_solar &
     , l_ir_source_quad, l_no_or_csr, l_spherical_solar &
     , tau_dir, tau, summ, diff, lambda, sec_0 & !, path_div &
     , gamma_up, gamma_down &
     , trans, reflect, trans_0_dir, trans_0 &
     , source_solar_up, source_solar_down &
     , source_ir_1d, source_ir_2d &
     )

  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE

! Control options:
  LOGICAL, INTENT(IN) :: l_ip_ir, l_ip_solar, l_ir_source_quad, &
    l_no_or_csr, l_spherical_solar

! Optical properties of the layer
  REAL (RealK), INTENT(IN) ::                                           &
      tau &
!       Optical depths of layers
    , tau_dir &
!       Unscaled optical depths of layers
    , summ &
!       Sum of alpha_1 and alpha_2
    , diff &
!       Difference of alpha_1 and alpha_2
    , lambda &
!       Lambda
    , sec_0 &
!       Secant of solar zenith angle
    !, path_div &
!       Path scaling for spherical geometry
    , gamma_up &
!       Basic solar coefficient for upward radiation
    , gamma_down
!       Basic solar coefficient for downward radiation

! Transmission and reflection coefficients and coefficients for
! source terms.
  REAL (RealK), INTENT(OUT) ::                                          &
      trans &
!       Diffuse transmission coefficient
    , reflect &
!       Diffuse reflection coefficient
    , trans_0 &
!       Direct transmission coefficient
    , trans_0_dir &
!       Direct transmission coefficient without scaling
    , source_solar_down, source_solar_up, source_ir_1d, source_ir_2d
!       Source coefficients


! Local variables
  REAL (RealK) ::                                                       &
      gamma1                                                             &
!       Gamma
    , exponential                                                       &
!       Exponential of scaled optical depth
    , gamma2                                                            &
!       Gamma squared
    , exponential2
!       Exponential squared

! Variables related to the treatment of ill-conditioning
  REAL (RealK), PARAMETER :: sq_eps_r=SQRT(EPSILON(sq_eps_r))
!       The square root of a real number that is negligible
!       compared to 1.0.
  REAL (RealK), PARAMETER :: tol=SQRT(sq_eps_r)
!       The tolerance used for switching to the asymptotic form
!       of the quadratic source function term.
  REAL (RealK) :: tmp_inv                 
!       Temporary work variable

! Determine the diffuse transmission and reflection coefficients.
  exponential=exp(-lambda*tau)
  exponential2=exponential*exponential
  gamma1=(summ-lambda) / (summ+lambda)
  gamma2=gamma1*gamma1
  tmp_inv=1.0_RealK / ( 1.0_RealK - exponential2*gamma2 )
  trans=exponential*(1.0_RealK-gamma2)*tmp_inv
  reflect=gamma1*(1.0_RealK-exponential2)*tmp_inv

  trans_0=0.0_RealK
  trans_0_dir=0.0_RealK
  source_ir_1d=0.0_RealK
  source_ir_2d=0.0_RealK
  source_solar_down=0.0_RealK
  source_solar_up=0.0_RealK

  IF (l_ip_solar) THEN

    !   Calculate the direct transmission and the source coefficients
    !   for the solar beam: in the solar case these are
    !   the coefficients which will multiply the direct flux at the
    !   top of the layer to give the source terms for the upward
    !   diffuse flux and the total downward flux.

    IF (l_spherical_solar) THEN
#if 0
      ! Untested
      trans_0 = exp(-tau*path_div)
      source_solar_up                           &
        = gamma_up-reflect*gamma_down             &
        - gamma_up*trans*trans_0
      source_solar_down=trans_0           &
        *(gamma_down-gamma_up*reflect)            &
        - gamma_down*trans
#endif
    ELSE
      trans_0 = exp(-tau*sec_0)
      IF (l_no_or_csr) THEN
        trans_0_dir = exp(-tau_dir*sec_0)
      END IF
      source_solar_up                           &
        =(gamma_up-reflect                              &
        *(1.0_RealK+gamma_down))                              &
        -gamma_up*trans*trans_0
      source_solar_down=trans_0           &
        *(1.0_RealK+gamma_down                                &
        -gamma_up*reflect)                              &
        -(1.0_RealK+gamma_down)*trans
    END IF

  ELSE IF (l_ip_ir) THEN

    !   In the case of infra-red radiation, the first source
    !   coefficient holds the multiplier for the first difference
    !   of the Planckian function across the layer, and the second
    !   that for the second difference.

    !       A tolerance is added to the numerator and the denomiator
    !       to avoid ill-conditioning at small optical depths.

    source_ir_1d=(1.0_RealK-trans         &
      +reflect+sq_eps_r)                                      &
      /(sq_eps_r+tau*summ)

    IF (l_ir_source_quad) THEN

      !     Quadratic correction to source function.
      !     This correction is very ill-conditioned for
      !     small optical depths so the asymptotic form is then used.

      ! Untested
      IF (tau > tol) THEN
        source_ir_2d                            &
          =-2.0_RealK                                               &
          *(1.0_RealK-trans-reflect+sq_eps_r)           &
          /(diff*tau+sq_eps_r)
      ELSE
        source_ir_2d                            &
          =-2.0_RealK+diff*tau
      END IF
      source_ir_2d                              &
        =-(1.0_RealK+reflect+trans                      &
        +source_ir_2d)                          &
        /(summ*tau+sq_eps_r)

    END IF

  END IF

END SUBROUTINE trans_source_coeff_scalar

SUBROUTINE solar_coefficient_basic_scalar(l_discrete_ord &
    , l_spherical_solar &
    , omega, asymmetry, sec_0 &
    !, path_div                                 &
    , summ, diff, lambda                                                &
    , gamma_up, gamma_down                                              &
    )

  USE realtype_rd, ONLY: RealK

  IMPLICIT NONE

! Control options:
  LOGICAL, INTENT(IN) :: l_discrete_ord, l_spherical_solar

  REAL (RealK), INTENT(IN) ::                                           &
      omega &
!       Albedo of single scattering
    , asymmetry &
!       Asymmetry
    , sec_0 !&
!       Secant of solar zenith angle
    !, path_div
!       Scaled path for spherical geometry

! Basic two-stream coefficients:
  REAL (RealK), INTENT(INOUT) ::                                        &
      summ &
!       Sum of two-stream coefficients
    , diff &
!       Difference of two-stream coefficients
    , lambda
!       Lambda
  REAL (RealK), INTENT(OUT) ::                                          &
      gamma_up &
!       Coefficient for upward radiation
    , gamma_down
!       Coefficient for downwad radiation

  REAL (RealK) ::                                                       &
      ksi_0 &
!       Difference in solar scattering fractions
    , factor
!       Temporary variable
  REAL (RealK) ::                                                       &
      root_3
!       Square root of 3
  PARAMETER(                                                            &
      root_3=1.7320508075688772e+00_RealK                               &
    )

! Variables related to the treatment of ill-conditioning
  REAL (RealK), PARAMETER :: &
    tol_perturb=3.2e+01_RealK*EPSILON(tol_perturb)
!       The tolerance used to judge where the two-stream
!       expressions for the solar source become ill-conditioned

! Determine the basic solar coefficients for the two-stream equations.
  IF (l_spherical_solar) THEN
#if 0
    ! Untested
    ! If LAMBDA is too close to SEC_0 it must be perturbed.
    IF (ABS(lambda-path_div) < tol_perturb) THEN
      summ=(1.0e+00_RealK+tol_perturb)*summ
      diff=(1.0e+00_RealK+tol_perturb)*diff
      lambda=(1.0e+00_RealK+tol_perturb)*lambda
    END IF
    factor=0.5e+00_RealK*omega*path_div                 &
      /((lambda-path_div)*(lambda+path_div))
    IF (l_discrete_ord) THEN
      ksi_0=root_3*asymmetry/path_div
    ELSE
      ksi_0=1.5e+00_RealK*asymmetry/path_div
    END IF
    gamma_up=factor*(summ-path_div                 &
      -ksi_0*(diff-path_div))
    gamma_down=factor*(summ+path_div               &
      +ksi_0*(diff+path_div))
#endif
  ELSE
    IF ((ABS(lambda-sec_0)) < tol_perturb) THEN
      summ=(1.0e+00_RealK+tol_perturb)*summ
      diff=(1.0e+00_RealK+tol_perturb)*diff
      lambda=(1.0e+00_RealK+tol_perturb)*lambda
    END IF
    factor=0.5e+00_RealK*omega*sec_0                       &
      /((lambda-sec_0)*(lambda+sec_0))
    IF (l_discrete_ord) THEN
      ! not tested
      ksi_0=root_3*asymmetry/sec_0
    ELSE
      ksi_0=1.5e+00_RealK*asymmetry/sec_0
    END IF
    gamma_up=factor*(summ-sec_0                       &
      -ksi_0*(diff-sec_0))
    gamma_down=factor*(summ+sec_0                     &
      +ksi_0*(diff+sec_0))
  END IF

END SUBROUTINE solar_coefficient_basic_scalar

SUBROUTINE two_coeff_basic_scalar( &
       i_2stream                                                        &
     , asymmetry, omega                                                 &
     , summ, diff                                                       &
     )

  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLY: ip_eddington, ip_elsasser, ip_discrete_ord, &
                     ip_pifm85, ip_2s_test, ip_hemi_mean, ip_pifm80
  USE diffusivity_factor, ONLY: elsasser_factor

  IMPLICIT NONE
  !$omp declare target

  INTEGER, INTENT(IN) ::                                                &
      i_2stream
!       Two stream scheme

! Optical properties of layer:
  REAL (RealK), INTENT(IN) ::                                           &
      asymmetry &
!       Asymmetry factor
    , omega
!       Albedo of single scattering

! coefficients in the two-stream equations:
  REAL (RealK), INTENT(OUT) ::                                          &
      summ &
!       Sum of alpha_1 and alpha_2
    , diff
!       Difference of alpha_1 and alpha_2

  REAL (RealK), PARAMETER ::                                            &
      root_3 = 1.7320508075688772e+00_RealK
!       Square root of 3

  IF (i_2stream == ip_eddington) THEN
    ! Untested
    summ=1.5e+00_RealK*(1.0e+00_RealK                          &
      -omega*asymmetry)
    diff=2.0e+00_RealK*(1.0e+00_RealK-omega)
  ELSE IF (i_2stream == ip_elsasser) THEN
    summ=elsasser_factor                                       &
      -1.5e+00_RealK*omega*asymmetry
    diff=elsasser_factor*(1.0e+00_RealK-omega)
  ELSE IF (i_2stream == ip_discrete_ord) THEN
    ! Untested
    summ=root_3*(1.0e+00_RealK                                 &
      -omega*asymmetry)
    diff=root_3*(1.0e+00_RealK-omega)
  ELSE IF (i_2stream == ip_pifm85) THEN
    ! Untested
    summ=2.0e+00_RealK                                         &
      -1.5e+00_RealK*omega*asymmetry
    diff=2.0e+00_RealK*(1.0e+00_RealK-omega)
  ELSE IF (i_2stream == ip_2s_test) THEN
    ! Untested
    summ=1.5e+00_RealK                                         &
      -1.5e+00_RealK*omega*asymmetry
    diff=1.5e+00_RealK*(1.0e+00_RealK-omega)
  ELSE IF (i_2stream == ip_hemi_mean) THEN
    ! Untested
    summ=2.0e+00_RealK                                         &
      *(1.0e+00_RealK-omega*asymmetry)
    diff=2.0e+00_RealK*(1.0e+00_RealK-omega)
  ELSE IF (i_2stream == ip_pifm80) THEN
    summ=2.0e+00_RealK                                         &
      -1.5e+00_RealK*omega*asymmetry                    &
      -0.5e+00_RealK*omega
    diff=2.0e+00_RealK*(1.0e+00_RealK-omega)
  ELSE
    ! Error
    summ=-HUGE(summ)
    diff=-HUGE(diff)
  END IF

END SUBROUTINE two_coeff_basic_scalar

SUBROUTINE two_coeff_cloud(ierr, control                                &
     , n_profile, i_layer_first, i_layer_last                           &
     , i_2stream, n_source_coeff                                        &
     , n_cloud_type, frac_cloud                                         &
     , phase_fnc_cloud, omega_cloud, tau_cloud_dir, tau_cloud           &
     , isolir, sec_0, sph                                               &
     , trans_cloud, reflect_cloud, trans_0_cloud_dir, trans_0_cloud     &
     , source_coeff_cloud                                               &
     , nd_profile, nd_layer, id_ct, nd_max_order                        &
     , nd_source_coeff, nd_cloud_type                                   &
     )

USE realtype_rd, ONLY: RealK
  USE def_spherical_geometry, ONLY: StrSphGeo
  USE rad_pcf, ONLY: &
    ip_2s_test, &
    ip_direct_csr_scaling, &
    ip_direct_noscaling, &
    ip_discrete_ord, &
    ip_eddington, &
    ip_elsasser, &
    ip_hemi_mean, &
    ip_infra_red, &
    ip_pifm80, &
    ip_pifm85, &
    ip_scf_ir_1d, &
    ip_scf_ir_2d, &
    ip_scf_solar_down, &
    ip_scf_solar_up, &
    ip_solar
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE def_control, ONLY: StrCtrl 

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl),      INTENT(IN)    :: control
! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer                                                          &
!       Maximum number of layers
    , id_ct                                                             &
!       Topmost declared potentially cloudy layer
    , nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_source_coeff                                                   &
!       Size allocated for source coefficients
    , nd_cloud_type
!       Maximum number of types of cloud


! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , i_layer_first                                                     &
!       First layer to consider
    , i_layer_last                                                      &
!       Last layer to consider
    , isolir                                                            &
!       Spectral region
    , n_cloud_type                                                      &
!       Number of types of clouds
    , i_2stream                                                         &
!       Two stream scheme
    , n_source_coeff
!       Number of source coefficients

! Optical properties of layer:
  REAL (RealK), INTENT(IN) ::                                           &
      frac_cloud(nd_profile, id_ct: nd_layer, nd_cloud_type)            &
!       Fractions of different types of clouds
    , phase_fnc_cloud(nd_profile, id_ct: nd_layer                       &
        , nd_max_order, nd_cloud_type)                                  &
!       Phase functions
    , omega_cloud(nd_profile, id_ct: nd_layer, nd_cloud_type)           &
!       Albedo of single scattering
    , tau_cloud(nd_profile, id_ct: nd_layer, nd_cloud_type)             &
!       Optical depth
    , tau_cloud_dir(nd_profile, id_ct: nd_layer, nd_cloud_type)
!       Optical depth for direct flux

! Solar beam
  REAL (RealK), INTENT(IN) ::                                           &
      sec_0(nd_profile)
!       Secant of zenith angle

  TYPE(StrSphGeo), INTENT(IN) :: sph
!       Spherical geometry fields


! Coefficients in the two-stream equations:
  REAL (RealK), INTENT(OUT) ::                                          &
      trans_cloud(nd_profile, nd_layer)                                 &
!       Mean diffuse transmission coefficient
    , reflect_cloud(nd_profile, nd_layer)                               &
!       Mean diffuse reflection coefficient
    , trans_0_cloud(nd_profile, nd_layer)                               &
!       Mean direct transmission coefficient
    , trans_0_cloud_dir(nd_profile, nd_layer)                           &
!       Mean direct transmission coefficient with direct tau
    , source_coeff_cloud(nd_profile, nd_layer, nd_source_coeff)
!       Mean source coefficients in two-stream equations


! Local variables.
  LOGICAL :: l_discrete_ord, l_ip_ir, l_ip_solar, l_ir_source_quad
  LOGICAL :: l_no_or_csr, l_spherical_solar

  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , j                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , l
!       Loop variable

  REAL (RealK) :: diff, gamma_down, gamma_up, lambda, reflect
  REAL (RealK) :: source_ir_1d, source_ir_2d, source_solar_down, source_solar_up
  REAL (RealK) :: summ, trans, trans_0, trans_0_dir

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='TWO_COEFF_CLOUD'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  IF ( (control%i_2stream == ip_eddington).OR.                          &
       (control%i_2stream == ip_elsasser).OR.                           &
       (control%i_2stream == ip_pifm85).OR.                             &
       (control%i_2stream == ip_2s_test).OR.                            &
       (control%i_2stream == ip_hemi_mean).OR.                          &
       (control%i_2stream == ip_pifm80) ) THEN

    l_discrete_ord=.FALSE.

  ELSE IF (control%i_2stream == ip_discrete_ord) THEN

    ! Untested
    l_discrete_ord=.TRUE.

  ELSE
    STOP 'ERROR'
  END IF

  l_ip_ir=(control%isolir == ip_infra_red)
  l_ip_solar=(control%isolir == ip_solar)
  l_ir_source_quad=(control%l_ir_source_quad)
  l_no_or_csr=(control%i_direct_tau == ip_direct_noscaling .OR. &
    control%i_direct_tau == ip_direct_csr_scaling)
  IF (control%l_spherical_solar) STOP __LINE__
  l_spherical_solar=.FALSE.

  !$omp target teams distribute parallel do simd collapse(2) &
  !$omp& private(diff, gamma_down, gamma_up, j, k, lambda, reflect) &
  !$omp& private(source_ir_1d, source_ir_2d, source_solar_down, source_solar_up) &
  !$omp& private(summ, trans, trans_0, trans_0_dir)
  DO i=i_layer_first, i_layer_last
    DO l=1, n_profile

      trans_cloud(l, i)=0.0e+00_RealK
      reflect_cloud(l, i)=0.0e+00_RealK
      DO j=1, n_source_coeff
        source_coeff_cloud(l, i, j)=0.0e+00_RealK
      END DO

      IF (isolir == ip_solar) THEN
        trans_0_cloud(l, i)=0.0e+00_RealK
        IF (l_no_or_csr) THEN
          trans_0_cloud_dir(l, i)=0.0e+00_RealK
        END IF
      END IF

      ! Calculate the transmission and reflection coefficients for
      ! each type of cloud and increment the totals, weighting with
      ! the cloud fraction.

      DO k=1, n_cloud_type
        IF (frac_cloud(l, i, k) >  0.0e+00_RealK) THEN

          ! Calculate the basic two-stream coefficients. (The single
          ! scattering albedo has already been perturbed away from 1 in
          ! SINGLE_SCATTERING.)
          CALL two_coeff_basic_scalar( &
            i_2stream &
            , phase_fnc_cloud(l, i, 1, k), omega_cloud(l, i, k) &
            , summ, diff &
            )
          lambda=sqrt(summ*diff)
          ! Calculate the basic coefficients for the solar source terms.
          IF (isolir == ip_solar) THEN
            !   LAMBDA may be perturbed by this routine to avoid
            !   ill-conditioning for the singular zenith angle.
            CALL solar_coefficient_basic_scalar(l_discrete_ord, l_spherical_solar &
              , omega_cloud(l, i, k), phase_fnc_cloud(l, i, 1, k) &
              , sec_0(l) &
              !, sph%common%path_div(l, i) & *** NOT ALLOCATED ***
              , summ, diff, lambda &
              , gamma_up, gamma_down &
              )
          END IF

          ! Determine the transmission and reflection coefficients.
          CALL trans_source_coeff_scalar(l_ip_ir, l_ip_solar &
            , l_ir_source_quad, l_no_or_csr, l_spherical_solar &
            , tau_cloud_dir(l, i, k), tau_cloud(l, i, k), summ, diff &
            , lambda, sec_0(l) &
            !, sph%common%path_div(l, i) & *** NOT ALLOCATED ***
            , gamma_up, gamma_down &
            , trans, reflect, trans_0_dir, trans_0 &
            , source_solar_up, source_solar_down &
            , source_ir_1d, source_ir_2d &
            )

          trans_cloud(l, i)=trans_cloud(l, i)+frac_cloud(l, i, k)*trans
          reflect_cloud(l, i)=reflect_cloud(l, i)+frac_cloud(l, i, k)*reflect
          IF (isolir == ip_solar) THEN
            trans_0_cloud(l, i)=trans_0_cloud(l, i)+frac_cloud(l, i, k)*trans_0
            IF (l_no_or_csr) THEN
              trans_0_cloud_dir(l, i)=trans_0_cloud_dir(l, i) &
                +frac_cloud(l, i, k)*trans_0_dir
            END IF
          END IF
          source_coeff_cloud(l, i, ip_scf_solar_up)= &
            source_coeff_cloud(l, i, ip_scf_solar_up) &
            +frac_cloud(l, i, k)*source_solar_up
          source_coeff_cloud(l, i, ip_scf_solar_down)= &
            source_coeff_cloud(l, i, ip_scf_solar_down) &
            +frac_cloud(l, i, k)*source_solar_down
          source_coeff_cloud(l, i, ip_scf_ir_1d)= &
            source_coeff_cloud(l, i, ip_scf_ir_1d) &
            +frac_cloud(l, i, k)*source_ir_1d
          source_coeff_cloud(l, i, ip_scf_ir_2d)= &
            source_coeff_cloud(l, i, ip_scf_ir_2d) &
            +frac_cloud(l, i, k)*source_ir_2d

        END IF
      END DO

    END DO
  END DO

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE two_coeff_cloud

END MODULE two_coeff_cloud_mod
