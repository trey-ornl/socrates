! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to apply a path-length scaling to the continuum.
!
! Method:
!   The scaling function is calculated. This is multpiled by a
!   suitable "amount" of continuum incorporating a broadening
!   density.
!
!- ---------------------------------------------------------------------
MODULE rescale_continuum_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'RESCALE_CONTINUUM_MOD'
CONTAINS
SUBROUTINE rescale_continuum(control, n_profile, n_layer, i_continuum   &
     , p, t, density, water_frac                                        &
     , amount_continuum                                                 &
     , i_fnc                                                            &
     , p_reference, t_reference, scale_parameter                        &
     , nd_profile, nd_layer                                             &
     , nd_scale_variable                                                &
     )


  USE realtype_rd, ONLY: RealK
  USE rad_pcf, ONLy: ip_frn_continuum, ip_self_continuum,               &
                     ip_scale_power_law, ip_scale_power_quad
  USE rad_ccf, ONLY: mol_weight_air, repsilon
  USE def_control, ONLY: StrCtrl
  USE vectlib_mod, ONLY : rtor_v
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Control options:
  TYPE(StrCtrl), INTENT(IN)  :: control

! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_scale_variable
!       Size allocated for scaling variables


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , i_continuum                                                       &
!       Continuum type
    , i_fnc
!       Scaling function
  REAL (RealK), INTENT(IN) ::                                           &
      water_frac(nd_profile, nd_layer)                                  &
!       Mass fraction of water
    , p(nd_profile, nd_layer)                                           &
!       Pressure
    , t(nd_profile, nd_layer)                                           &
!       Temperature
    , density(nd_profile, nd_layer)                                     &
!       Overall density
    , p_reference                                                       &
!       Reference pressure
    , t_reference                                                       &
!       Reference pressure
    , scale_parameter(nd_scale_variable)
!       Scaling paramters
  REAL (RealK), INTENT(OUT) ::                                          &
      amount_continuum(nd_profile, nd_layer)
!       Amount of continuum

! Local variables.
  LOGICAL :: l_mixing_ratio
  INTEGER :: l, i
!       Loop variables
  REAL (RealK) :: molar_density_water
!       Molar density of water vapour
  REAL (RealK) :: pwk_in
  REAL (RealK) :: pwk
  REAL (RealK) :: twk_in
  REAL (RealK) :: twk
  REAL (RealK) :: sp1
  REAL (RealK) :: sp2

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  CHARACTER(LEN=*), PARAMETER :: RoutineName='RESCALE_CONTINUUM'


  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  l_mixing_ratio = control%l_mixing_ratio

  !$omp target teams distribute parallel do simd collapse(2) &
  !$omp& private(molar_density_water, pwk, pwk_in, sp1, sp2, twk, twk_in) 
  DO i=1, n_layer
    DO l=1, n_profile

      sp1=scale_parameter(1)
      sp2=scale_parameter(2)
      pwk_in=p(l, i)/p_reference
      pwk = pwk_in**sp1

      IF (i_fnc == ip_scale_power_law) THEN
        twk_in=t(l, i)/t_reference
        twk = twk_in**sp2
        amount_continuum(l, i)=pwk*twk
      ELSE IF(i_fnc == ip_scale_power_quad) THEN
        amount_continuum(l, i)                                          &
          =pwk                                                    &
          *(1.0e+00+scale_parameter(2)*(t(l, i)                        &
          /t_reference-1.0e+00)                                        &
          +scale_parameter(3)*(t(l, i)                                 &
          /t_reference-1.0e+00)**2)
      END IF

      IF (i_continuum == ip_self_continuum) THEN
        molar_density_water=density(l, i)                         &
          *water_frac(l, i)/(repsilon*mol_weight_air)
        amount_continuum(l, i)=amount_continuum(l, i)                   &
          *molar_density_water*water_frac(l, i)
      ELSE IF (i_continuum == ip_frn_continuum) THEN
        IF (l_mixing_ratio) THEN
          ! In this case density and mass are for the dry component
          amount_continuum(l, i)=amount_continuum(l, i)                 &
            *density(l, i)*water_frac(l, i)/mol_weight_air
        ELSE
          amount_continuum(l, i)=amount_continuum(l, i)                 &
            *density(l, i)*water_frac(l, i)*(1.0_RealK-water_frac(l, i))&
            /mol_weight_air
        END IF
      END IF
    END DO
  END DO

  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE rescale_continuum
END MODULE rescale_continuum_mod
