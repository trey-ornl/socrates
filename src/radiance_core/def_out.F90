! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure for output fields.
!
! Description:
!   This module contains the declaration of the structure
!   used to store output fields for the radiation code.
!
!------------------------------------------------------------------------------
MODULE def_out

USE realtype_rd, ONLY: RealK

IMPLICIT NONE


TYPE StrOut

! Fluxes and radiances calculated
  REAL (RealK), ALLOCATABLE :: flux_direct(:, :, :)
!   Direct flux
  REAL (RealK), ALLOCATABLE :: flux_down(:, :, :)
!   Total downward flux
  REAL (RealK), ALLOCATABLE :: flux_up(:, :, :)
!   Upward flux
  REAL (RealK), ALLOCATABLE :: flux_div(:, :, :)
!   Flux divergence
  REAL (RealK), ALLOCATABLE :: flux_direct_clear(:, :, :)
!   Clear-sky direct flux
  REAL (RealK), ALLOCATABLE :: flux_down_clear(:, :, :)
!   Clear-sky downward flux
  REAL (RealK), ALLOCATABLE :: flux_up_clear(:, :, :)
!   Clear-sky upward flux
  REAL (RealK), ALLOCATABLE :: flux_div_clear(:, :, :)
!   Clear-sky flux divergence
  REAL (RealK), ALLOCATABLE :: flux_direct_div(:, :, :)
!   Direct flux divergence
  REAL (RealK), ALLOCATABLE :: flux_direct_sph(:, :, :)
!   Direct flux for spherical geometry
  REAL (RealK), ALLOCATABLE :: flux_direct_clear_div(:, :, :)
!   Clear-sky direct flux divergence
  REAL (RealK), ALLOCATABLE :: flux_direct_clear_sph(:, :, :)
!   Clear-sky direct flux for spherical geometry
  REAL (RealK), ALLOCATABLE :: radiance(:, :, :, :)
!   Radiances
  REAL (RealK), ALLOCATABLE :: photolysis(:, :, :)
!   Rates of photolysis
  REAL (RealK), ALLOCATABLE :: solar_tail_flux(:)
!   Solar tail flux considered in LW region
  REAL (RealK), ALLOCATABLE :: contrib_funci(:, :, :)
!   Contribution function (intensity)
  REAL (RealK), ALLOCATABLE :: contrib_funcf(:, :, :)
!   Contribution function (flux)

! Diagnostic fluxes
  REAL (RealK), ALLOCATABLE :: flux_up_tile(:, :, :)
!   Upward fluxes at tiled surface points
  REAL (RealK), ALLOCATABLE :: flux_up_blue_tile(:, :, :)
!   Upward blue fluxes at tiled surface points
  REAL (RealK), ALLOCATABLE :: flux_direct_blue_surf(:)
!   Direct blue flux at the surface
  REAL (RealK), ALLOCATABLE :: flux_down_blue_surf(:)
!   Total downward blue flux at the surface
  REAL (RealK), ALLOCATABLE :: flux_up_blue_surf(:)
!   Upward blue flux at the surface
  REAL (RealK), ALLOCATABLE :: flux_direct_band(:, :, :)
!   Direct flux per band
  REAL (RealK), ALLOCATABLE :: flux_direct_div_band(:, :, :)
!   Direct flux divergence per band
  REAL (RealK), ALLOCATABLE :: flux_direct_sph_band(:, :, :)
!   Direct flux for spherical geometry per band
  REAL (RealK), ALLOCATABLE :: flux_down_band(:, :, :)
!   Total (or diffuse for spherical geometry) downward flux per band
  REAL (RealK), ALLOCATABLE :: flux_up_band(:, :, :)
!   Upward flux per band
  REAL (RealK), ALLOCATABLE :: flux_div_band(:, :, :)
!   Flux divergence per band
  REAL (RealK), ALLOCATABLE :: flux_direct_clear_band(:, :, :)
!   Clear-sky direct flux per band
  REAL (RealK), ALLOCATABLE :: flux_direct_clear_div_band(:, :, :)
!   Clear-sky direct flux divergence per band
  REAL (RealK), ALLOCATABLE :: flux_direct_clear_sph_band(:, :, :)
!   Clear-sky direct flux for spherical geometry per band
  REAL (RealK), ALLOCATABLE :: flux_down_clear_band(:, :, :)
!   Clear-sky downward flux per band
  REAL (RealK), ALLOCATABLE :: flux_up_clear_band(:, :, :)
!   Clear-sky upward flux per band
  REAL (RealK), ALLOCATABLE :: flux_div_clear_band(:, :, :)
!   Clear-sky flux divergence per band
  REAL (RealK), ALLOCATABLE :: contrib_funci_band(:, :, :)
!   Contribution function per band (intensity)
  REAL (RealK), ALLOCATABLE :: contrib_funcf_band(:, :, :)
!   Contribution function per band (flux)

! Photolysis diagnostics
  REAL (RealK), ALLOCATABLE :: actinic_flux(:, :, :)
!   Actinic flux
  REAL (RealK), ALLOCATABLE :: actinic_flux_clear(:, :, :)
!   Clear-sky actinic flux
  REAL (RealK), ALLOCATABLE :: actinic_flux_band(:, :, :)
!   Actinic flux per band
  REAL (RealK), ALLOCATABLE :: actinic_flux_clear_band(:, :, :)
!   Clear-sky actinic flux per band
  REAL (RealK), ALLOCATABLE :: photolysis_rate(:, :, :, :)
!   Photolysis rate for each reaction pathway
  REAL (RealK), ALLOCATABLE :: photolysis_rate_clear(:, :, :, :)
!   Clear-sky photolysis rate for each reaction pathway
  REAL (RealK), ALLOCATABLE :: photolysis_div(:, :, :, :)
!   Flux divergence for photolysis for each reaction pathway
  REAL (RealK), ALLOCATABLE :: photolysis_div_clear(:, :, :, :)
!   Clear-sky flux divergence for photolysis for each reaction pathway
 
! Cloud diagnostics
  REAL (RealK), ALLOCATABLE :: tot_cloud_cover(:)
!   Total cloud cover
  REAL (RealK), ALLOCATABLE :: cloud_absorptivity(:, :)
!   Absorptivity of cloud weighted by cloud fraction
!   and upward clear-sky infra-red flux.
  REAL (RealK), ALLOCATABLE :: cloud_weight_absorptivity(:, :)
!   Weights to be applied to absorptivies.
  REAL (RealK), ALLOCATABLE :: ls_cloud_absorptivity(:, :)
!   Absorptivity of layer cloud weighted by cloud fraction
!   and upward clear-sky infra-red flux.
  REAL (RealK), ALLOCATABLE :: ls_cloud_weight_absorptivity(:, :)
!   Weights to be applied to layer cloud absorptivies.
  REAL (RealK), ALLOCATABLE :: cnv_cloud_absorptivity(:, :)
!   Absorptivity of convective cloud weighted by cloud fraction
!   and upward clear-sky infra-red flux.
  REAL (RealK), ALLOCATABLE :: cnv_cloud_weight_absorptivity(:, :)
!   Weights to be applied to convective cloud absorptivies.
  REAL (RealK), ALLOCATABLE :: cloud_extinction(:, :)
!   Extinction of cloud weighted by cloud fraction
!   and downward clear-sky solar flux.
  REAL (RealK), ALLOCATABLE :: cloud_weight_extinction(:, :)
!   Weights to be applied to extinctions.
  REAL (RealK), ALLOCATABLE :: ls_cloud_extinction(:, :)
!   Extinction of layer cloud weighted by cloud fraction
!   and downward clear-sky solar flux.
  REAL (RealK), ALLOCATABLE :: ls_cloud_weight_extinction(:, :)
!   Weights to be applied to layer cloud extinctions.
  REAL (RealK), ALLOCATABLE :: cnv_cloud_extinction(:, :)
!   Extinction of convective cloud weighted by cloud fraction
!   and downward clear-sky solar flux.
  REAL (RealK), ALLOCATABLE :: cnv_cloud_weight_extinction(:, :)
!   Weights to be applied to convective cloud extinctions.

! Aerosol diagnostics
  REAL (RealK), ALLOCATABLE :: aerosol_absorption_band(:, :, :)
!   Total aerosol absorption per band
  REAL (RealK), ALLOCATABLE :: aerosol_scattering_band(:, :, :)
!   Total aerosol scattering per band
  REAL (RealK), ALLOCATABLE :: aerosol_asymmetry_band(:, :, :)
!   Total aerosol asymmetry (weighted by scattering) per band

! Spherical geometry diagnostics
  REAL (RealK), ALLOCATABLE :: spherical_path(:, :, :)
!   Path length for direct beam through spherical layers

END TYPE StrOut


CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE allocate_out(radout, control, dimen, sp)

USE def_control,  ONLY: StrCtrl
USE def_dimen,    ONLY: StrDim
USE def_spectrum, ONLY: StrSpecData

IMPLICIT NONE

TYPE (StrOut),      INTENT(INOUT) :: radout
TYPE (StrCtrl),     INTENT(IN)    :: control
TYPE (StrDim),      INTENT(IN)    :: dimen
TYPE (StrSpecData), INTENT(IN)    :: sp

IF (.NOT. ALLOCATED(radout%flux_direct)) THEN
  ALLOCATE(radout%flux_direct                  ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))
  !$omp target enter data map(alloc: radout%flux_direct)
END IF

IF (.NOT. ALLOCATED(radout%flux_down)) THEN
  ALLOCATE(radout%flux_down                    ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))
  !$omp target enter data map(alloc: radout%flux_down)
END IF

IF (.NOT. ALLOCATED(radout%flux_up)) THEN
  ALLOCATE(radout%flux_up                      ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))
  !$omp target enter data map(alloc: radout%flux_up)
END IF

IF (control%l_flux_div) THEN
  IF (.NOT. ALLOCATED(radout%flux_div)) THEN
    ALLOCATE(radout%flux_div                   ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 dimen%nd_channel            ))
  END IF
END IF

IF (.NOT. ALLOCATED(radout%flux_direct_clear)) THEN
  ALLOCATE(radout%flux_direct_clear            ( dimen%nd_2sg_profile,         &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))
  !$omp target enter data map(alloc: radout%flux_direct_clear)
END IF

IF (.NOT. ALLOCATED(radout%flux_down_clear)) THEN
  ALLOCATE(radout%flux_down_clear              ( dimen%nd_2sg_profile,         &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))
  !$omp target enter data map(alloc: radout%flux_down_clear)
END IF

IF (.NOT. ALLOCATED(radout%flux_up_clear)) THEN
  ALLOCATE(radout%flux_up_clear                ( dimen%nd_2sg_profile,         &
                                                 0: dimen%nd_layer,            &
                                                 dimen%nd_channel            ))
  !$omp target enter data map(alloc: radout%flux_up_clear)
END IF

IF (control%l_flux_div) THEN
  IF (.NOT. ALLOCATED(radout%flux_div_clear)) THEN
    ALLOCATE(radout%flux_div_clear             ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 dimen%nd_channel            ))
  END IF
END IF

IF (.NOT. ALLOCATED(radout%flux_direct_div)) THEN
  ALLOCATE(radout%flux_direct_div              ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 dimen%nd_channel            ))
END IF

IF (.NOT. ALLOCATED(radout%flux_direct_sph)) THEN
  ALLOCATE(radout%flux_direct_sph              ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer+1,          &
                                                 dimen%nd_channel            ))
END IF

IF (.NOT. ALLOCATED(radout%flux_direct_clear_div)) THEN
  ALLOCATE(radout%flux_direct_clear_div        ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 dimen%nd_channel            ))
END IF

IF (.NOT. ALLOCATED(radout%flux_direct_clear_sph)) THEN
  ALLOCATE(radout%flux_direct_clear_sph        ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer+1,          &
                                                 dimen%nd_channel            ))
END IF

IF (.NOT. ALLOCATED(radout%radiance)) THEN
  ALLOCATE(radout%radiance                     ( dimen%nd_radiance_profile,    &
                                                 dimen%nd_viewing_level,       &
                                                 dimen%nd_direction,           &
                                                 dimen%nd_channel            ))
END IF

IF (.NOT. ALLOCATED(radout%photolysis)) THEN
  ALLOCATE(radout%photolysis                   ( dimen%nd_j_profile,           &
                                                 dimen%nd_viewing_level,       &
                                                 dimen%nd_channel            ))
END IF

IF (.NOT. ALLOCATED(radout%solar_tail_flux)) THEN
  ALLOCATE(radout%solar_tail_flux              ( dimen%nd_profile            ))
END IF


IF (.NOT. ALLOCATED(radout%flux_up_tile)) THEN
  ALLOCATE(radout%flux_up_tile                 ( dimen%nd_point_tile,          &
                                                 dimen%nd_tile,                &
                                                 dimen%nd_channel            ))
END IF

IF (.NOT. ALLOCATED(radout%flux_up_blue_tile)) THEN
  ALLOCATE(radout%flux_up_blue_tile            ( dimen%nd_point_tile,          &
                                                 dimen%nd_tile,                &
                                                 dimen%nd_channel            ))
END IF

IF (.NOT. ALLOCATED(radout%flux_direct_blue_surf)) THEN
  ALLOCATE(radout%flux_direct_blue_surf        ( dimen%nd_flux_profile       ))
END IF

IF (.NOT. ALLOCATED(radout%flux_down_blue_surf)) THEN
  ALLOCATE(radout%flux_down_blue_surf          ( dimen%nd_flux_profile       ))
END IF

IF (.NOT. ALLOCATED(radout%flux_up_blue_surf)) THEN
  ALLOCATE(radout%flux_up_blue_surf            ( dimen%nd_flux_profile       ))
END IF


IF (control%l_flux_direct_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_direct_band)) THEN
    ALLOCATE(radout%flux_direct_band           ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_flux_direct_div_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_direct_div_band)) THEN
    ALLOCATE(radout%flux_direct_div_band       ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_flux_direct_sph_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_direct_sph_band)) THEN
    ALLOCATE(radout%flux_direct_sph_band       ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer+1,          &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_flux_down_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_down_band)) THEN
    ALLOCATE(radout%flux_down_band             ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_flux_up_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_up_band)) THEN
    ALLOCATE(radout%flux_up_band               ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer,            &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_flux_div_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_div_band)) THEN
    ALLOCATE(radout%flux_div_band              ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_flux_direct_clear_band .OR.                                      &
     (.NOT.control%l_spherical_solar .AND.                                     &
       ( control%l_cloud_extinction .OR.                                       &
         control%l_ls_cloud_extinction .OR.                                    &
         control%l_cnv_cloud_extinction ) ) ) THEN
  IF (.NOT. ALLOCATED(radout%flux_direct_clear_band)) THEN
    ALLOCATE(radout%flux_direct_clear_band     ( dimen%nd_2sg_profile,         &
                                                 0: dimen%nd_layer,            &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_flux_direct_clear_div_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_direct_clear_div_band)) THEN
    ALLOCATE(radout%flux_direct_clear_div_band ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_flux_direct_clear_sph_band .OR.                                  &
     (control%l_spherical_solar .AND.                                          &
       ( control%l_cloud_extinction .OR.                                       &
         control%l_ls_cloud_extinction .OR.                                    &
         control%l_cnv_cloud_extinction ) ) ) THEN
  IF (.NOT. ALLOCATED(radout%flux_direct_clear_sph_band)) THEN
    ALLOCATE(radout%flux_direct_clear_sph_band ( dimen%nd_flux_profile,        &
                                                 0: dimen%nd_layer+1,          &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_flux_down_clear_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_down_clear_band)) THEN
    ALLOCATE(radout%flux_down_clear_band       ( dimen%nd_2sg_profile,         &
                                                 0: dimen%nd_layer,            &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_flux_up_clear_band .OR.                                          &
    control%l_cloud_absorptivity .OR.                                          &
    control%l_ls_cloud_absorptivity .OR.                                       &
    control%l_cnv_cloud_absorptivity) THEN
  IF (.NOT. ALLOCATED(radout%flux_up_clear_band)) THEN
    ALLOCATE(radout%flux_up_clear_band         ( dimen%nd_2sg_profile,         &
                                                 0: dimen%nd_layer,            &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_flux_div_clear_band) THEN
  IF (.NOT. ALLOCATED(radout%flux_div_clear_band)) THEN
    ALLOCATE(radout%flux_div_clear_band        ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_actinic_flux) THEN
  IF (.NOT. ALLOCATED(radout%actinic_flux)) THEN
    ALLOCATE(radout%actinic_flux               ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 dimen%nd_channel            ))
    !$omp target enter data map(alloc: radout%actinic_flux)
  END IF
END IF

IF (control%l_actinic_flux_clear) THEN
  IF (.NOT. ALLOCATED(radout%actinic_flux_clear)) THEN
    ALLOCATE(radout%actinic_flux_clear         ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 dimen%nd_channel            ))
  END IF
END IF

IF (control%l_actinic_flux_band) THEN
  IF (.NOT. ALLOCATED(radout%actinic_flux_band)) THEN
    ALLOCATE(radout%actinic_flux_band          ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_actinic_flux_clear_band) THEN
  IF (.NOT. ALLOCATED(radout%actinic_flux_clear_band)) THEN
    ALLOCATE(radout%actinic_flux_clear_band    ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_photolysis_rate) THEN
  IF (.NOT. ALLOCATED(radout%photolysis_rate)) THEN
    ALLOCATE(radout%photolysis_rate            ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_pathway,            &
                                                 dimen%nd_channel            ))
  END IF
END IF

IF (control%l_photolysis_rate_clear) THEN
  IF (.NOT. ALLOCATED(radout%photolysis_rate_clear)) THEN
    ALLOCATE(radout%photolysis_rate_clear      ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_pathway,            &
                                                 dimen%nd_channel            ))
  END IF
END IF

IF (control%l_photolysis_div) THEN
  IF (.NOT. ALLOCATED(radout%photolysis_div)) THEN
    ALLOCATE(radout%photolysis_div             ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_pathway,            &
                                                 dimen%nd_channel            ))
  END IF
END IF

IF (control%l_photolysis_div_clear) THEN
  IF (.NOT. ALLOCATED(radout%photolysis_div_clear)) THEN
    ALLOCATE(radout%photolysis_div_clear       ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_pathway,            &
                                                 dimen%nd_channel            ))
  END IF
END IF

IF (.NOT. ALLOCATED(radout%tot_cloud_cover)) THEN
  ALLOCATE(radout%tot_cloud_cover              ( dimen%nd_profile            ))
  !$omp target enter data map(alloc: radout%tot_cloud_cover)
END IF

IF (control%l_cloud_absorptivity) THEN
  IF (.NOT. ALLOCATED(radout%cloud_absorptivity)) THEN
    ALLOCATE(radout%cloud_absorptivity         ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  END IF
  IF (.NOT. ALLOCATED(radout%cloud_weight_absorptivity)) THEN
    ALLOCATE(radout%cloud_weight_absorptivity  ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  END IF
END IF

IF (control%l_ls_cloud_absorptivity) THEN
  IF (.NOT. ALLOCATED(radout%ls_cloud_absorptivity)) THEN
    ALLOCATE(radout%ls_cloud_absorptivity      ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  END IF
  IF (.NOT. ALLOCATED(radout%ls_cloud_weight_absorptivity)) THEN
    ALLOCATE(radout%ls_cloud_weight_absorptivity ( dimen%nd_profile,           &
                                                   dimen%nd_layer            ))
  END IF
END IF

IF (control%l_cnv_cloud_absorptivity) THEN
  IF (.NOT. ALLOCATED(radout%cnv_cloud_absorptivity)) THEN
    ALLOCATE(radout%cnv_cloud_absorptivity     ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  END IF
  IF (.NOT. ALLOCATED(radout%cnv_cloud_weight_absorptivity)) THEN
    ALLOCATE(radout%cnv_cloud_weight_absorptivity ( dimen%nd_profile,          &
                                                    dimen%nd_layer           ))
  END IF
END IF

IF (control%l_cloud_extinction) THEN
  IF (.NOT. ALLOCATED(radout%cloud_extinction)) THEN
    ALLOCATE(radout%cloud_extinction           ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  END IF
  IF (.NOT. ALLOCATED(radout%cloud_weight_extinction)) THEN
    ALLOCATE(radout%cloud_weight_extinction    ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  END IF
END IF

IF (control%l_ls_cloud_extinction) THEN
  IF (.NOT. ALLOCATED(radout%ls_cloud_extinction)) THEN
    ALLOCATE(radout%ls_cloud_extinction        ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  END IF
  IF (.NOT. ALLOCATED(radout%ls_cloud_weight_extinction)) THEN
    ALLOCATE(radout%ls_cloud_weight_extinction ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  END IF
END IF

IF (control%l_cnv_cloud_extinction) THEN
  IF (.NOT. ALLOCATED(radout%cnv_cloud_extinction)) THEN
    ALLOCATE(radout%cnv_cloud_extinction       ( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  END IF
  IF (.NOT. ALLOCATED(radout%cnv_cloud_weight_extinction)) THEN
    ALLOCATE(radout%cnv_cloud_weight_extinction( dimen%nd_profile,             &
                                                 dimen%nd_layer              ))
  END IF
END IF

IF (control%l_aerosol_absorption_band) THEN
  IF (.NOT. ALLOCATED(radout%aerosol_absorption_band)) THEN
    ALLOCATE(radout%aerosol_absorption_band    ( dimen%nd_profile,             &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_aerosol_scattering_band) THEN
  IF (.NOT. ALLOCATED(radout%aerosol_scattering_band)) THEN
    ALLOCATE(radout%aerosol_scattering_band    ( dimen%nd_profile,             &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_aerosol_asymmetry_band) THEN
  IF (.NOT. ALLOCATED(radout%aerosol_asymmetry_band)) THEN
    ALLOCATE(radout%aerosol_asymmetry_band     ( dimen%nd_profile,             &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
  END IF
END IF

IF (control%l_spherical_path_diag) THEN
  IF (.NOT. ALLOCATED(radout%spherical_path)) THEN
    ALLOCATE(radout%spherical_path             ( dimen%nd_profile,             &
                                                 dimen%nd_layer,               &
                                                 0:dimen%nd_layer+1          ))
  END IF
END IF

IF (control%l_contrib_func) THEN
  IF (.NOT. ALLOCATED(radout%contrib_funci)) THEN
    ALLOCATE(radout%contrib_funci              ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 dimen%nd_channel            ))
  END IF
  IF (.NOT. ALLOCATED(radout%contrib_funcf)) THEN
    ALLOCATE(radout%contrib_funcf              ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 dimen%nd_channel            ))
  END IF
END IF

IF (control%l_contrib_func_band) THEN
  IF (.NOT. ALLOCATED(radout%contrib_funci_band)) THEN
    ALLOCATE(radout%contrib_funci_band         ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
  END IF
  IF (.NOT. ALLOCATED(radout%contrib_funcf_band)) THEN
    ALLOCATE(radout%contrib_funcf_band         ( dimen%nd_flux_profile,        &
                                                 dimen%nd_layer,               &
                                                 sp%dim%nd_band              ))
  END IF
END IF

END SUBROUTINE allocate_out
!------------------------------------------------------------------------------
SUBROUTINE deallocate_out(radout)

IMPLICIT NONE

TYPE (StrOut), INTENT(INOUT) :: radout

IF (ALLOCATED(radout%contrib_funcf_band)) THEN
    DEALLOCATE(radout%contrib_funcf_band)
END IF
IF (ALLOCATED(radout%contrib_funci_band)) THEN
    DEALLOCATE(radout%contrib_funci_band)
END IF
IF (ALLOCATED(radout%contrib_funcf)) THEN
    DEALLOCATE(radout%contrib_funcf)
END IF
IF (ALLOCATED(radout%contrib_funci)) THEN
    DEALLOCATE(radout%contrib_funci)
END IF
IF (ALLOCATED(radout%spherical_path)) THEN
    DEALLOCATE(radout%spherical_path)
END IF
IF (ALLOCATED(radout%aerosol_asymmetry_band)) THEN
    DEALLOCATE(radout%aerosol_asymmetry_band)
END IF
IF (ALLOCATED(radout%aerosol_scattering_band)) THEN
    DEALLOCATE(radout%aerosol_scattering_band)
END IF
IF (ALLOCATED(radout%aerosol_absorption_band)) THEN
    DEALLOCATE(radout%aerosol_absorption_band)
END IF
IF (ALLOCATED(radout%cnv_cloud_weight_extinction)) THEN
    DEALLOCATE(radout%cnv_cloud_weight_extinction)
END IF
IF (ALLOCATED(radout%cnv_cloud_extinction)) THEN
    DEALLOCATE(radout%cnv_cloud_extinction)
END IF
IF (ALLOCATED(radout%ls_cloud_weight_extinction)) THEN
    DEALLOCATE(radout%ls_cloud_weight_extinction)
END IF
IF (ALLOCATED(radout%ls_cloud_extinction)) THEN
    DEALLOCATE(radout%ls_cloud_extinction)
END IF
IF (ALLOCATED(radout%cloud_weight_extinction)) THEN
    DEALLOCATE(radout%cloud_weight_extinction)
END IF
IF (ALLOCATED(radout%cloud_extinction)) THEN
    DEALLOCATE(radout%cloud_extinction)
END IF
IF (ALLOCATED(radout%cnv_cloud_weight_absorptivity)) THEN
    DEALLOCATE(radout%cnv_cloud_weight_absorptivity)
END IF
IF (ALLOCATED(radout%cnv_cloud_absorptivity)) THEN
    DEALLOCATE(radout%cnv_cloud_absorptivity)
END IF
IF (ALLOCATED(radout%ls_cloud_weight_absorptivity)) THEN
    DEALLOCATE(radout%ls_cloud_weight_absorptivity)
END IF
IF (ALLOCATED(radout%ls_cloud_absorptivity)) THEN
    DEALLOCATE(radout%ls_cloud_absorptivity)
END IF
IF (ALLOCATED(radout%cloud_weight_absorptivity)) THEN
    DEALLOCATE(radout%cloud_weight_absorptivity)
END IF
IF (ALLOCATED(radout%cloud_absorptivity)) THEN
    DEALLOCATE(radout%cloud_absorptivity)
END IF
IF (ALLOCATED(radout%tot_cloud_cover)) THEN
  !$omp target exit data map(delete: radout%tot_cloud_cover)
  DEALLOCATE(radout%tot_cloud_cover)
END IF
IF (ALLOCATED(radout%photolysis_div_clear)) THEN
    DEALLOCATE(radout%photolysis_div_clear)
END IF
IF (ALLOCATED(radout%photolysis_div)) THEN
    DEALLOCATE(radout%photolysis_div)
END IF
IF (ALLOCATED(radout%photolysis_rate_clear)) THEN
    DEALLOCATE(radout%photolysis_rate_clear)
END IF
IF (ALLOCATED(radout%photolysis_rate)) THEN
    DEALLOCATE(radout%photolysis_rate)
END IF
IF (ALLOCATED(radout%actinic_flux_clear_band)) THEN
    DEALLOCATE(radout%actinic_flux_clear_band)
END IF
IF (ALLOCATED(radout%actinic_flux_band)) THEN
    DEALLOCATE(radout%actinic_flux_band)
END IF
IF (ALLOCATED(radout%actinic_flux_clear)) THEN
    DEALLOCATE(radout%actinic_flux_clear)
END IF
IF (ALLOCATED(radout%actinic_flux)) THEN
    !$omp target exit data map(delete: radout%actinic_flux)
    DEALLOCATE(radout%actinic_flux)
END IF
IF (ALLOCATED(radout%flux_div_clear_band)) THEN
    DEALLOCATE(radout%flux_div_clear_band)
END IF
IF (ALLOCATED(radout%flux_up_clear_band)) THEN
    DEALLOCATE(radout%flux_up_clear_band)
END IF
IF (ALLOCATED(radout%flux_down_clear_band)) THEN
    DEALLOCATE(radout%flux_down_clear_band)
END IF
IF (ALLOCATED(radout%flux_direct_clear_sph_band)) THEN
    DEALLOCATE(radout%flux_direct_clear_sph_band)
END IF
IF (ALLOCATED(radout%flux_direct_clear_div_band)) THEN
    DEALLOCATE(radout%flux_direct_clear_div_band)
END IF
IF (ALLOCATED(radout%flux_direct_clear_band)) THEN
    DEALLOCATE(radout%flux_direct_clear_band)
END IF
IF (ALLOCATED(radout%flux_div_band)) THEN
    DEALLOCATE(radout%flux_div_band)
END IF
IF (ALLOCATED(radout%flux_up_band)) THEN
    DEALLOCATE(radout%flux_up_band)
END IF
IF (ALLOCATED(radout%flux_down_band)) THEN
    DEALLOCATE(radout%flux_down_band)
END IF
IF (ALLOCATED(radout%flux_direct_sph_band)) THEN
    DEALLOCATE(radout%flux_direct_sph_band)
END IF
IF (ALLOCATED(radout%flux_direct_div_band)) THEN
    DEALLOCATE(radout%flux_direct_div_band)
END IF
IF (ALLOCATED(radout%flux_direct_band)) THEN
    DEALLOCATE(radout%flux_direct_band)
END IF
IF (ALLOCATED(radout%flux_up_blue_surf)) THEN
    DEALLOCATE(radout%flux_up_blue_surf)
END IF
IF (ALLOCATED(radout%flux_down_blue_surf)) THEN
    DEALLOCATE(radout%flux_down_blue_surf)
END IF
IF (ALLOCATED(radout%flux_direct_blue_surf)) THEN
    DEALLOCATE(radout%flux_direct_blue_surf)
END IF
IF (ALLOCATED(radout%flux_up_blue_tile)) THEN
    DEALLOCATE(radout%flux_up_blue_tile)
END IF
IF (ALLOCATED(radout%flux_up_tile)) THEN
    DEALLOCATE(radout%flux_up_tile)
END IF
IF (ALLOCATED(radout%solar_tail_flux)) THEN
    DEALLOCATE(radout%solar_tail_flux)
END IF
IF (ALLOCATED(radout%photolysis)) THEN
    DEALLOCATE(radout%photolysis)
END IF
IF (ALLOCATED(radout%radiance)) THEN
    DEALLOCATE(radout%radiance)
END IF
IF (ALLOCATED(radout%flux_direct_clear_sph)) THEN
    DEALLOCATE(radout%flux_direct_clear_sph)
END IF
IF (ALLOCATED(radout%flux_direct_clear_div)) THEN
    DEALLOCATE(radout%flux_direct_clear_div)
END IF
IF (ALLOCATED(radout%flux_direct_sph)) THEN
    DEALLOCATE(radout%flux_direct_sph)
END IF
IF (ALLOCATED(radout%flux_direct_div)) THEN
    DEALLOCATE(radout%flux_direct_div)
END IF
IF (ALLOCATED(radout%flux_div_clear)) THEN
    DEALLOCATE(radout%flux_div_clear)
END IF
IF (ALLOCATED(radout%flux_up_clear)) THEN
    !$omp target exit data map(delete: radout%flux_up_clear)
    DEALLOCATE(radout%flux_up_clear)
END IF
IF (ALLOCATED(radout%flux_down_clear)) THEN
    !$omp target exit data map(delete: radout%flux_down_clear)
    DEALLOCATE(radout%flux_down_clear)
END IF
IF (ALLOCATED(radout%flux_direct_clear)) THEN
    !$omp target exit data map(delete: radout%flux_direct_clear)
    DEALLOCATE(radout%flux_direct_clear)
END IF
IF (ALLOCATED(radout%flux_div)) THEN
    DEALLOCATE(radout%flux_div)
END IF
IF (ALLOCATED(radout%flux_up)) THEN
    !$omp target exit data map(delete: radout%flux_up)
    DEALLOCATE(radout%flux_up)
END IF
IF (ALLOCATED(radout%flux_down)) THEN
    !$omp target exit data map(delete: radout%flux_down)
    DEALLOCATE(radout%flux_down)
END IF
IF (ALLOCATED(radout%flux_direct)) THEN
    !$omp target exit data map(delete: radout%flux_direct)
    DEALLOCATE(radout%flux_direct)
END IF

END SUBROUTINE deallocate_out
!------------------------------------------------------------------------------

END MODULE def_out
