! WRF:MODEL_LAYER:PHYSICS
 
    MODULE module_sf_3Dflux

!    USE module_bc, only: set_physical_bc3d
!    USE module_state_description, only: p_m23, p_m13, p_m22, p_m33, p_r23, p_r13, p_r12, p_m12, p_m11
!    USE module_big_step_utilities_em, only: grid_config_rec_type, param_first_scalar, p_qv, p_qi, p_qc
    USE module_model_constants, only: &
       &karman, g, p1000mb, &
       &cp, r_d, r_v, rcp, xlv, xlf, xls, &
       &svp1, svp2, svp3, svpt0, ep_1, ep_2, rvovrd, &
       &cpv, cliq, cice

    IMPLICIT NONE

    REAL, PARAMETER :: xlvcp=xlv/cp, xlscp=(xlv+xlf)/cp, ev=xlv, rd=r_d, &
       &rk=cp/rd, svp11=svp1*1.e3, p608=ep_1, ep_3=1.-ep_2

    REAL, PARAMETER :: tref=300.0     ! reference temperature (K)
    REAL, PARAMETER :: TKmin=253.0    ! for total water conversion, Tripoli and Cotton (1981)
    REAL, PARAMETER :: tv0=p608*tref, tv1=(1.+p608)*tref, gtr=g/tref

    CONTAINS

!=======================================================================
!=======================================================================

    SUBROUTINE Epifanio( u, v, usfc, vsfc,                     &
                             ids, ide, jds, jde, kds, kde,     &
                             ims, ime, jms, jme, kms, kme,     &
                             its, ite, jts, jte, kts, kte      )

! History:     Dec 2016  Module created by Joe Olson and Yi Li.
!
! Purpose:     Full 3D representation of the surface stress.
!
! References:  Epifanio (2007, MWR)
!
! INPUT: U, V, W, vertical eddy viscosity, 
!        HGT (terrain height), height of first model level, 
!        map-scale factors 
!
! OUPUT: 3D stress (m11,m22,m33,m12,m13,m23)
!        1D stress (more mesoscale coupling - enhanced u*) 
!-----------------------------------------------------------------------
! Comments:
! 16-Dec-2016   -  Adding Epifanio's non-local surface velocity calculation
!                  for serial use only. Need alternative for parallel computing.
!
!-----------------------------------------------------------------------
! Begin declarations.

    IMPLICIT NONE

!    TYPE( grid_config_rec_type ), INTENT( IN )  &
!    :: config_flags

    INTEGER, INTENT( IN )  &
    :: ids, ide, jds, jde, kds, kde, &
       ims, ime, jms, jme, kms, kme, &
       its, ite, jts, jte, kts, kte

 !   REAL, INTENT( IN )  &
 !   :: rdx, rdy, cf1, cf2, cf3

 !   REAL, DIMENSION( kms:kme ), INTENT( IN )  &
 !   :: fnm, fnp, dn, dnw, u_base, v_base

!    REAL, DIMENSION( ims:ime , jms:jme ),  INTENT( IN )  &
!    :: msfux, msfuy, msfvx, msfvy, msftx, msfty

!    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )  &
!    ::  u, v, w, zx, zy, rdz, rdzw

!    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( INOUT )  &
!    :: m11, m22, m33, m12, m13, m23 

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN ) :: u
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN ) :: v 

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT( OUT )  &
     :: usfc
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT( OUT )  &
     :: vsfc
   
    INTEGER  &
    :: i, j, k, ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

! End declarations.
!=======================================================================

    i_start = MAX( ids+1, its )
    i_end   = MIN( ite, ide-1 )
    j_start = MAX( jds+1, jts )
    j_end   = MIN( jte, jde-1 )

! Square the map scale factor.
    k=1

    DO j = j_start, j_end
    DO i = i_start, i_end
      usfc(i,j) = u(i,k,j+1) * v(i-1,k,j)
      vsfc(i,j) = u(i-1,k,j+1) * v(i-1,k,j+1)
    END DO
    END DO


END SUBROUTINE EPIFANIO

!-----------------------------------------------------------------------      
!=======================================================================
!=======================================================================

    END MODULE module_sf_3Dflux

!=======================================================================
!=======================================================================
