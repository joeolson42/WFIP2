! WRF:MODEL_LAYER:PHYSICS
 
    MODULE module_3Dsfcflux

    USE module_bc, only: set_physical_bc3d
    USE module_state_description, only: p_m23, p_m13, p_m22, p_m33, p_r23, p_r13, p_r12, p_m12, p_m11
    USE module_big_step_utilities_em, only: grid_config_rec_type, param_first_scalar, p_qv, p_qi, p_qc
    USE module_model_constants, only: &
       &karman, g, p1000mb, &
       &cp, r_d, r_v, rcp, xlv, xlf, xls, &
       &svp1, svp2, svp3, svpt0, ep_1, ep_2, rvovrd, &
       &cpv, cliq, cice

      real g,pii,alat,rgas,cp,gamma,gcp,rcp,rcpinv,cv,cvr,rcv,cpinv
      parameter(g = 9.806,pii = 3.141592)
      parameter(alat = 2.5e6,rgas = 287.04,cp = 1004.,gamma = 1.4)
      parameter(gcp = g/cp,rcp = rgas/cp,rcpinv = cp/rgas)
      parameter(cv = 717.0,cvr = cv/rgas,rcv = rgas/cv)
      parameter(cpinv = 1./cp)


    IMPLICIT NONE

    REAL, PARAMETER :: xlvcp=xlv/cp, xlscp=(xlv+xlf)/cp, ev=xlv, rd=r_d, &
       &rk=cp/rd, svp11=svp1*1.e3, p608=ep_1, ep_3=1.-ep_2

    REAL, PARAMETER :: tref=300.0     ! reference temperature (K)
    REAL, PARAMETER :: TKmin=253.0    ! for total water conversion, Tripoli and Cotton (1981)
    REAL, PARAMETER :: tv0=p608*tref, tv1=(1.+p608)*tref, gtr=g/tref

    CONTAINS

!=======================================================================
!=======================================================================

    SUBROUTINE Epifanio( config_flags, u, v, w, HGT, dz       &
                                   m11, m22, m33,             &
                                   m12, m13, m23,             &
                                   u_base, v_base, msfux, msfuy,     &
                                   msfvx, msfvy, msftx, msfty,       &
                                   rdx, rdy, dn, dnw, rdz, rdzw,     &
                                   fnm, fnp, cf1, cf2, cf3, zx, zy,  &
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

    TYPE( grid_config_rec_type ), INTENT( IN )  &
    :: config_flags

    INTEGER, INTENT( IN )  &
    :: ids, ide, jds, jde, kds, kde, &
       ims, ime, jms, jme, kms, kme, &
       its, ite, jts, jte, kts, kte

    REAL, INTENT( IN )  &
    :: rdx, rdy, cf1, cf2, cf3

    REAL, DIMENSION( kms:kme ), INTENT( IN )  &
    :: fnm, fnp, dn, dnw, u_base, v_base

    REAL, DIMENSION( ims:ime , jms:jme ),  INTENT( IN )  &
    :: msfux, msfuy, msfvx, msfvy, msftx, msfty

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( IN )  &
    ::  u, v, w, zx, zy, rdz, rdzw

    REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT( INOUT )  &
    :: m11, m22, m33, m12, m13, m23 

!!!!!!!      include 'grid.h'
      integer ::  nsint, nsreal, nx1d, ny1d, nz1d, nxy2d, nxz2d, & 
                  nyz2d, nxyz3d, ivrsn 
      INTEGER, parameter :: nsint  = 18, & ! number of integer constants
                            nsreal = 60, & ! number of real constants
                            nx1d   = 2,  & ! number of 1-d arrays in x
                            ny1d   = 2,  & ! number of 1-d arrays in y
                            nz1d   = 23, & ! number of 1-d arrays in z
                            nxy2d  = 3,  & ! number of 2-d xy arrays
                            nxz2d  = 0,  & ! number of 2-d xz arrays
                            nyz2d  = 0,  & ! number of 2-d yz arrays
                            nxyz3d = 23, & ! number of 3-d xyz arrays
                            ivrsn = 0      ! version of model
                                           ! 0:  dry model  
                                           ! 1:  wet model
                                           ! 2:  ice model
!!!!!!!!      include 'comblk.h'
c Include the scratch file to declare tmp space for needed solver data
c----------------------------------------------------------------------
!      include 'scratch.h'
      integer nmax !,nxmax,nymax,nzmax
!      parameter (nxmax=201,nymax=201,nzmax=161)
!             nxmax = ime, nymax=jme, nzmax=kme 
      parameter (nmax=201) 
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c Logical constants
c----------------------------------------------------------------------
      logical top_radiation,twod,mixing,bouss,hydrostatic,          &
              spongz,spongxy,viscous,outflow_only,outflowx0,outflowxl,&
              outflowy0,outflowyl,nested,symx0,symxl,symy0,symyl,   &
              perx,pery,constk
      common /lflags/ top_radiation,twod,mixing,bouss,hydrostatic,  &
              spongz,spongxy,viscous,outflow_only,outflowx0,        &
              outflowxl,outflowy0,outflowyl,nested,symx0,symxl,     &
              symy0,symyl,perx,pery,constk
c----------------------------------------------------------------------
c Integer constants
c----------------------------------------------------------------------
      integer nsm,bcx,bcy,bcz,opts_spec,type_terrain,idiag,hflux_opt,&
              drag_opt
      common /iconstlist/ nsm,bcx,bcy,bcz,opts_spec,type_terrain,   &
           idiag,hflux_opt,drag_opt
c----------------------------------------------------------------------
c Real constants 
c----------------------------------------------------------------------
      real gx0,gy0,gx1,gy1,dt,dx,dy,dz,dzbot,um,vm,dxt,eps0,eps1,tf0, &
           tf1,fdamph,alowk0,rayd_hgt,rayd_mag,dmpxy_l,dmpxy_a,rprndl,&
           Cm,ahighk,psfc,dts,rdt,rdx,rdy,rdz,dtbydx,dtbydy,dtbydz,   &
           coriol,lmda1,lmda2,lmda3,thets,pis,Nsrf,shxc,shyc,shlx,    &
           shly,shlz,shamp,shfreq,shramp,rdiag,rdiag2,bet,dusrf,cdmax,&
           cdmin,sdxc,sdyc,sdlx,sdly
      common /rconstlist/ gx0,gy0,gx1,gy1,dt,dx,dy,dz,dzbot,um,vm,    &
           dxt,eps0,eps1,tf0,tf1,fdamph,alowk0,rayd_hgt,rayd_mag,     &
           dmpxy_l,dmpxy_a,rprndl,Cm,ahighk,psfc,dts,rdt,rdx,rdy,rdz, &
           dtbydx,dtbydy,dtbydz,coriol,lmda1,lmda2,lmda3,thets,pis,   &
           Nsrf,shxc,shyc,shlx,shly,shlz,shamp,shfreq,shramp,rdiag,   &
           rdiag2,bet,dusrf,cdmax,cdmin,sdxc,sdyc,sdlx,sdly
c----------------------------------------------------------------------
c X1D arrays
c----------------------------------------------------------------------
      real sgx(nmax),ugx(nmax),xbnd(nmax),xrbnd(nmax),xlbnd(nmax)
      common /x1dlist/ sgx,ugx,xbnd,xrbnd,xlbnd
c----------------------------------------------------------------------
c Y1D arrays
c----------------------------------------------------------------------
      real sgy(nmax),vgy(nmax),ybnd(nmax),yrbnd(nmax),ylbnd(nmax)
      common /y1dlist/ sgy,vgy,ybnd,yrbnd,ylbnd
c----------------------------------------------------------------------
c Z1D arrays
c----------------------------------------------------------------------
      real sgz(kme),wgz(kme),mfc(kme),mfe(kme),zbnd(kme),             &
           zlbnd(kme),zrbnd(kme),zlbndw(kme),zrbndw(kme),             &
           alowk(kme),delta(kme),kmh(kme),kmz(kme),                   &
           atri(kme),btri(kme),ctri(kme),atrb(kme),                   &
           btrb(kme),ctrb(kme),gam(kme),aterm(kme),                   &
           rayd(kme),uinit(kme),vinit(kme),winit(kme),                &
           sinit(kme,7),piinit(kme),rhinit(kme),thalf(kme),           &
           cs2(kme),dzc(-2:kme+3),dalfac(-2:kme+3),                   &
           darfac(-2:kme+3),km_min(kme),km_max(kme),                  &
           csmths(kme,5),csmthw(kme,5)
      common /z1dlist/ sgz,wgz,mfc,mfe,zbnd,zlbnd,zrbnd,zlbndw,zrbndw,&
           alowk,delta,kmh,kmz,atri,btri,ctri,atrb,btrb,ctrb,gam,     &
           aterm,rayd,uinit,vinit,winit,sinit,piinit,rhinit,thalf,cs2,&
           dzc,dalfac,darfac,km_min,km_max,csmths,csmthw
c----------------------------------------------------------------------
c Terrain coordinate transformation arrays:
c   1: therm pts    2: u pts     3: v pts
c----------------------------------------------------------------------

      real dqdz(ime,jme,3),dzdq(ime,jme,3)
      real dqdx_xy(ime,jme,3),dqdy_xy(ime,jme,3)
      real dqdx_zt(kme),dqdy_zt(kme),dqdx_zw(kme),dqdy_zw(kme)
      common /metrics/ dqdz,dzdq,dqdx_xy,dqdy_xy,dqdx_zt,dqdy_zt,     &
                       dqdx_zw,dqdy_zw
c----------------------------------------------------------------------
c Arrays for surface stress definitions
c----------------------------------------------------------------------
      real xsdu(ime,jme),ysdu(ime,jme),dsu(ime,jme)
      real xsdv(ime,jme),ysdv(ime,jme),dsv(ime,jme)
      common /stress/ xsdu,ysdu,dsu,xsdv,ysdv,dsv
c----------------------------------------------------------------------
c Sponge layer arrays
c----------------------------------------------------------------------
      real rdmpxy(ime,jme,3)
      common /sponge/ rdmpxy


      common/info/cstep,dxf,dyf,dtf
      integer cstep
      real dxf,dyf,dtf
!----------------------------------------------------------------------
      integer iconst(nsint),nx,ny,nz
      real rconst(nsreal),x1d(nx,nx1d),y1d(ny,ny1d),z1d(nz,nz1d)
      real zbot(nx,ny,3)
!----------------------------------------------------------------------
      real tvinit(kme) 
!      real    pltmp(nmax,nmax)
      real fdamph0,alpha,timefilt,pr
      integer i,j,k,n,km1,iconv,maxit
      real dtl,dxd,dzd,dz1,avg0,avg1,gammai,zdif,z0,zfac,dz2,avg2
      real dztmp,ztop,dz0,x,y,z,xd0,yd0,xdl,ydl,del
      real x_0,y_0,z_0,x_scale,y_scale,z_scale,x_scale2,y_scale2,y_len
      real xsrf,ysrf,zsrf,h0,n1,n2,n3,gamm,tol,res,h_surf

! Local variables.

    INTEGER  &
    :: i, j, k, ktf, ktes1, ktes2, i_start, i_end, j_start, j_end

    REAL  &
    :: tmp, tmpzx, tmpzy, tmpzeta_z, cft1, cft2

    REAL, DIMENSION( its:ite, jts:jte )  &
    :: mm, zzavg, zeta_zd12

    REAL, DIMENSION( its-2:ite+2, kts:kte, jts-2:jte+2 )  &
    :: tmp1, hat, hatavg

! End declarations.
!=======================================================================
! In the following section, calculate 3-d divergence and the first three
! (defor11, defor22, defor33) of six deformation terms.

    ktes1   = kte-1
    ktes2   = kte-2

    cft2    = - 0.5 * dnw(ktes1) / dn(ktes1)
    cft1    = 1.0 - cft2

    ktf     = MIN( kte, kde-1 )

    i_start = its
    i_end   = MIN( ite, ide-1 )
    j_start = jts
    j_end   = MIN( jte, jde-1 )

! Square the map scale factor.

    DO j = j_start, j_end
    DO i = i_start, i_end
      mm(i,j) = msftx(i,j) * msfty(i,j)
    END DO
    END DO

    IF ( config_flags%open_xs .OR. config_flags%specified .OR. &
          config_flags%nested) i_start = MAX( ids+1, its )
    IF ( config_flags%open_ys .OR. config_flags%specified .OR. &
         config_flags%nested) j_start = MAX( jds+1, jts )
    IF ( config_flags%periodic_y ) j_end = MIN( jte, jde )
    IF ( config_flags%periodic_x ) i_start = its
    IF ( config_flags%periodic_x ) i_end = MIN( ite, ide-1 )

!-----------------------------------------------------------------------------
! Coordinate metrics.
!    1: therm pts
!    2: u pts
!    3: v pts
! Start with dqdz, dzdq (q is transformed distance in vertical).
!---------------------------------------------------------------------------- 
      do j=j_start,j_end
      do i=i_start,i_end
        dqdz(i,j,1) = ztop/(ztop-zbot(i,j,1))
        dqdz(i,j,2) = ztop/(ztop-zbot(i,j,2))
        dqdz(i,j,3) = ztop/(ztop-zbot(i,j,3))
        dzdq(i,j,1) = 1./dqdz(i,j,1)
        dzdq(i,j,2) = 1./dqdz(i,j,2)
        dzdq(i,j,3) = 1./dqdz(i,j,3)
      enddo
      enddo 
!----------------------------------------------------------------------------
! x and y dependent parts of dqdx.  Boundary, then interior.  We use 
!   all the terrain pts here (i.e., u, v, and therm).  Note that 
!   dqdx_xy(i,j,3) is normally not used during model integrations. 
!---------------------------------------------------------------------------
      do j=j_start,j_end
        dqdx_xy(1,j,1) = ztop/(ztop-zbot(1,j,1))     &
                     * rdx*(zbot(2,j,2)-zbot(1,j,2))
        dqdx_xy(1,j,2) = ztop/(ztop-zbot(1,j,2))     &
                      * 2.*rdx*(zbot(1,j,1)-zbot(1,j,2))
        dqdx_xy(1,j,3) = ztop/(ztop-zbot(1,j,3))     &
                      * rdx*(zbot(2,j,3)-zbot(1,j,3))
        dqdx_xy(i_end,j,2) = ztop/(ztop-zbot(i_end,j,2))   &
                      * 2.*rdx*(zbot(i_end,j,2)-zbot(i_end-1,j,1))
        dqdx_xy(i_end-1,j,1) = ztop/(ztop-zbot(i_end-1,j,1)) &
                      * rdx*(zbot(i_end,j,2)-zbot(i_end-1,j,2))
        dqdx_xy(i_end-1,j,2) = ztop/(ztop-zbot(i_end-1,j,2)) &
                * (1./6.)*rdx*( 8.*(zbot(i_end-1,j,1)-zbot(i_end-2,j,1)) &
                      - (zbot(i_end,j,2)-zbot(i_end-2,j,2)) )
        dqdx_xy(i_end-1,j,3) = ztop/(ztop-zbot(i_end-1,j,3)) &
                      * rdx*(zbot(i_end-1,j,3)-zbot(i_end-2,j,3))
      enddo

      do j=j_start,j_end
      do i=i_start+1,i_end-2
        dqdx_xy(i,j,1) = ztop/(ztop-zbot(i,j,1))     &
                * (1./6.)*rdx*( 8.*(zbot(i+1,j,2)-zbot(i,j,2)) &
                      - (zbot(i+1,j,1)-zbot(i-1,j,1)) )
        dqdx_xy(i,j,2) = ztop/(ztop-zbot(i,j,2))     &
                * (1./6.)*rdx*( 8.*(zbot(i,j,1)-zbot(i-1,j,1)) &
                      - (zbot(i+1,j,2)-zbot(i-1,j,2)) )
        dqdx_xy(i,j,3) = ztop/(ztop-zbot(i,j,3))     &
                      * 0.5*rdx*(zbot(i+1,j,3)-zbot(i-1,j,3))
      enddo
      enddo 
!-----------------------------------------------------------------------------
! x and y dependent parts of dqdy.  Boundary, then interior.  We use 
!   all the terrain pts here (i.e., u, v, and therm).  Note that
!   dqdy_xy(i,j,2) is normally not used during model integrations.
!----------------------------------------------------------------------------
      do i=i_start,i_end
        dqdy_xy(i,1,1) = ztop/(ztop-zbot(i,1,1)) &
                      * rdy*(zbot(i,2,3)-zbot(i,1,3))
        dqdy_xy(i,1,2) = ztop/(ztop-zbot(i,1,2)) &
                      * rdy*(zbot(i,2,2)-zbot(i,1,2))
        dqdy_xy(i,1,3) = ztop/(ztop-zbot(i,1,3)) &
                      * 2.*rdy*(zbot(i,1,1)-zbot(i,1,3))
        dqdy_xy(i,j_end,3) = ztop/(ztop-zbot(i,j_end,3)) &
                      * 2.*rdy*(zbot(i,j_end,3)-zbot(i,j_end-1,1))
        dqdy_xy(i,j_end-1,1) = ztop/(ztop-zbot(i,j_end-1,1)) &
                      * rdy*(zbot(i,j_end,3)-zbot(i,j_end-1,3))
        dqdy_xy(i,j_end-1,2) = ztop/(ztop-zbot(i,j_end-1,2)) &
                      * rdy*(zbot(i,j_end-1,2)-zbot(i,j_end-2,2))
        dqdy_xy(i,j_end-1,3) = ztop/(ztop-zbot(i,j_end-1,3)) &
                * (1./6.)*rdy*( 8.*(zbot(i,j_end-1,1)-zbot(i,j_end-2,1)) &
                      - (zbot(i,j_end,3)-zbot(i,j_end-2,3)) ) 
      enddo

      do j=j_start+1,j_end-2
      do i=i_start,i_end
        dqdy_xy(i,j,1) = ztop/(ztop-zbot(i,j,1)) &
                * (1./6.)*rdy*( 8.*(zbot(i,j+1,3)-zbot(i,j,3)) &
                      - (zbot(i,j+1,1)-zbot(i,j-1,1)) )
        dqdy_xy(i,j,2) = ztop/(ztop-zbot(i,j,2)) &
                      * 0.5*rdy*(zbot(i,j+1,2)-zbot(i,j-1,2)) &
        dqdy_xy(i,j,3) = ztop/(ztop-zbot(i,j,3)) &
                * (1./6.)*rdy*( 8.*(zbot(i,j,1)-zbot(i,j-1,1)) &
                      - (zbot(i,j+1,3)-zbot(i,j-1,3)) )
      enddo
      enddo
!----------------------------------------------------------------------------
! z dependent parts of dqdx and dqdy (therm heights first, then w heights)
!----------------------------------------------------------------------------
      do k=1,nz
        dqdx_zt(k) = (sgz(k)-ztop)/ztop
        dqdy_zt(k) = (sgz(k)-ztop)/ztop
        dqdx_zw(k) = (wgz(k)-ztop)/ztop
        dqdy_zw(k) = (wgz(k)-ztop)/ztop
      enddo
!----------------------------------------------------------------------------
! set up interpolations for surface stress definitions.  This involves
!   finding the intersection of the terrain normal with the first
!   full coordinate level.  For simplicity, we solve this by bisection.
!----------------------------------------------------------------------------
      if (viscous) then 

        x_scale  = rconst(32)
        y_scale  = rconst(33)
        z_scale  = rconst(34)
        x_0      = rconst(35)
        y_0      = rconst(36)
        z_0      = rconst(37)
        x_scale2 = rconst(38)
        y_scale2 = rconst(39)
        y_len    = rconst(40)
        tol = 0.01*dzbot
        maxit = 12 
!
! First the u pts.  
!
        do j=j_start,j_end-1
        do i=i_start,i_end

          xsrf = gx0 + (float(i)-1.)*dx  
          ysrf = gy0 + (float(j)-0.5)*dy
          zsrf = zbot(i,j,2)
          gamm = dqdz(i,j,2)*sqrt( (dqdx_xy(i,j,2)*dzdq(i,j,2))**2 &
                     + (dqdy_xy(i,j,2)*dzdq(i,j,2))**2 + 1. )
          n1 = dqdx_xy(i,j,2)*dqdx_zw(1)/gamm
          n2 = dqdy_xy(i,j,2)*dqdy_zw(1)/gamm
          n3 = dqdz(i,j,2)/gamm 
 
          if ((abs(n1).lt.1.e-4).and.(abs(n2).lt.1.e-4)) then

            x = xsrf
            y = ysrf
            dsu(i,j) = dzdq(i,j,2)*dzbot

          else

            z = zsrf + dzbot
            del = dzbot
            iconv = 0
            n = 0          
 750                if ((iconv.eq.0).and.(n.lt.maxit)) then 
              z = z - del/2.
              x = (z-zsrf)*n1/n3 + xsrf
              y = (z-zsrf)*n2/n3 + ysrf
              h0 = h_surf(x-x_0,y-y_0,x_scale,y_scale,z_scale, &
                    x_scale2,y_scale2,y_len,z_0,type_terrain)
              res = (z-h0)/(ztop-h0)*ztop - dzbot
              if (abs(res).lt.tol) then
                iconv = 1
              else 
                if (res.lt.0.) z = z + del/2.
                del = del/2.
                n = n+1
              end if
              goto 750
            endif
            if (iconv.eq.0) then
              write(6,*) 'Error:  surface stress intercept solver ', &
                             'failed to converge.'
              write(6,*) 'Stopping.'
              stop
            end if
            dsu(i,j) = sqrt((x-xsrf)**2+(y-ysrf)**2+(z-zsrf)**2)

          end if

          xsdu(i,j) = x
          ysdu(i,j) = y

        enddo
        enddo
!
! Now the v pts.  
!
        do j=j_start,j_end
        do i=1_start,i_end-1

          xsrf = gx0 + (float(i)-0.5)*dx  
          ysrf = gy0 + (float(j)-1.)*dy
          zsrf = zbot(i,j,3)
          gamm = dqdz(i,j,3)*sqrt( (dqdx_xy(i,j,3)*dzdq(i,j,3))**2 &
                     + (dqdy_xy(i,j,3)*dzdq(i,j,3))**2 + 1. )
          n1 = dqdx_xy(i,j,3)*dqdx_zw(1)/gamm
          n2 = dqdy_xy(i,j,3)*dqdy_zw(1)/gamm
          n3 = dqdz(i,j,3)/gamm 
 
          if ((abs(n1).lt.1.e-4).and.(abs(n2).lt.1.e-4)) then

            x = xsrf
            y = ysrf
            dsv(i,j) = dzdq(i,j,3)*dzbot

          else

            z = zsrf + dzbot
            del = dzbot
            iconv = 0
            n = 0          
 752                if ((iconv.eq.0).and.(n.lt.maxit)) then 
              z = z - del/2.
              x = (z-zsrf)*n1/n3 + xsrf
              y = (z-zsrf)*n2/n3 + ysrf
              h0 = h_surf(x-x_0,y-y_0,x_scale,y_scale,z_scale, &
                    x_scale2,y_scale2,y_len,z_0,type_terrain)
              res = (z-h0)/(ztop-h0)*ztop - dzbot
              if (abs(res).lt.tol) then
                iconv = 1
              else 
                if (res.lt.0.) z = z + del/2.
                del = del/2.
                n = n+1
              end if
              goto 752
            endif
            if (iconv.eq.0) then
              write(6,*) 'Error:  surface stress intercept solver ', &
                             'failed to converge.'
              write(6,*) 'Stopping.'
              stop
            end if
            dsv(i,j) = sqrt((x-xsrf)**2+(y-ysrf)**2+(z-zsrf)**2)

          end if

          xsdv(i,j) = x
          ysdv(i,j) = y

        enddo
        enddo
!
! done with surface intercept calculations
!
      end if

    END SUBROUTINE Epifanio

!=======================================================================

      subroutine slip_coeff(ucof,vcof,nx,ny)

      include 'comblk.h'
!       
! compute coefficients for tangential stress compuations at lower
!   boundary.  Some notes:
!
!     dqdx_xy*dzdq = dhdx;  dqdy_xy*dzdq = dhdy;  dqdz = zt/(zt-zs)
!
!   ucof are at u pts and vcof are at v pts.  We'll get outer values
!   through one of:  a) periodicity; b) twod; or c) zero terrain at
!   boundaries.  So we do not worry about outer values here.
!
      integer nx,ny
      real ucof(nx,ny,10),vcof(nx,ny,10)

      integer i,j
      real n1,n2,n3,gamm
!
! u pt coeffs first
!
      do j=2,ny-2
      do i=2,nx-1

        gamm = dqdz(i,j,2)*sqrt( (dqdx_xy(i,j,2)*dzdq(i,j,2))**2
     >             + (dqdy_xy(i,j,2)*dzdq(i,j,2))**2 + 1. )
        n1 = dqdx_xy(i,j,2)*dqdx_zw(1)/gamm
        n2 = dqdy_xy(i,j,2)*dqdy_zw(1)/gamm
        n3 = dqdz(i,j,2)/gamm

        ucof(i,j,1) = rdx*n1*n3
     >      + 0.5*rdx*(n3*n3-n1*n1)*dqdx_xy(i+1,j,2)*dzdq(i+1,j,2)
        ucof(i,j,2) = -rdx*n1*n3
     >      - 0.5*rdx*(n3*n3-n1*n1)*dqdx_xy(i-1,j,2)*dzdq(i-1,j,2)
        ucof(i,j,3) = 0.5*rdy*n2*n3
     >      - 0.5*rdy*n2*n1*dqdx_xy(i,j+1,2)*dzdq(i,j+1,2)
        ucof(i,j,4) = -0.5*rdy*n2*n3
     >      + 0.5*rdy*n2*n1*dqdx_xy(i,j-1,2)*dzdq(i,j-1,2)
        ucof(i,j,5) = -rdz*mfe(1)*gamm*n3
     >      + rdz*mfc(1)*gamm*n1*dqdx_xy(i,j,2)*dzdq(i,j,2)
        ucof(i,j,6) = 0.5*rdx*n2*n3
     >      + 0.5*rdx*(n3*n3-n1*n1)*dqdy_xy(i,j,3)*dzdq(i,j,3)
     >        + 0.5*rdy*n2*n1*dqdy_xy(i,j,3)*dzdq(i,j,3)
     >          + 0.25*rdz*mfc(1)*gamm*n1*dqdy_xy(i,j,2)*dzdq(i,j,2)
        ucof(i,j,7) = 0.5*rdx*n2*n3
     >      + 0.5*rdx*(n3*n3-n1*n1)*dqdy_xy(i,j+1,3)*dzdq(i,j+1,3)
     >        - 0.5*rdy*n2*n1*dqdy_xy(i,j+1,3)*dzdq(i,j+1,3)
     >          + 0.25*rdz*mfc(1)*gamm*n1*dqdy_xy(i,j,2)*dzdq(i,j,2)
        ucof(i,j,8) = -0.5*rdx*n2*n3
     >      - 0.5*rdx*(n3*n3-n1*n1)*dqdy_xy(i-1,j+1,3)*dzdq(i-1,j+1,3)
     >        - 0.5*rdy*n2*n1*dqdy_xy(i-1,j+1,3)*dzdq(i-1,j+1,3)
     >          + 0.25*rdz*mfc(1)*gamm*n1*dqdy_xy(i,j,2)*dzdq(i,j,2)
        ucof(i,j,9) = -0.5*rdx*n2*n3
     >      - 0.5*rdx*(n3*n3-n1*n1)*dqdy_xy(i-1,j,3)*dzdq(i-1,j,3)
     >        + 0.5*rdy*n2*n1*dqdy_xy(i-1,j,3)*dzdq(i-1,j,3)
     >          + 0.25*rdz*mfc(1)*gamm*n1*dqdy_xy(i,j,2)*dzdq(i,j,2)

      enddo
      enddo
!
! now v pts
!
      do j=2,ny-1
      do i=2,nx-2

        gamm = dqdz(i,j,3)*sqrt( (dqdx_xy(i,j,3)*dzdq(i,j,3))**2
     >             + (dqdy_xy(i,j,3)*dzdq(i,j,3))**2 + 1. )
        n1 = dqdx_xy(i,j,3)*dqdx_zw(1)/gamm
        n2 = dqdy_xy(i,j,3)*dqdy_zw(1)/gamm
        n3 = dqdz(i,j,3)/gamm

        vcof(i,j,1) = 0.5*rdx*n1*n3
     >      - 0.5*rdx*n1*n2*dqdy_xy(i+1,j,3)*dzdq(i+1,j,3)
        vcof(i,j,2) = -0.5*rdx*n1*n3
     >      + 0.5*rdx*n1*n2*dqdy_xy(i-1,j,3)*dzdq(i-1,j,3)
        vcof(i,j,3) = rdy*n3*n2
     >      + 0.5*rdy*(n3*n3-n2*n2)*dqdy_xy(i,j+1,3)*dzdq(i,j+1,3)
        vcof(i,j,4) = -rdy*n3*n2
     >      - 0.5*rdy*(n3*n3-n2*n2)*dqdy_xy(i,j-1,3)*dzdq(i,j-1,3)
        vcof(i,j,5) = -rdz*mfe(1)*gamm*n3
     >      + rdz*mfc(1)*gamm*n2*dqdy_xy(i,j,3)*dzdq(i,j,3)
        vcof(i,j,6) = -0.5*rdy*n1*n3
     >      - 0.5*rdy*(n3*n3-n2*n2)*dqdx_xy(i+1,j-1,2)*dzdq(i+1,j-1,2)                                                                                     
     >        - 0.5*rdx*n1*n2*dqdx_xy(i+1,j-1,2)*dzdq(i+1,j-1,2)
     >          + 0.25*rdz*mfc(1)*gamm*n2*dqdx_xy(i,j,3)*dzdq(i,j,3)
        vcof(i,j,7) = 0.5*rdy*n1*n3
     >      + 0.5*rdy*(n3*n3-n2*n2)*dqdx_xy(i+1,j,2)*dzdq(i+1,j,2)
     >        - 0.5*rdx*n1*n2*dqdx_xy(i+1,j,2)*dzdq(i+1,j,2)
     >          + 0.25*rdz*mfc(1)*gamm*n2*dqdx_xy(i,j,3)*dzdq(i,j,3)
        vcof(i,j,8) = 0.5*rdy*n1*n3
     >      + 0.5*rdy*(n3*n3-n2*n2)*dqdx_xy(i,j,2)*dzdq(i,j,2)
     >        + 0.5*rdx*n1*n2*dqdx_xy(i,j,2)*dzdq(i,j,2)
     >          + 0.25*rdz*mfc(1)*gamm*n2*dqdx_xy(i,j,3)*dzdq(i,j,3)
        vcof(i,j,9) = - 0.5*rdy*n1*n3
     >      - 0.5*rdy*(n3*n3-n2*n2)*dqdx_xy(i,j-1,2)*dzdq(i,j-1,2)
     >        + 0.5*rdx*n1*n2*dqdx_xy(i,j-1,2)*dzdq(i,j-1,2)
     >          + 0.25*rdz*mfc(1)*gamm*n2*dqdx_xy(i,j,3)*dzdq(i,j,3)

      enddo
      enddo

      return
      end

!-------------------------------------------------------------------------

      subroutine slip_solve(us,vs,resu,resv,pcru,pcrv,lpu,lpv,lru,&
                       lrv,dragx,dragy,ucof,vcof,u,v,w,u0,v0,km,&
                       zbot,time,nx,ny,nz)

      include 'comblk.h'
!
! linear solver to obtain slip/drag condition at bottom boundary.
!   Solves for us and vs as a function of the interior values of
!   u and v plus any specified drag.
!
      integer nx,ny,nz
      real us(nx,ny),vs(nx,ny),resu(nx,ny),resv(nx,ny)
      real pcru(nx,ny),pcrv(nx,ny),lpu(nx,ny),lpv(nx,ny)
      real lru(nx,ny),lrv(nx,ny),dragx(nx,ny),dragy(nx,ny)
      real ucof(nx,ny,10),vcof(nx,ny,10)
      real u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)
      real u0(nx,ny,nz),v0(nx,ny,nz),km(nx,ny,nz)
      real zbot(nx,ny,3),time

      integer i,j,k,m,ii,im1,jj,jm1,istop,maxits,iperx,ipery
      real eps,resmax,restol,beta,alpha,rhsmax
      real gamm,n1,n2,n3,kmsrf,dudz,dvdz
      double precision lplp,rlp,lplr
      parameter (maxits=500,eps=5.e-5)

      iperx = 0
      if (perx) iperx = 1
      ipery = 0
      if ((pery).and.(.not.twod)) ipery = 1
!
! set up the surface drag arrays.  First we compute any surface drag
!   associated with the background viscosity acting on the basic state
!   shear.  Then we add in any specified surface drag on top of that.
!
      do j=1,ny-1
      do i=1,nx
        gamm = dqdz(i,j,2)*sqrt( (dqdx_xy(i,j,2)*dzdq(i,j,2))**2
     >             + (dqdy_xy(i,j,2)*dzdq(i,j,2))**2 + 1. )
        n1 = dqdx_xy(i,j,2)*dqdx_zw(1)/gamm
        n2 = dqdy_xy(i,j,2)*dqdy_zw(1)/gamm
        n3 = dqdz(i,j,2)/gamm
        if (drag_opt.eq.0) then
          dragx(i,j) = -alowk0*(n3**2-n1**2)/sqrt(n3**2+n1**2)
     >        * rdz*mfe(2)*(u0(i,j,2)-u0(i,j,1))*dqdz(i,j,2)
        else if (drag_opt.eq.1) then
          dragx(i,j) = 0.
        else
          dragx(i,j) = -alowk0*(n3**2-n1**2)/sqrt(n3**2+n1**2)
     >                              * dusrf
        end if
      enddo
      enddo

      do j=1,ny
        jj = min0(j,ny-1)
        jm1 = max0(j-1,1)
      do i=1,nx-1
        gamm = dqdz(i,j,3)*sqrt( (dqdx_xy(i,j,3)*dzdq(i,j,3))**2
     >             + (dqdy_xy(i,j,3)*dzdq(i,j,3))**2 + 1. )
        n1 = dqdx_xy(i,j,3)*dqdx_zw(1)/gamm
        n2 = dqdy_xy(i,j,3)*dqdy_zw(1)/gamm
        n3 = dqdz(i,j,3)/gamm
        if (drag_opt.eq.0) then
          dragy(i,j) = alowk0*n1*n2/sqrt(n3**2+n2**2)
     >    * 0.25*rdz*mfe(2)*( u0(i+1,jm1,2)+u0(i+1,jj,2)+u0(i,jj,2)
     >          + u0(i,jm1,2)-u0(i+1,jm1,1)-u0(i+1,jj,1)-u0(i,jj,1)
     >          - u0(i,jm1,1) ) * dqdz(i,j,3)
        else if (drag_opt.eq.1) then
          dragy(i,j) = 0.
        else
          dragy(i,j) = alowk0*n1*n2/sqrt(n3**2+n2**2)*dusrf
        end if
      enddo
      enddo

      if (drag_opt.gt.2)
     >     call surf_drag(dragx,dragy,u,v,w,zbot,time,nx,ny,nz)
!
! set initial guess for us and vs by neglecting the terrain
!
      do j=1,ny-1
      do i=1,nx
        ii = min0(i,nx-1)
        im1 = max0(i-1,1)
        kmsrf = 0.5*(km(ii,j,1)+km(im1,j,1))
        dudz = -dragx(i,j)/kmsrf
        us(i,j) = u(i,j,1) - dudz/(rdz*mfe(1))*dzdq(i,j,2)
      enddo
      enddo
      do j=1,ny
      do i=1,nx-1
        jj = min0(j,ny-1)
        jm1 = max0(j-1,1)
        kmsrf = 0.5*(km(i,jj,1)+km(i,jm1,1))
        dvdz = -dragy(i,j)/kmsrf
        vs(i,j) = v(i,j,1) - dvdz/(rdz*mfe(1))*dzdq(i,j,3)
      enddo
      enddo
!
! compute rhs.  This includes any specified drag.
!
      rhsmax = 0.

      do j=2,ny-2
      do i=2,nx-1

        gamm = dqdz(i,j,2)*sqrt( (dqdx_xy(i,j,2)*dzdq(i,j,2))**2
     >             + (dqdy_xy(i,j,2)*dzdq(i,j,2))**2 + 1. )
        n1 = dqdx_xy(i,j,2)*dqdx_zw(1)/gamm
        n2 = dqdy_xy(i,j,2)*dqdy_zw(1)/gamm
        n3 = dqdz(i,j,2)/gamm
!       kmsrf = 0.75*(km(i,j,1)+km(i-1,j,1))
!    >                  - 0.25*(km(i,j,2)+km(i-1,j,2))
        kmsrf = 0.5*(km(i,j,1)+km(i-1,j,1))
        ucof(i,j,10) = -dragx(i,j)/kmsrf
     >                   * sqrt(n3**2+n1**2)
     >             - rdz*mfe(1)*gamm*n3*u(i,j,1)
     >        + 0.5*rdz*mfc(1)*gamm*n1*(w(i,j,2)+w(i-1,j,2))
        if (abs(ucof(i,j,10)).gt.rhsmax) rhsmax = abs(ucof(i,j,10))

      enddo
      enddo

      do j=2,ny-1
      do i=2,nx-2

        gamm = dqdz(i,j,3)*sqrt( (dqdx_xy(i,j,3)*dzdq(i,j,3))**2
     >             + (dqdy_xy(i,j,3)*dzdq(i,j,3))**2 + 1. )
        n1 = dqdx_xy(i,j,3)*dqdx_zw(1)/gamm
        n2 = dqdy_xy(i,j,3)*dqdy_zw(1)/gamm
        n3 = dqdz(i,j,3)/gamm
!       kmsrf = 0.75*(km(i,j,1)+km(i,j-1,1))
!    >               - 0.25*(km(i,j,2)+km(i,j-1,2))
        kmsrf = 0.5*(km(i,j,1)+km(i,j-1,1))
        vcof(i,j,10) = -dragy(i,j)/kmsrf
     >                     * sqrt(n3**2+n2**2)
     >             - rdz*mfe(1)*gamm*n3*v(i,j,1)
     >        + 0.5*rdz*mfc(1)*gamm*n2*(w(i,j,2)+w(i,j-1,2))
        if (abs(vcof(i,j,10)).gt.rhsmax) rhsmax = abs(vcof(i,j,10))

      enddo
      enddo

      restol = eps*rhsmax
!
! make sure initial us and vs satisfy bcs
!
      call set_bcs_slip(us,vs,u,v,km,dragx,dragy,nx,ny,nz,1)
!
! compute initial residual, u pts and then v pts
!
      resmax = 0.

      do j=2+ipery,ny-2-ipery
      do i=2+2*iperx,nx-1-iperx
        resu(i,j) = ucof(i,j,1)*us(i+1,j) + ucof(i,j,2)*us(i-1,j)
     >            + ucof(i,j,3)*us(i,j+1) + ucof(i,j,4)*us(i,j-1)
     >            + ucof(i,j,5)*us(i,j) + ucof(i,j,6)*vs(i,j)
     >            + ucof(i,j,7)*vs(i,j+1) + ucof(i,j,8)*vs(i-1,j+1)
     >            + ucof(i,j,9)*vs(i-1,j) - ucof(i,j,10)
        if (abs(resu(i,j)).gt.resmax) resmax = abs(resu(i,j))
      enddo
      enddo

      do j=2+2*ipery,ny-1-ipery
      do i=2+iperx,nx-2-iperx
        resv(i,j) = vcof(i,j,1)*vs(i+1,j) + vcof(i,j,2)*vs(i-1,j)
     >            + vcof(i,j,3)*vs(i,j+1) + vcof(i,j,4)*vs(i,j-1)
     >            + vcof(i,j,5)*vs(i,j) + vcof(i,j,6)*us(i+1,j-1)
     >            + vcof(i,j,7)*us(i+1,j) + vcof(i,j,8)*us(i,j)
     >            + vcof(i,j,9)*us(i,j-1) - vcof(i,j,10)
        if (abs(resv(i,j)).gt.resmax) resmax = abs(resv(i,j))
      enddo
      enddo

      if (resmax.le.restol) return
!
! set bcs on initial residual
!                                                                                                                                                          
      call set_bcs_slip(resu,resv,u,v,km,dragx,dragy,nx,ny,nz,0)
!
! initialize lpu,lpv,pcru,pcrv
!                                                                                                                                                          
      do j=2+ipery,ny-2-ipery
      do i=2+2*iperx,nx-1-iperx
        pcru(i,j) = resu(i,j)
        lpu(i,j) =
     >      ucof(i,j,1)*resu(i+1,j) + ucof(i,j,2)*resu(i-1,j)
     >    + ucof(i,j,3)*resu(i,j+1) + ucof(i,j,4)*resu(i,j-1)
     >    + ucof(i,j,5)*resu(i,j) + ucof(i,j,6)*resv(i,j)
     >    + ucof(i,j,7)*resv(i,j+1) + ucof(i,j,8)*resv(i-1,j+1)
     >    + ucof(i,j,9)*resv(i-1,j)
      enddo
      enddo

      do j=2+2*ipery,ny-1-ipery
      do i=2+iperx,nx-2-iperx
        pcrv(i,j) = resv(i,j)
        lpv(i,j) =
     >      vcof(i,j,1)*resv(i+1,j) + vcof(i,j,2)*resv(i-1,j)
     >    + vcof(i,j,3)*resv(i,j+1) + vcof(i,j,4)*resv(i,j-1)
     >    + vcof(i,j,5)*resv(i,j) + vcof(i,j,6)*resu(i+1,j-1)
     >    + vcof(i,j,7)*resu(i+1,j) + vcof(i,j,8)*resu(i,j)
     >    + vcof(i,j,9)*resu(i,j-1)
      enddo
      enddo
!
! done initializing; start main iteration loop
!                                                                                                                                                          
      istop = 0
      n = 1

 997    if ((n.le.maxits).and.(istop.eq.0)) then
! 
! compute beta
!                                                                                                                                                          
        rlp = 0.
        lplp = 0.

        do j=2+ipery,ny-2-ipery
        do i=2+2*iperx,nx-1-iperx
          rlp = rlp + resu(i,j)*lpu(i,j)
          lplp = lplp + lpu(i,j)*lpu(i,j)
        enddo
        enddo

        do j=2+2*ipery,ny-1-ipery
        do i=2+iperx,nx-2-iperx
          rlp = rlp + resv(i,j)*lpv(i,j)
          lplp = lplp + lpv(i,j)*lpv(i,j)
        enddo
        enddo

        beta = -real(rlp)/real(lplp)
!
! advance (us,vs) and (resu,resv)
!                                                                                                                                                          
        resmax = 0.

        do j=2+ipery,ny-2-ipery
        do i=2+2*iperx,nx-1-iperx
          us(i,j) = us(i,j) + beta*pcru(i,j)
          resu(i,j) = resu(i,j) + beta*lpu(i,j)
          if (abs(resu(i,j)).gt.resmax) resmax = abs(resu(i,j))
        enddo
        enddo

        do j=2+2*ipery,ny-1-ipery
        do i=2+iperx,nx-2-iperx
          vs(i,j) = vs(i,j) + beta*pcrv(i,j)
          resv(i,j) = resv(i,j) + beta*lpv(i,j)
          if (abs(resv(i,j)).gt.resmax) resmax = abs(resv(i,j))
        enddo
        enddo
!
! set bcs on res
!                                                                                                                                                          
        call set_bcs_slip(resu,resv,u,v,km,dragx,dragy,nx,ny,nz,0)
!
! check for convergence and continue
!                                                                                                                                                          
        if (resmax.lt.restol) then

          istop = 1

        else
!
! compute lru,lrv
!                                                                                                                                                          
          do j=2+ipery,ny-2-ipery
          do i=2+2*iperx,nx-1-iperx
            lru(i,j) =
     >        ucof(i,j,1)*resu(i+1,j) + ucof(i,j,2)*resu(i-1,j)
     >      + ucof(i,j,3)*resu(i,j+1) + ucof(i,j,4)*resu(i,j-1)
     >      + ucof(i,j,5)*resu(i,j) + ucof(i,j,6)*resv(i,j)
     >      + ucof(i,j,7)*resv(i,j+1) + ucof(i,j,8)*resv(i-1,j+1)
     >      + ucof(i,j,9)*resv(i-1,j)
          enddo
          enddo

          do j=2+2*ipery,ny-1-ipery
          do i=2+iperx,nx-2-iperx
            lrv(i,j) =
     >        vcof(i,j,1)*resv(i+1,j) + vcof(i,j,2)*resv(i-1,j)
     >      + vcof(i,j,3)*resv(i,j+1) + vcof(i,j,4)*resv(i,j-1)
     >      + vcof(i,j,5)*resv(i,j) + vcof(i,j,6)*resu(i+1,j-1)
     >      + vcof(i,j,7)*resu(i+1,j) + vcof(i,j,8)*resu(i,j)
     >      + vcof(i,j,9)*resu(i,j-1)
          enddo
          enddo
!
! compute alpha
!                                                                                                                                                          
          lplr = 0.

          do j=2+ipery,ny-2-ipery
          do i=2+2*iperx,nx-1-iperx
            lplr = lplr + lpu(i,j)*lru(i,j)
          enddo
          enddo

          do j=2+2*ipery,ny-1-ipery
          do i=2+iperx,nx-2-iperx
            lplr = lplr + lpv(i,j)*lrv(i,j)
          enddo
          enddo

          alpha = -real(lplr)/real(lplp)
!
! advance (pcru,pcrv) and (lpu,lpv)
!
          do j=2+ipery,ny-2-ipery
          do i=2+2*iperx,nx-1-iperx
            pcru(i,j) = resu(i,j) + alpha*pcru(i,j)
            lpu(i,j) = lru(i,j) + alpha*lpu(i,j)
          enddo
          enddo

          do j=2+2*ipery,ny-1-ipery
          do i=2+iperx,nx-2-iperx
            pcrv(i,j) = resv(i,j) + alpha*pcrv(i,j)
            lpv(i,j) = lrv(i,j) + alpha*lpv(i,j)
          enddo
          enddo

        end if

        n = n+1
        goto 997

      end if
!
! check for success
!
      if (istop.eq.0) then
        write(6,*)
        write(6,*) 'ERROR:  iterative free-slip solver failed.'
 767        FORMAT(' max frac residual = ',F8.3)
        write(6,767) resmax/restol
        write(6,*) 'Stopping.'
        stop
      end if
!
! set bcs on us and vs and set w at surface
!
      call set_bcs_slip(us,vs,u,v,km,dragx,dragy,nx,ny,nz,1)

      do j=1,ny-1
      do i=1,nx-1
        w(i,j,1) =
     >     0.5*(us(i,j)+us(i+1,j))*dqdx_xy(i,j,1)*dzdq(i,j,1)
     >   + 0.5*(vs(i,j)+vs(i,j+1))*dqdy_xy(i,j,1)*dzdq(i,j,1)
      enddo
      enddo

      return
      end

!-----------------------------------------------------------------------      

      subroutine set_bcs_slip(fus,fvs,u,v,km,dragx,dragy,nx,ny,nz,ius)

!      include 'comblk.h'
! 
! set lateral boundary conditions for slip/drag lower boundary solver.
!   We assume the boundaries are one of:  i) periodic; ii) twod; or
!   iii) no terrain at boundary.  fus and fvs are either us and vs or 
!   the corresponding residuals.  The bcs on the residuals are    
!   homogenous versions of the bcs on us and vs.  
!                                                                                                                                                          
      integer nx,ny,nz,ius
      real fus(nx,ny),fvs(nx,ny),u(nx,ny,nz),v(nx,ny,nz),km(nx,ny,nz)
      real dragx(nx,ny),dragy(nx,ny)

      integer i,j
      real dudz,dvdz,kmsrf

      if (twod) then
        do i=1,nx
          fus(i,1) = fus(i,2)
          fus(i,ny-1) = fus(i,ny-2)
        enddo
        do i=1,nx-1
          fvs(i,1) = fvs(i,2)
          fvs(i,ny) = fvs(i,ny-1)
        enddo
      else if (pery) then
        do i=1,nx
          fus(i,1) = fus(i,ny-4)
          fus(i,2) = fus(i,ny-3)
          fus(i,ny-2) = fus(i,3)
          fus(i,ny-1) = fus(i,4)
        enddo
        do i=1,nx-1
          fvs(i,1) = fvs(i,ny-4)
          fvs(i,2) = fvs(i,ny-3)
          fvs(i,3) = fvs(i,ny-2)
          fvs(i,ny-1) = fvs(i,4)
          fvs(i,ny) = fvs(i,5)
        enddo
      else
        do i=1,nx
          ii = min0(i,nx-1)
          im1 = max0(i-1,1)
          kmsrf = 0.5*(km(ii,1,1)+km(im1,1,1))
          dudz = -dragx(i,1)/kmsrf
          fus(i,1) = float(ius)*(u(i,1,1) - dudz/(rdz*mfe(1)))
          kmsrf = 0.5*(km(ii,ny-1,1)+km(im1,ny-1,1))
          dudz = -dragx(i,ny-1)/kmsrf
          fus(i,ny-1) = float(ius)*(u(i,ny-1,1) - dudz/(rdz*mfe(1)))
        enddo
        do i=1,nx-1
          kmsrf = km(i,1,1)
          dvdz = -dragy(i,1)/kmsrf
          fvs(i,1) = float(ius)*(v(i,1,1) - dvdz/(rdz*mfe(1)))
          kmsrf = km(i,ny-1,1)
          dvdz = -dragy(i,ny)/kmsrf
          fvs(i,ny) = float(ius)*(v(i,ny,1) - dvdz/(rdz*mfe(1)))
        enddo
      end if

      if (perx) then
        do j=1,ny-1
          fus(1,j) = fus(nx-4,j)
          fus(2,j) = fus(nx-3,j)
          fus(3,j) = fus(nx-2,j)
          fus(nx-1,j) = fus(4,j)
          fus(nx,j) = fus(5,j)
        enddo
        do j=1,ny
          fvs(1,j) = fvs(nx-4,j)
          fvs(2,j) = fvs(nx-3,j)
          fvs(nx-2,j) = fvs(3,j)
          fvs(nx-1,j) = fvs(4,j)
        enddo
      else
        do j=1,ny-1
          kmsrf = km(1,j,1)
          dudz = -dragx(1,j)/kmsrf
          fus(1,j) = float(ius)*(u(1,j,1) - dudz/(rdz*mfe(1)))
          kmsrf = km(nx-1,j,1)
          dudz = -dragx(nx,j)/kmsrf
          fus(nx,j) = float(ius)*(u(nx,j,1) - dudz/(rdz*mfe(1)))
        enddo
        do j=1,ny
          jj = min0(j,ny-1)
          jm1 = max0(j-1,1)
          kmsrf = 0.5*(km(1,jj,1)+km(1,jm1,1))
          dvdz = -dragy(1,j)/kmsrf
          fvs(1,j) = float(ius)*(v(1,j,1) - dvdz/(rdz*mfe(1)))
          kmsrf = 0.5*(km(nx-1,jj,1)+km(nx-1,jm1,1))
          dvdz = -dragy(nx-1,j)/kmsrf
          fvs(nx-1,j) = float(ius)*(v(nx-1,j,1) - dvdz/(rdz*mfe(1)))
        enddo
      end if

      return
      end subroutine slip

!=======================================================================

      subroutine surf_drag(dragx,dragy,u,v,w,zbot,time,nx,ny,nz)
!
! surf_drag returns the specified surface drag across the lower boundary
!   of the model domain.  This drag is then later imposed in the 
!   subroutines that compute the viscous stresses at the boundary.
! Note:  specified surface drags should always correspond to drag_opt > 2
!   (options 0 through 2 are already taken elsewhere and correspond to
!   zero specified surface drag)
!
      include 'comblk.h'

      integer nx,ny,nz
      real dragx(nx,ny),dragy(nx,ny),zbot(nx,ny,3),time
      real u(nx,ny,nz),v(nx,ny,nz),w(nx,ny,nz)

      integer i,j,id,jd
      real x,y,x0,y0,Cd,utan,vtan,uabs
      real udrg,vdrg,wdrg,un,unx,uny,unz,gamm,n1,n2,n3,drag_coeff
!
! First the upts.
!
      do j=1,ny-1
      do i=1,nx
!
! interpolate u, v and w to the intercept of the terrain normal
!  with the first full coordinate level.  
!
        x = xsdu(i,j)
        y = ysdu(i,j)
        id = max0(min0(int((x-gx0)/dx+1.),nx-1),1)
        jd = max0(min0(int((y-gy0-0.5*dy)/dy+1.),ny-2),1) 
        x0 = gx0 + (float(id)-1.)*dx
        y0 = gy0 + (float(jd)-0.5)*dy
        udrg = 0.5*rdx*rdy*(
     >      (x-x0-dx)*(y-y0-dy)*(u(id,jd,1)+u(id,jd,2))
     >        - (x-x0)*(y-y0-dy)*(u(id+1,jd,1)+u(id+1,jd,2))
     >           - (x-x0-dx)*(y-y0)*(u(id,jd+1,1)+u(id,jd+1,2))
     >               + (x-x0)*(y-y0)*(u(id+1,jd+1,1)+u(id+1,jd+1,2)) )
        id = max0(min0(int((x-gx0-0.5*dx)/dx+1.),nx-2),1)
        jd = max0(min0(int((y-gy0)/dy+1.),ny-1),1) 
        x0 = gx0 + (float(id)-0.5)*dx
        y0 = gy0 + (float(jd)-1.)*dy
        vdrg = 0.5*rdx*rdy*(
     >    (x-x0-dx)*(y-y0-dy)*(v(id,jd,1)+v(id,jd,2))
     >        - (x-x0)*(y-y0-dy)*(v(id+1,jd,1)+v(id+1,jd,2))
     >           - (x-x0-dx)*(y-y0)*(v(id,jd+1,1)+v(id,jd+1,2))
     >               + (x-x0)*(y-y0)*(v(id+1,jd+1,1)+v(id+1,jd+1,2)) )
        id = max0(min0(int((x-gx0-0.5*dx)/dx+1.),nx-2),1)
        jd = max0(min0(int((y-gy0-0.5*dy)/dy+1.),ny-2),1)
        x0 = gx0 + (float(id)-0.5)*dx
        y0 = gy0 + (float(jd)-0.5)*dy
        wdrg = rdx*rdy*(
     >          (x-x0-dx)*(y-y0-dy)*w(id,jd,2)
     >             - (x-x0)*(y-y0-dy)*w(id+1,jd,2)
     >                 - (x-x0-dx)*(y-y0)*w(id,jd+1,2)
     >                     + (x-x0)*(y-y0)*w(id+1,jd+1,2) )
!
! next compute the tangential winds
!
        gamm = dqdz(i,j,2)*sqrt( (dqdx_xy(i,j,2)*dzdq(i,j,2))**2
     >               + (dqdy_xy(i,j,2)*dzdq(i,j,2))**2 + 1. )
        n1 = dqdx_xy(i,j,2)*dqdx_zw(1)/gamm
        n2 = dqdy_xy(i,j,2)*dqdy_zw(1)/gamm
        n3 = dqdz(i,j,2)/gamm
        utan = (udrg*n3-wdrg*n1)/sqrt(n3**2+n1**2)
        un = udrg*n1+vdrg*n2+wdrg*n3
        unx = un*n1
        uny = un*n2
        unz = un*n3
        uabs = sqrt((udrg-unx)**2+(vdrg-uny)**2+(wdrg-unz)**2)
!
! compute the drag coefficient
!
        x = gx0 + (float(i)-1.)*dx - sdxc
        y = gy0 + (float(j)-0.5)*dy - sdyc
        Cd = drag_coeff(x,y,dsu(i,j),cdmax,cdmin,sdlx,sdly,drag_opt)
!
! finally the drag
!
        dragx(i,j) = dragx(i,j) - Cd*uabs*utan 

      enddo
      enddo
!
! Now for the v pts.
!
      do j=1,ny
      do i=1,nx-1
!
! interpolate u, v and w to the intercept points
!
        x = xsdv(i,j)
        y = ysdv(i,j)
        id = max0(min0(int((x-gx0)/dx+1.),nx-1),1)
        jd = max0(min0(int((y-gy0-0.5*dy)/dy+1.),ny-2),1) 
        x0 = gx0 + (float(id)-1.)*dx
        y0 = gy0 + (float(jd)-0.5)*dy
        udrg = 0.5*rdx*rdy*(
     >    (x-x0-dx)*(y-y0-dy)*(u(id,jd,1)+u(id,jd,2))
     >        - (x-x0)*(y-y0-dy)*(u(id+1,jd,1)+u(id+1,jd,2))
     >           - (x-x0-dx)*(y-y0)*(u(id,jd+1,1)+u(id,jd+1,2))
     >               + (x-x0)*(y-y0)*(u(id+1,jd+1,1)+u(id+1,jd+1,2)) )
        id = max0(min0(int((x-gx0-0.5*dx)/dx+1.),nx-2),1)
        jd = max0(min0(int((y-gy0)/dy+1.),ny-1),1) 
        x0 = gx0 + (float(id)-0.5)*dx
        y0 = gy0 + (float(jd)-1.)*dy
        vdrg = 0.5*rdx*rdy*(
     >    (x-x0-dx)*(y-y0-dy)*(v(id,jd,1)+v(id,jd,2))
     >        - (x-x0)*(y-y0-dy)*(v(id+1,jd,1)+v(id+1,jd,2))
     >           - (x-x0-dx)*(y-y0)*(v(id,jd+1,1)+v(id,jd+1,2))
     >               + (x-x0)*(y-y0)*(v(id+1,jd+1,1)+v(id+1,jd+1,2)) )
        id = max0(min0(int((x-gx0-0.5*dx)/dx+1.),nx-2),1)
        jd = max0(min0(int((y-gy0-0.5*dy)/dy+1.),ny-2),1)
        x0 = gx0 + (float(id)-0.5)*dx
        y0 = gy0 + (float(jd)-0.5)*dy
        wdrg = rdx*rdy*(
     >          (x-x0-dx)*(y-y0-dy)*w(id,jd,2)
     >             - (x-x0)*(y-y0-dy)*w(id+1,jd,2)
     >                 - (x-x0-dx)*(y-y0)*w(id,jd+1,2)
     >                     + (x-x0)*(y-y0)*w(id+1,jd+1,2) )
!
! next the tangential winds
!
        gamm = dqdz(i,j,3)*sqrt( (dqdx_xy(i,j,3)*dzdq(i,j,3))**2
     >             + (dqdy_xy(i,j,3)*dzdq(i,j,3))**2 + 1. )
        n1 = dqdx_xy(i,j,3)*dqdx_zw(1)/gamm
        n2 = dqdy_xy(i,j,3)*dqdy_zw(1)/gamm
        n3 = dqdz(i,j,3)/gamm
        vtan = (vdrg*n3-wdrg*n2)/sqrt(n3**2+n2**2)
        un = udrg*n1+vdrg*n2+wdrg*n3
        unx = un*n1
        uny = un*n2
        unz = un*n3
        uabs = sqrt((udrg-unx)**2+(vdrg-uny)**2+(wdrg-unz)**2)
!
! drag coefficient
!
        x = gx0 + (float(i)-0.5)*dx - sdxc
        y = gy0 + (float(j)-1.)*dy - sdyc
        Cd = drag_coeff(x,y,dsv(i,j),cdmax,cdmin,sdlx,sdly,drag_opt)
!
! and finally the drag
!
        dragy(i,j) = dragy(i,j) - Cd*uabs*vtan

      enddo
      enddo

      return
      end
     
!--------------------------------------------------------------------- 

      real function drag_coeff(x,y,ds,cdmax,cdmin,sdlx,sdly,drag_opt)
!
!  returns the drag coefficient according to the drag option and
!    associated parameters.  Here x and y are relative to a 
!    coordinate system with origin at (sdxc,sdyc).
!  
      integer drag_opt
      real x,y,ds,cdmax,cdmin,sdlx,sdly

      real rad,sdlx2,pi,srgh

      pi = 2.*asin(1.)
      sdlx2 = sdly
     
      if (drag_opt.eq.3) then
 
        rad = sqrt(x**2+y**2)
        if (rad.le.sdlx) then
          drag_coeff = cdmax
        else if (rad.le.(sdlx+sdlx2)) then
          drag_coeff = cdmin + 0.5*(cdmax-cdmin)
     >                       * (1.+cos(pi*(rad-sdlx)/sdlx2))
        else
          drag_coeff = cdmin
        end if

      else if (drag_opt.eq.4) then

        rad = sqrt(x**2+y**2)
        if (rad.le.sdlx) then
          srgh = cdmax
        else if (rad.le.(sdlx+sdlx2)) then
          srgh = cdmin + 0.5*(cdmax-cdmin)
     >                       * (1.+cos(pi*(rad-sdlx)/sdlx2))
        else
          srgh = cdmin
        end if
        if (srgh.lt.1.e-15) then
          drag_coeff = 0.
        else
          drag_coeff = 0.16/log(ds/srgh)**2
        end if

      else if (drag_opt.eq.5) then

        rad = abs(x)
        if (rad.le.sdlx) then
          srgh = cdmax
        else if (rad.le.(sdlx+sdlx2)) then
          srgh = cdmin + 0.5*(cdmax-cdmin) &
                            * (1.+cos(pi*(rad-sdlx)/sdlx2))
        else
          srgh = cdmin
        end if
        if (srgh.lt.1.e-15) then
          drag_coeff = 0.
        else
          drag_coeff = 0.16/log(ds/srgh)**2
        end if

      else

        write(6,*) 'Error:  drag option ',drag_opt,&
                        ' not recognized.'
        write(6,*) 'Stopping.'
        stop
    
      end if

      return
      end


!=======================================================================
!=======================================================================

    END MODULE module_3Dsfcflux

!=======================================================================
!=======================================================================
