MODULE avg_RA 
   !!==============================================================================
   !!                       ***  MODULE avg_RA   ***
   !! spatial averaging of input field
   !!==============================================================================
   !! History :        !  08-17  (R. Abel) module based on shapiro from
                      !  J.M.Molines
   !!----------------------------------------------------------------------
   !! * Modules used
   USE in_out_manager
   USE dom_oce         ! ocean space and time domain
   USE wrk_nemo        ! work arrays
   USE timing          ! Timing
   USE lbclnk

   
   IMPLICIT NONE
   PRIVATE

   PUBLIC avg_2D  ! use by sbcblk_core  and sbcssr

  CONTAINS

  SUBROUTINE avg_2D(ptabin, cd_overlap, ptabout) !GIG
      !!----------------------------------------------------------------------
      !!                  ***  routine Shapiro_1D  ***
      !!
      !! ** Purpose :  Multiple application (kiter) of a shapiro filter
      !!               on ptabin to produce ptabout.
      !!
      !! ** Method  :   
      !!
      !! ** Action  :   ptabout filtered output from ptabin
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(IN)  :: ptabin     ! input array
      CHARACTER(len=*),         INTENT(IN)  :: cd_overlap ! = one of MERCA_GLOB, REGULAR_GLOB, ORCA_GLOB (??)
      REAL(wp), DIMENSION(:,:), INTENT(OUT) :: ptabout    ! output array

      ! * Local variable
      INTEGER                               :: ji, jj, jn ! dummy loop index
      REAL(wp), POINTER, DIMENSION(:,:)     :: zvarout    ! working array
      REAL(wp)                              :: grid_sum   ! working sum
      REAL(wp), PARAMETER                   :: rp_aniso_diff_XY=2.25 !  anisotrope case (???)
                                                          ! Empirical value for 140 iterations
                                                          ! for an anisotropic ratio of 1.5.
                                                          ! (re ???)
      REAL(wp)                              :: zalphax    ! weight coeficient (x direction)
      REAL(wp)                              :: zalphay    ! weight coeficient (y direction)
      REAL(wp)                              :: znum       ! numerator
      REAL(wp)                              :: zden       ! denominator
!
!------------------------------------------------------------------------------
!
      IF( nn_timing == 1 )  CALL timing_start('avg_RA')
      !
      CALL wrk_alloc( jpi,jpj, zvarout )

!     Global ocean case
      IF (( cd_overlap == 'MERCA_GLOB' )   .OR.   &
          ( cd_overlap == 'REGULAR_GLOB' ) .OR.   &
          ( cd_overlap == 'ORCA_GLOB' )) THEN
       ptabout(:,:) = ptabin(:,:)
       zvarout(:,:) = ptabout(:,:)  ! ptabout intent out ???

       zalphax=1./2.
       zalphay=1./2.

!  Dx/Dy=rp_aniso_diff_XY  , D_ = vitesse de diffusion
!  140 passes du fitre, Lx/Ly=1.5, le rp_aniso_diff_XY correspondant est:

       IF ( rp_aniso_diff_XY >=  1. ) zalphay=zalphay/rp_aniso_diff_XY
       IF ( rp_aniso_diff_XY <   1. ) zalphax=zalphax*rp_aniso_diff_XY

!            DO jj = 2,jpjm1
!               DO ji = 2,jpim1
!                  ! We crop on the coast
!                   znum = zvarout(ji,jj)   &
!                          + 0.25*zalphax*(zvarout(ji-1,jj  )-zvarout(ji,jj))*tmask(ji-1,jj  ,1)  &
!                          + 0.25*zalphax*(zvarout(ji+1,jj  )-zvarout(ji,jj))*tmask(ji+1,jj  ,1)  &
!                          + 0.25*zalphay*(zvarout(ji  ,jj-1)-zvarout(ji,jj))*tmask(ji  ,jj-1,1)  &
!                          + 0.25*zalphay*(zvarout(ji  ,jj+1)-zvarout(ji,jj))*tmask(ji  ,jj+1,1)
!                   ptabout(ji,jj)=znum*tmask(ji,jj,1)+ptabin(ji,jj)*(1.-tmask(ji,jj,1))
!                ENDDO  ! end loop jj
!            ENDDO  ! end loop ji

           ! Averaging tailored for 4x4 boxes !Note: take care about CPU arrangement
            DO jj = nldj,nlej,4
               !WRITE(numout,*) 'nldj',nldj
               !WRITE(numout,*) 'nlej',nlej
               !WRITE(numout,*) 'jj',jj
               DO ji = nldi,nlei,4
                  ! We crop on the coast
                   grid_sum = 0.0
                   grid_sum = SUM( e1e2t(ji:ji+3,jj:jj+3) * tmask(ji:ji+3,jj:jj+3,1) )
                   !WRITE(numout,*) 'grid_sum', grid_sum
                   IF (grid_sum /= 0) THEN
                      ptabout(ji:ji+3,jj:jj+3)=SUM(zvarout(ji:ji+3,jj:jj+3)* tmask(ji:ji+3,jj:jj+3,1) &
                                                * e1e2t(ji:ji+3,jj:jj+3)  )/grid_sum 
                   ELSE
                      ptabout(ji:ji+3,jj:jj+3)=0
                   ENDIF

                                       !   / SUM(tmask(ji:ji+3,jj:jj+3,1)* e1e2t(ji:ji+3,jj:jj+3) ) 
!                   znum = zvarout(ji,jj)   &
!                          + 0.125*zvarout(ji-2,jj)*tmask(ji-2,jj  ,1)  &
!                          + 0.125*zvarout(ji-1,jj)*tmask(ji-1,jj  ,1)  &
!                          + 0.125*zvarout(ji+1,jj)*tmask(ji+1,jj  ,1)  &
!                          + 0.125*zvarout(ji+2,jj)*tmask(ji+2,jj  ,1)  &
!                          + 0.125*zvarout(ji,jj-2)*tmask(ji  ,jj-2,1)  &
!                          + 0.125*zvarout(ji,jj-1)*tmask(ji  ,jj-1,1)  &
!                          + 0.125*zvarout(ji,jj+1)*tmask(ji  ,jj+1,1)
!                          + 0.125*zvarout(ji,jj+2)*tmask(ji  ,jj+2,1)
!                   ptabout(ji,jj)=znum*tmask(ji,jj,1)+ptabin(ji,jj)*(1.-tmask(ji,jj,1))
                ENDDO  ! end loop jj
            ENDDO  ! end loop ji
!var_out[j-fact:j,i-fact:i]=np.nanmean(var[j-fact:j,i-fact:i].flatten(),0)
!
!
!           Periodical condition in case of cd_overlap (global ocean)
!           - on a mercator projection grid we consider that singular point at poles
!             are a mean of the values at points of the previous latitude
!           - on ORCA and regular grid we copy the values at points of the previous latitude
            IF ( cd_overlap == 'MERCAT_GLOB' ) THEN
!GIG case unchecked  ! JMM for sure not valid in MPP (BUG)
               ptabout(1,1) = SUM(ptabout(:,2)) / jpi
               ptabout(jpi,jpj) = SUM(ptabout(:,jpj-1)) / jpi
            ELSE
               CALL lbc_lnk(ptabout, 'T', 1.) ! Boundary condition
            ENDIF
            zvarout(:,:) = ptabout(:,:)
      ENDIF

      CALL wrk_dealloc( jpi,jpj, zvarout )
      IF( nn_timing == 1 )  CALL timing_stop('avg_RA')
!
END SUBROUTINE avg_2D     

END MODULE avg_RA
