!=======================================================================
MODULE FixedBottomSupportStructure


   ! This MODULE stores variables and routines used in these time domain
   !   hydrodynamic loading for the fixed-bottom offshore support
   !   structures.

   USE                             NWTC_Library


CONTAINS
!=======================================================================
   SUBROUTINE MorisonTwrLd ( JNode, TwrDiam, TwrCA, TwrCD, X, XD, ZTime, TwrAM, TwrFt, WaveDat, ErrStat )
!bjj: input X is not used

!JASON: WE MUST CORRECT THE CALCULATIONS IN THIS ROUTINE!: CREATE AND STORE EXTRA ARRAYS OF WAVE KINEMATICS AT THE FREE SURFACE AND SEABED.  USE THESE ARRAYS TO COMPUTE THE WAVE KINEMATICS FOR ELEMENTS THAT ARE ONLY PARTIALLY IMMERSED IN FLUID (I.E., WHERE THE WAVE ELEVATION IS JUST BELOW THE CURRENT TOWER NODE OR THE SEABED IS JUST ABOVE THE CURRENT TOWER NODE).  FOR GH BLADED WAVE DATA, FORM THESE ARRAYS USING THE HIGHEST OR LOWEST AVAILABLE WAVE DATA.  WITHOUT THIS FIX, ONE SHOULD GET ARROUND THE PROBLEM BY USING A FINE RESOLUTION OF TOWER ELEMENTS OVER THE PORTION OF THE TOWER THAT MAY END UP BEING COVERED OR UNCOVERED BY THE INCIDENT WAVES.
!JASON: USE THE UNDEFLECTED / UNDISPLACED TOWER ELEMENT LOCATION AND ORIENTATION TO BEGIN WITH; INCLUDE THE EFFECTS OF THE INSTANTANEOUS ELEMENT LOCATION AND ORIENTATION LATER!!!
       ! This routine is used to implement Morison's equation for the
       ! hydrodynamic loading on a monopile.


   USE                             Waves


   IMPLICIT                        NONE


      ! Passed Variables:

   REAL(ReKi),           INTENT(  OUT)      :: TwrAM  (6,6)                       ! Added mass matrix per unit length of current tower element (kg/m, kg-m/m, kg-m^2/m)
   REAL(ReKi),           INTENT(IN   )      :: TwrCA                              ! Normalized hydrodynamic added mass   coefficient of current tower element (-)
   REAL(ReKi),           INTENT(IN   )      :: TwrCD                              ! Normalized hydrodynamic viscous drag coefficient of current tower element (-)
   REAL(ReKi),           INTENT(IN   )      :: TwrDiam                            ! Diameter of current tower element (meters)
   REAL(ReKi),           INTENT(  OUT)      :: TwrFt    (6)                       ! The surge/xi (1), sway/yi (2), and heave/zi (3)-components of the portion of the tower force per unit length (in N/m) at the current tower element and the roll/xi (4), pitch/yi (5), and yaw/zi (6)-components of the portion of the tower moment per unit length (in N-m/m) acting at the current tower element associated with everything but the added-mass effects; positive forces are in the direction of motion.
   REAL(ReKi),           INTENT(IN   )      :: X        (6)                       ! The 3 components of the translational displacement (in m  ) of the current tower node and the 3 components of the rotational displacement       (in rad  ) of the current tower element relative to the inertial frame origin at ground level [onshore] or MSL [offshore].
   REAL(ReKi),           INTENT(IN   )      :: XD       (6)                       ! The 3 components of the translational velocity     (in m/s) of the current tower node and the 3 components of the rotational (angular) velocity (in rad/s) of the current tower element relative to the inertial frame origin at ground level [onshore] or MSL [offshore].
   REAL(DbKi),           INTENT(IN   )      :: ZTime                              ! Current simulation time (sec)

   TYPE(Waves_DataType), INTENT(INOUT)      :: WaveDat                            ! wave data
   INTEGER,              INTENT(OUT)        :: ErrStat                            ! A non-zero value indicates an error occurred
   INTEGER,              INTENT(IN )        :: JNode                              ! The number of the current tower node / element (-) [1 to TwrNodes]

   

      ! Local Variables:

   REAL(ReKi)                               :: DZFract                            ! The fraction of the current tower element that is below the free surface of the incident wave and above the seabed (0.0 <= DZFract  <= 1.0): 0.0 = the element is entirely above the free surface, 1.0 = element is entirely below the free surface and above the seabed (-)
   REAL(ReKi)                               :: DZFractS                           ! The fraction of the current tower element that is                                                 above the seabed (0.0 <= DZFractS <= 1.0): 0.0 = the element is entirely below the seabed      , 1.0 = element is entirely                            above the seabed (-)
   REAL(ReKi)                               :: DZFractW                           ! The fraction of the current tower element that is below the free surface of the incident wave                      (0.0 <= DZFractW <= 1.0): 0.0 = the element is entirely above the free surface, 1.0 = element is entirely below the free surface                      (-)
   REAL(ReKi)                               :: InertiaForce     (2)               ! Wave inertia force in the xi- (1) and yi- (2) directions, respectively, on the current tower element at the current time (N)
   REAL(ReKi)                               :: MagVRel                            ! The magnitude of the horizontal incident wave velocity relative to the current tower node at the current time (m/s)
   REAL(ReKi)                               :: MomArm                             ! Moment arm in the vertical direction from the current tower node to the center of pressure of the wave load on the current tower element (meters)
   REAL(ReKi)                               :: TowerAM                            ! Force -translation                     component of TwrAM (kg    /m)
   REAL(ReKi)                               :: TowerAMM                           ! Force -rotation and moment-translation component of TwrAM (kg-m  /m)
   REAL(ReKi)                               :: TowerAMM2                          !                     Moment-rotation    component of TwrAM (kg-m^2/m)
   REAL(ReKi)                               :: TwrArea                            ! Cross-sectional area of current tower element (m^2)
   REAL(ReKi)                               :: TwrVelocity      (2)               ! Velocity of the center of pressure of the wave load on the current tower element in the xi- (1) and yi- (2) directions, respectively, at the current time (m/s)
   REAL(ReKi)                               :: ViscousForce     (2)               ! Viscous drag force in the xi- (1) and yi- (2) directions, respectively, on the current tower element at the current time (N)
   REAL(ReKi)                               :: WaveAcceleration0(2)               ! Acceleration of incident waves in the xi- (1) and yi- (2) directions, respectively, at the current tower node and time (m/s^2)
   REAL(ReKi)                               :: WaveElevation0                     ! Elevation of incident waves at the platform reference point and current time (meters)
   REAL(ReKi)                               :: WaveVelocity0    (2)               ! Velocity     of incident waves in the xi- (1) and yi- (2) directions, respectively, at the current tower node and time (m/s  )

   INTEGER                                  :: K                                  ! Generic index


      ! Initialize the error status
   ErrStat = 0      

      ! Initialize the added mass matrix per unit length of the current tower
      !   element, TwrAM, and the portion of the current tower element load per
      !   unit length associated with everything but the added mass effects,
      !   TwrFt, to zero:

   TwrAM(1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
   TwrAM(2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
   TwrAM(3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
   TwrAM(4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
   TwrAM(5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
   TwrAM(6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)

   TwrFt(1)   = 0.0
   TwrFt(2)   = 0.0
   TwrFt(3)   = 0.0
   TwrFt(4)   = 0.0
   TwrFt(5)   = 0.0
   TwrFt(6)   = 0.0



      ! Find the fraction of the current tower element that is below the free
      !   surface of the incident wave and above the seabed:

!bjj: this could be better if we used temp variables?

   IF ( WaveDat%WaveStMod == 0 )  THEN   ! .TRUE. if we have no stretching; therefore, integrate up to the MSL, regardless of the instantaneous free surface elevation.

      IF (     ( WaveDat%WaveKinzi0(JNode) - 0.5*WaveDat%DZNodes(JNode) ) >= 0.0            )  THEN ! .TRUE. if the current tower element lies entirely above the MSL.
         DZFractW = 0.0
      ELSEIF ( ( WaveDat%WaveKinzi0(JNode) + 0.5*WaveDat%DZNodes(JNode) ) <= 0.0            )  THEN ! .TRUE. if the current tower element lies entirely below the MSL.
         DZFractW = 1.0
      ELSE                                                                                          ! The free surface of the incident wave must fall somewhere along the current tower element; thus, interpolate.
         DZFractW = ( ( 0.0              - ( WaveDat%WaveKinzi0(JNode) - 0.5*WaveDat%DZNodes(JNode) ) )/WaveDat%DZNodes(JNode) )
      ENDIF

   ELSE                          ! We must have some sort of stretching.

      WaveElevation0 = WaveElevation ( 1, ZTime, WaveDat, ErrStat )
      IF ( ErrStat /= 0 ) RETURN

      IF (     ( WaveDat%WaveKinzi0(JNode) - 0.5*WaveDat%DZNodes(JNode) ) >= WaveElevation0 )  THEN ! .TRUE. if the current tower element lies entirely above the free surface of the incident wave.
         DZFractW = 0.0
      ELSEIF ( ( WaveDat%WaveKinzi0(JNode) + 0.5*WaveDat%DZNodes(JNode) ) <= WaveElevation0 )  THEN ! .TRUE. if the current tower element lies entirely below the free surface of the incident wave.
         DZFractW = 1.0
      ELSE                                                                          ! The free surface of the incident wave must fall somewhere along the current tower element; thus, interpolate.
         DZFractW = ( ( WaveElevation0 - ( WaveDat%WaveKinzi0(JNode) - 0.5*WaveDat%DZNodes(JNode) ) )/WaveDat%DZNodes(JNode) )
      ENDIF

   ENDIF

   IF (        ( WaveDat%WaveKinzi0(JNode) - 0.5*WaveDat%DZNodes(JNode) ) >= -WaveDat%WtrDpth )  THEN ! .TRUE. if the current tower element lies entirely above the seabed.
         DZFractS = 1.0
   ELSEIF (    ( WaveDat%WaveKinzi0(JNode) + 0.5*WaveDat%DZNodes(JNode) ) <= -WaveDat%WtrDpth )  THEN ! .TRUE. if the current tower element lies entirely below the seabed.
         DZFractS = 0.0
   ELSE                                                                             ! The seabed must fall somewhere along the current tower element; thus, interpolate.
         DZFractS = ( ( (WaveDat%WaveKinzi0(JNode) + 0.5*WaveDat%DZNodes(JNode)) - ( -WaveDat%WtrDpth ) )/WaveDat%DZNodes(JNode) )
   ENDIF

   DZFract = DZFractW*DZFractS



      ! Compute the hydrodynamic loads using Morison's equation for the portion of
      !   the current tower element that lies below the free surface of the
      !   incident wave and above the seabed:

   IF ( DZFract > 0.0 )  THEN ! .TRUE. if a portion of the current tower element lies below the free surface of the incident wave.


      ! Compute the moment arm in the vertical direction between the current tower
      !   node and the center of pressure of the wave load on the current tower
      !   element:

      MomArm = 0.5*WaveDat%DZNodes(JNode)*( DZFractW - DZFractS )   ! NOTE: MomArm = 0.0 when the entire element is submerged in the fluid; consequently, the roll and pitch components of the load are zero when the entire element is submerged in the fluid


      ! Compute the velocity and acceleration of the incident waves in the xi- (1)
      !   and yi- (2) directions, respectively, at the current tower node and
      !   time:

      DO K = 1,2     ! Loop through the xi- (1) and yi- (2) directions
         WaveVelocity0    (K) = WaveVelocity     ( JNode, K, ZTime, WaveDat, ErrStat )
         IF ( ErrStat /= 0 ) RETURN
         
         WaveAcceleration0(K) = WaveAcceleration ( JNode, K, ZTime, WaveDat, ErrStat )
         IF ( ErrStat /= 0 ) RETURN
      ENDDO          ! K - The xi- (1) and yi- (2) directions


      ! Compute the velocity of the center of pressure of the wave load on the
      !   current tower element in the xi- (1) and yi- (2) directions,
      !   respectively, at the current time:

      TwrVelocity(1) = XD(1) + XD(5)*MomArm
      TwrVelocity(2) = XD(2) - XD(4)*MomArm


      ! Compute the magnitude of the horizontal incident wave velocity relative to
      !   the center of pressure of the wave load on the current tower element at
      !   the current time:

      MagVRel = SQRT(   ( WaveVelocity0(1) - TwrVelocity(1) )**2 &
                      + ( WaveVelocity0(2) - TwrVelocity(2) )**2   )


      ! Compute the cross-sectional area of the current tower element:

      TwrArea = PiOvr4*TwrDiam*TwrDiam


      ! Compute the added mass matrix per unit length of the current tower
      !   element:

      TowerAM    = TwrCA*WaveDat%WtrDens*TwrArea*DZFract    ! force -translation                     component
      TowerAMM   = TowerAM *MomArm                          ! force -rotation and moment-translation component
      TowerAMM2  = TowerAMM*MomArm                          !                     moment-rotation    component

      TwrAM(1,1) = TwrAM(1,1) + TowerAM   ! surge-surge component
      TwrAM(2,2) = TwrAM(2,2) + TowerAM   ! sway -sway  component
      TwrAM(4,4) = TwrAM(4,4) + TowerAMM2 ! roll -roll  component
      TwrAM(5,5) = TwrAM(5,5) + TowerAMM2 ! pitch-pitch component
      TwrAM(2,4) = TwrAM(2,4) - TowerAMM  ! sway -roll  component
      TwrAM(4,2) = TwrAM(4,2) - TowerAMM  ! roll -sway  component
      TwrAM(1,5) = TwrAM(1,5) + TowerAMM  ! surge-pitch component
      TwrAM(5,1) = TwrAM(5,1) + TowerAMM  ! pitch-surge component


      ! Compute the portions of the current tower element load per unit length
      !   associated with the incident wave acceleration and the viscous drag:

      DO K = 1,2     ! Loop through the xi- (1) and yi- (2) directions
         InertiaForce(K) = ( 1.0 + TwrCA )*WaveDat%WtrDens*TwrArea*WaveAcceleration0(K)*DZFract
         ViscousForce(K) = 0.5*TwrCD*WaveDat%WtrDens*TwrDiam*( WaveVelocity0(K) - TwrVelocity(K) )*MagVRel*DZFract
      ENDDO          ! K - The xi- (1) and yi- (2) directions

      TwrFt(1  ) = TwrFt(1  ) +   InertiaForce(1) + ViscousForce(1)           ! surge component
      TwrFt(2  ) = TwrFt(2  ) +   InertiaForce(2) + ViscousForce(2)           ! sway  component
      TwrFt(4  ) = TwrFt(4  ) - ( InertiaForce(2) + ViscousForce(2) )*MomArm  ! roll  component
      TwrFt(5  ) = TwrFt(5  ) + ( InertiaForce(1) + ViscousForce(1) )*MomArm  ! pitch component


   ENDIF



   RETURN
   END SUBROUTINE MorisonTwrLd
!=======================================================================
   SUBROUTINE FB_Terminate( )
   
   ! Deallocate arrays
   
   
   
   ! close files
   
   END SUBROUTINE FB_Terminate

!=======================================================================
END MODULE FixedBottomSupportStructure
!=======================================================================
