!=======================================================================
MODULE FloatingPlatform


   ! This MODULE stores variables and routines used in these time domain
   !   hydrodynamic loading and mooring system dynamics routines for the
   !   floating platform.

   USE                                  NWTC_Library

!BJJ start of proposed change v1.00.00a-bjj
   USE                                  UserLineModule

!REMOVE THE BELOW
!rm REAL(ReKi)                   :: HdroAdMsI (6,6)                                 ! Infinite-frequency limit of the frequency-dependent hydrodynamic added mass matrix from the radiation problem (kg, kg-m, kg-m^2 )
!rm REAL(ReKi)                   :: HdroSttc  (6,6)                                 ! Linear hydrostatic restoring matrix from waterplane area and the center-of-buoyancy (kg/s^2, kg-m/s^2, kg-m^2/s^2)
!rm REAL(ReKi), ALLOCATABLE      :: LAnchHTe  (:)                                   ! Effective horizontal tension at the anchor   of each mooring line (N)
!rm REAL(ReKi), ALLOCATABLE      :: LAnchVTe  (:)                                   ! Effective vertical   tension at the anchor   of each mooring line (N)
!rm REAL(ReKi), ALLOCATABLE      :: LAnchxi   (:)                                   ! xi-coordinate of each anchor   in the inertial frame        coordinate system (meters)
!rm REAL(ReKi), ALLOCATABLE      :: LAnchyi   (:)                                   ! yi-coordinate of each anchor   in the inertial frame        coordinate system (meters)
!rm REAL(ReKi), ALLOCATABLE      :: LAnchzi   (:)                                   ! zi-coordinate of each anchor   in the inertial frame        coordinate system (meters)
!rm REAL(ReKi), ALLOCATABLE      :: LEAStff   (:)                                   ! Extensional stiffness of each mooring line (N)
!rm REAL(ReKi), ALLOCATABLE      :: LFairHTe  (:)                                   ! Effective horizontal tension at the fairlead of each mooring line (N)
!rm REAL(ReKi), ALLOCATABLE      :: LFairVTe  (:)                                   ! Effective vertical   tension at the fairlead of each mooring line (N)
!rm REAL(ReKi), ALLOCATABLE      :: LFairxt   (:)                                   ! xt-coordinate of each fairlead in the tower base / platform coordinate system (meters)
!rm REAL(ReKi), ALLOCATABLE      :: LFairyt   (:)                                   ! yt-coordinate of each fairlead in the tower base / platform coordinate system (meters)
!rm REAL(ReKi), ALLOCATABLE      :: LFairzt   (:)                                   ! zt-coordinate of each fairlead in the tower base / platform coordinate system (meters)
!rm REAL(ReKi), ALLOCATABLE      :: LFldWght  (:)                                   ! Weight of each mooring line in fluid per unit length (N/m)
!rm REAL(ReKi), ALLOCATABLE      :: LNodesPi  (:,:,:)                               ! xi- (1), yi- (2), and zi (3) -coordinates in the inertial frame of the position of each node of each line where the line position and tension can be output (meters)
!rm REAL(ReKi), ALLOCATABLE      :: LNodesTe  (:,:)                                 ! Effective line tensions                                                         at each node of each line where the line position and tension can be output (N     )

!rm REAL(ReKi), ALLOCATABLE      :: LNodesX   (:)                                   ! X -coordinates in the local coordinate system of the current line (this coordinate system lies at the current anchor, Z being vertical, and X directed from the current anchor to the current fairlead) of each node where the line position and tension can be output (meters)
!rm REAL(ReKi), ALLOCATABLE      :: LNodesZ   (:)                                   ! Z -coordinates in the local coordinate system of the current line (this coordinate system lies at the current anchor, Z being vertical, and X directed from the current anchor to the current fairlead) of each node where the line position and tension can be output (meters)

!rm REAL(ReKi), ALLOCATABLE      :: LSeabedCD (:)                                   ! Coefficient of seabed static friction drag of each mooring line (a negative value indicates no seabed) (-)
!rm REAL(ReKi), ALLOCATABLE      :: LSNodes   (:,:)                                 ! Unstretched arc distance along mooring line from anchor to each node where the line position and tension can be output (meters)
!rm REAL(ReKi), ALLOCATABLE      :: LTenTol   (:)                                   ! Convergence tolerance within Newton-Raphson iteration of each mooring line specified as a fraction of tension (-)
!rm REAL(ReKi), ALLOCATABLE      :: LUnstrLen (:)                                   ! Unstretched length of each mooring line (meters)
!rm REAL(ReKi)                   :: PtfmCD                                          ! Effective platform normalized hydrodynamic viscous drag coefficient in calculation of viscous drag term from Morison's equation (-)
!rm REAL(ReKi)                   :: PtfmDiam                                        ! Effective platform diameter in calculation of viscous drag term from Morison's equation (meters)
!rm REAL(ReKi)                   :: PtfmVol0                                        ! Displaced volume of water when the platform is in its undisplaced position (m^3)
!rm REAL(ReKi)                   :: RdtnDT                                          ! Time step for wave radiation kernel calculations (sec)
!rm REAL(ReKi)                   :: RdtnTMax                                        ! Analysis time for wave radiation kernel calculations (sec)
!rm REAL(ReKi), ALLOCATABLE      :: RdtnKrnl  (:,:,:)                               ! Instantaneous values of the wave radiation kernel (kg/s^2, kg-m/s^2, kg-m^2/s^2)
!rm REAL(ReKi), ALLOCATABLE      :: WaveExctn (:,:)                                 ! Instantaneous values of the total excitation force on the support platfrom from incident waves (N, N-m)
!rm REAL(ReKi), ALLOCATABLE      :: XDHistory (:,:)                                 ! The time history of the 3 components of the translational velocity        (in m/s)        of the platform reference and the 3 components of the rotational (angular) velocity  (in rad/s)        of the platform relative to the inertial frame

!rm INTEGER                      :: LineMod                                         ! Mooring line model switch {0: none, 1: standard quasi-static, 2: user-defined from routine UserLine} (switch)
!rm INTEGER                      :: LineNodes                                       ! Number of nodes per line where the mooring line position and tension can be output (-)
!rm INTEGER                      :: NStepRdtn                                       ! Total number of frequency components = total number of time steps in the wave radiation kernel (-)
!rm INTEGER                      :: NStepRdtn1                                      ! = NStepRdtn + 1 (-)
!rm INTEGER                      :: NumLines                                        ! Number of mooring lines (-)
!rm
!rm LOGICAL                      :: UseRdtn                                         ! Flag for determining whether or not the to model wave radiation damping (flag)
!rm
!rm CHARACTER(1024)              :: DirRoot
!REMOVE THE ABOVE

TYPE, PUBLIC :: Line_InitDataType
   REAL(ReKi)                           :: LAnchxi                                ! xi-coordinate of the anchor   in the inertial frame        coordinate system (meters)
   REAL(ReKi)                           :: LAnchyi                                ! yi-coordinate of the anchor   in the inertial frame        coordinate system (meters)
   REAL(ReKi)                           :: LAnchzi                                ! zi-coordinate of the anchor   in the inertial frame        coordinate system (meters)
   REAL(ReKi)                           :: LDiam                                  ! Effective diameter of this mooring line for calculation of the line buoyancy (meters)
   REAL(ReKi)                           :: LEAStff                                ! Extensional stiffness of each mooring line (N)
   REAL(ReKi)                           :: LFairxt                                ! xt-coordinate of fairlead in the tower base / platform coordinate system (meters)
   REAL(ReKi)                           :: LFairyt                                ! yt-coordinate of fairlead in the tower base / platform coordinate system (meters)
   REAL(ReKi)                           :: LFairzt                                ! zt-coordinate of fairlead in the tower base / platform coordinate system (meters)
   REAL(ReKi)                           :: LMassDen                               ! Mass density of this mooring line (kg/m)
   REAL(ReKi)                           :: LSeabedCD                              ! Coefficient of seabed static friction drag of this mooring line (a negative value indicates no seabed) (-)
   REAL(ReKi), ALLOCATABLE              :: LSNodes    (:)                         ! Unstretched arc distance along mooring line from anchor to each node where the line position and tension can be output (meters)
   REAL(ReKi)                           :: LTenTol                                ! Convergence tolerance within Newton-Raphson iteration of this mooring line specified as a fraction of tension (-)
   REAL(ReKi)                           :: LUnstrLen                              ! Unstretched length of this mooring line (meters)

   INTEGER                              :: LineNodes              = 0             ! Number of nodes in this line where the mooring line position and tension can be output (-)
END TYPE Line_InitDataType


TYPE, PUBLIC :: FltPtfm_InitDataType
   REAL(ReKi)                           :: PtfmCD                 = 0.0           ! Effective platform normalized hydrodynamic viscous drag coefficient in calculation of viscous drag term from Morison's equation (-)
   REAL(ReKi)                           :: PtfmDiam               = 0.0           ! Effective platform diameter in calculation of viscous drag term from Morison's equation (meters)
   REAL(ReKi)                           :: PtfmVol0               = 0.0           ! Displaced volume of water when the platform is in its undisplaced position (m^3)
   REAL(ReKi)                           :: RdtnDT                 = 0.0           ! Time step for wave radiation kernel calculations (sec)
   REAL(ReKi)                           :: RdtnTMax               = 0.0           ! Analysis time for wave radiation kernel calculations; the actual analysis time may be larger than this value in order for the maintain an effecient (co)sine transform (sec)
   REAL(ReKi)                           :: WAMITULEN              = 1.0           ! Characteristic body length scale used to redimensionalize WAMIT output (meters)
!bjj rm v1.00.00a-bjj:   REAL(ReKi)                           :: X0         (6)                         ! The 3 components of the initial translational displacement (in m) of the platform reference and the 3 components of the initial rotational displacement (in rad) of the platform relative to the inertial frame

   TYPE(Line_InitDataType), ALLOCATABLE :: MooringLine(:)
   INTEGER                              :: LineMod                = 0             ! Mooring line model switch {0: none, 1: standard quasi-static, 2: user-defined from routine UserLine} (switch)
   INTEGER                              :: NumLines               = 0             ! Number of mooring lines (-)

   CHARACTER(1024)                      :: DirRoot                = ""
   CHARACTER(1024)                      :: WAMITFile              = ""            ! Root name of WAMIT output files containing the linear, nondimensionalized, hydrostatic restoring matrix (.hst extension), frequency-dependent hydrodynamic added mass matrix and damping matrix (.1 extension), and frequency- and direction-dependent wave excitation force vector per unit wave amplitude (.3 extension).
END TYPE FltPtfm_InitDataType


TYPE, PRIVATE :: Line_DataType
   REAL(ReKi)                           :: LAnchxi                                ! xi-coordinate of this anchor   in the inertial frame        coordinate system (meters)
   REAL(ReKi)                           :: LAnchyi                                ! yi-coordinate of this anchor   in the inertial frame        coordinate system (meters)
   REAL(ReKi)                           :: LAnchzi                                ! zi-coordinate of this anchor   in the inertial frame        coordinate system (meters)
   REAL(ReKi)                           :: LEAStff                                ! Extensional stiffness of this mooring line (N)
   REAL(ReKi)                           :: LFairxt                                ! xt-coordinate of this fairlead in the tower base / platform coordinate system (meters)
   REAL(ReKi)                           :: LFairyt                                ! yt-coordinate of this fairlead in the tower base / platform coordinate system (meters)
   REAL(ReKi)                           :: LFairzt                                ! zt-coordinate of this fairlead in the tower base / platform coordinate system (meters)
   REAL(ReKi)                           :: LFldWght                               ! Weight of this mooring line in fluid per unit length (N/m)
   REAL(ReKi), ALLOCATABLE              :: LNodesPi  (:,:)                        ! xi- (1), yi- (2), and zi (3) -coordinates in the inertial frame of the position of each node of this line where the line position and tension can be output (meters)
   REAL(ReKi)                           :: LSeabedCD                              ! Coefficient of seabed static friction drag of this mooring line (a negative value indicates no seabed) (-)
   REAL(ReKi), ALLOCATABLE              :: LSNodes   (:)                          ! Unstretched arc distance along mooring line from anchor to each node where the line position and tension can be output (meters)
   REAL(ReKi)                           :: LTenTol                                ! Convergence tolerance within Newton-Raphson iteration of this mooring line specified as a fraction of tension (-)
   REAL(ReKi)                           :: LUnstrLen                              ! Unstretched length of this mooring line (meters)

   REAL(ReKi)                           :: LAnchHTe               = 0.0           ! Effective horizontal tension at the anchor   of this mooring line (N)
   REAL(ReKi)                           :: LAnchVTe               = 0.0           ! Effective vertical   tension at the anchor   of this mooring line (N)
   REAL(ReKi)                           :: LFairHTe               = 0.0           ! Effective horizontal tension at the fairlead of this mooring line (N)
   REAL(ReKi)                           :: LFairVTe               = 0.0           ! Effective vertical   tension at the fairlead of this mooring line (N)
   REAL(ReKi), ALLOCATABLE              :: LNodesTe  (:)                          ! Effective line tensions                                                         at each node of this line where the line position and tension can be output (N     )
   REAL(ReKi), ALLOCATABLE              :: LNodesX   (:)                          ! X -coordinates in the local coordinate system of the current line (this coordinate system lies at the current anchor, Z being vertical, and X directed from the current anchor to the current fairlead) of each node where the line position and tension can be output (meters)
   REAL(ReKi), ALLOCATABLE              :: LNodesZ   (:)                          ! Z -coordinates in the local coordinate system of the current line (this coordinate system lies at the current anchor, Z being vertical, and X directed from the current anchor to the current fairlead) of each node where the line position and tension can be output (meters)

   INTEGER                              :: LineNodes              = 0             ! Number of nodes in this line where the mooring line position and tension can be output (-)
END TYPE Line_DataType

TYPE, PUBLIC :: FltPtfm_DataType
   PRIVATE

   REAL(ReKi)                           :: HdroAdMsI (6,6)                        ! Infinite-frequency limit of the frequency-dependent hydrodynamic added mass matrix from the radiation problem (kg, kg-m, kg-m^2 )
   REAL(ReKi)                           :: HdroSttc  (6,6)                        ! Linear hydrostatic restoring matrix from waterplane area and the center-of-buoyancy (kg/s^2, kg-m/s^2, kg-m^2/s^2)
   REAL(ReKi)                           :: PtfmCD                                 ! Effective platform normalized hydrodynamic viscous drag coefficient in calculation of viscous drag term from Morison's equation (-)
   REAL(ReKi)                           :: PtfmDiam                               ! Effective platform diameter in calculation of viscous drag term from Morison's equation (meters)
   REAL(ReKi)                           :: PtfmVol0                               ! Displaced volume of water when the platform is in its undisplaced position (m^3)
   REAL(ReKi)                           :: RdtnDT                                 ! Time step for wave radiation kernel calculations (sec)
   REAL(ReKi)                           :: RdtnTMax                               ! Analysis time for wave radiation kernel calculations (sec)
   REAL(ReKi), ALLOCATABLE              :: RdtnKrnl  (:,:,:)                      ! Instantaneous values of the wave radiation kernel (kg/s^2, kg-m/s^2, kg-m^2/s^2)
   REAL(ReKi), ALLOCATABLE              :: WaveExctn (:,:)                        ! Instantaneous values of the total excitation force on the support platfrom from incident waves (N, N-m)
   REAL(ReKi), ALLOCATABLE              :: XDHistory (:,:)                        ! The time history of the 3 components of the translational velocity        (in m/s)        of the platform reference and the 3 components of the rotational (angular) velocity  (in rad/s)        of the platform relative to the inertial frame

   REAL(ReKi)                           :: LastTime               = 0.0           ! Last time the values in XDHistory were saved (sec)

   TYPE(Line_DataType), ALLOCATABLE     :: MooringLine (:)
   TYPE(UserLine_DataType)              :: UserLine_Data                          ! Stores data for the UserLineModule
   
   INTEGER                              :: LineMod                                ! Mooring line model switch {0: none, 1: standard quasi-static, 2: user-defined from routine UserLine} (switch)
   INTEGER                              :: NStepRdtn                              ! Total number of frequency components = total number of time steps in the wave radiation kernel (-)
   INTEGER                              :: NStepRdtn1                             ! = NStepRdtn + 1 (-)
   INTEGER                              :: NumLines               = 0             ! Number of mooring lines (-)

   INTEGER                              :: LastIndRdtn                            ! Index into the radiation     arrays saved from the last call as a starting point for current call.
   INTEGER                              :: LastIndRdtn2                           ! Index into the radiation     arrays saved from the last call as a starting point for current call.
   INTEGER                              :: LastIndWave            = 1             ! Index into the incident wave arrays saved from the last call as a starting point for current call.

   LOGICAL                              :: CalculateFirstGuess    = .TRUE.        ! Flag for determining if we're on the first time step and therefore need an initial quess for the Newton Raphson (can't use previous time step)
   LOGICAL                              :: UseRdtn                = .FALSE.       ! Flag for determining whether or not the to model wave radiation damping (flag)   
END TYPE FltPtfm_DataType

!TYPE, PUBLIC :: UserLine_DataType
!   PRIVATE
!   REAL(ReKi), ALLOCATABLE      :: LAnchHTe  (:)                                   ! Effective horizontal tension at the anchor   of each mooring line (N)
!   REAL(ReKi), ALLOCATABLE      :: LAnchVTe  (:)                                   ! Effective vertical   tension at the anchor   of each mooring line (N)
!   REAL(ReKi), ALLOCATABLE      :: LFairHTe  (:)                                   ! Effective horizontal tension at the fairlead of each mooring line (N)
!   REAL(ReKi), ALLOCATABLE      :: LFairVTe  (:)                                   ! Effective vertical   tension at the fairlead of each mooring line (N)
!   INTEGER                      :: LineNodes                                       ! Number of nodes per line where the mooring line position and tension can be output (-)
!   REAL(ReKi), ALLOCATABLE      :: LNodesTe  (:,:)                                 ! Effective line tensions                                                         at each node of each line where the line position and tension can be output (N     )
!   REAL(ReKi), ALLOCATABLE      :: LNodesX   (:)                                   ! X -coordinates in the local coordinate system of the current line (this coordinate system lies at the current anchor, Z being vertical, and X directed from the current anchor to the current fairlead) of each node where the line position and tension can be output (meters)
!   REAL(ReKi), ALLOCATABLE      :: LNodesZ   (:)                                   ! Z -coordinates in the local coordinate system of the current line (this coordinate system lies at the current anchor, Z being vertical, and X directed from the current anchor to the current fairlead) of each node where the line position and tension can be output (meters)
!END TYPE UserLine_DataType
!BJJ end of proposed change v1.00.00a-bjj

   REAL(ReKi), PARAMETER, PRIVATE      :: OnePlusEps  = 1.0 + EPSILON(OnePlusEps)   ! The number slighty greater than unity in the precision of ReKi.

CONTAINS
!=======================================================================
   SUBROUTINE AnchorTension ( ILine, AnchTe, AnchTeAng, FP_Data, ErrStat )


      ! This SUBROUTINE is used to return the instantaneous effective
      ! tension in, and the vertical angle of, the line at the anchor for
      ! mooring line ILine to the calling program.


   IMPLICIT                            NONE


      ! Passed Variables:

   REAL(ReKi), INTENT(OUT)            :: AnchTe                                   ! Instantaneous effective tension in the line at the anchor (N  )
   REAL(ReKi), INTENT(OUT)            :: AnchTeAng                                ! Instantaneous vertical angle    of the line at the anchor (rad)

   INTEGER,    INTENT(IN )            :: ILine                                    ! Mooring line number (-)

!bjj start of proposed change v1.00.00a-bjj
   TYPE(FltPtfm_DataType), INTENT(IN) :: FP_Data
   INTEGER,    INTENT(OUT)            :: ErrStat                                  ! a non-zero value indicates an error has occurred

   ErrStat = 0
!bjj end of proposed change


      ! Abort if the mooring line parameters have not been computed yet or if
      !   ILine is not one of the existing mooring lines:

!bjj start of proposed change v1.00.00a-bjj
!this must be reorganized with the new data types
!   IF ( .NOT. ALLOCATED ( LAnchHTe )                )  THEN
!      CALL ProgAbort ( ' Routine InitFltngPtfmLd() must be called before routine AnchorTension().', TrapErrors = .TRUE.)
!      ErrStat = 1
!      RETURN
!   ELSEIF ( ( ILine < 1 ) .OR. ( ILine > FP_Data%NumLines ) )  THEN
!      CALL ProgAbort ( ' Mooring line '//TRIM( Int2LStr( ILine ) )//' has not been analyzed.', TrapErrors = .TRUE.)
!      ErrStat = 1
!      RETURN
!   ENDIF


   IF ( ( ILine < 1 ) .OR. ( ILine > FP_Data%NumLines ) )  THEN
      CALL ProgAbort ( ' Mooring line '//TRIM( Int2LStr( ILine ) )//' has not been analyzed.', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( .NOT. ALLOCATED ( FP_Data%MooringLine )     )  THEN
      CALL ProgAbort ( ' Routine InitFltngPtfmLd() must be called before routine AnchorTension().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ENDIF
!bjj end of proposed change v1.00.00a-bjj


      ! Return the instantaneous effective tension and angle:

   AnchTe       = SQRT(  FP_Data%MooringLine(ILine)%LAnchHTe**2 + FP_Data%MooringLine(ILine)%LAnchVTe**2 )
   IF ( AnchTe == 0.0 )  THEN ! .TRUE. if the effective tension at the anchor is zero so that ATAN2() will be ill-conditioned; return zero instead
      AnchTeAng = 0.0
   ELSE
      AnchTeAng = ATAN2( FP_Data%MooringLine(ILine)%LAnchVTe    , FP_Data%MooringLine(ILine)%LAnchHTe    )
   ENDIF



   RETURN
   END SUBROUTINE AnchorTension
!=======================================================================
   SUBROUTINE FairleadTension ( ILine, FairTe, FairTeAng, FP_Data, ErrStat )


      ! This SUBROUTINE is used to return the instantaneous effective
      ! tension in, and the vertical angle of, the line at the fairlead for
      ! mooring line ILine to the calling program.


   IMPLICIT                            NONE


      ! Passed Variables:

   REAL(ReKi), INTENT(OUT)            :: FairTe                                   ! Instantaneous effective tension in the line at the fairlead (N  )
   REAL(ReKi), INTENT(OUT)            :: FairTeAng                                ! Instantaneous vertical angle    of the line at the fairlead (rad)

   INTEGER,    INTENT(IN )            :: ILine                                    ! Mooring line number (-)

!bjj start of proposed change v1.00.00a-bjj
   TYPE(FltPtfm_DataType), INTENT(IN) :: FP_Data
   INTEGER,    INTENT(OUT)            :: ErrStat                                  ! a non-zero value indicates an error has occurred

   ErrStat = 0
!bjj end of proposed change


      ! Abort if the mooring line parameters have not been computed yet or if
      !   ILine is not one of the existing mooring lines:

!bjj start of proposed change v1.00.00a-bjj
!bjj: this must be reordered for the new data types:
!   IF ( .NOT. ALLOCATED ( LFairHTe )                )  THEN
!      CALL ProgAbort ( ' Routine InitFltngPtfmLd() must be called before routine FairleadTension().', TrapErrors = .TRUE.)
!      ErrStat = 1
!      RETURN
!   ELSEIF ( ( ILine < 1 ) .OR. ( ILine > FP_Data%NumLines ) )  THEN
!      CALL ProgAbort ( ' Mooring line '//TRIM( Int2LStr( ILine ) )//' has not been analyzed.', TrapErrors = .TRUE.)
!      ErrStat = 1
!      RETURN
!   ENDIF

   IF ( ( ILine < 1 ) .OR. ( ILine > FP_Data%NumLines ) )  THEN
      CALL ProgAbort ( ' Mooring line '//TRIM( Int2LStr( ILine ) )//' has not been analyzed.', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( .NOT. ALLOCATED ( FP_Data%MooringLine )     )  THEN
      CALL ProgAbort ( ' Routine InitFltngPtfmLd() must be called before routine FairleadTension().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ENDIF

!bjj end of proposed change v1.00.00a-bjj


      ! Return the instantaneous effective tension and angle:

   FairTe       = SQRT(  FP_Data%MooringLine(ILine)%LFairHTe**2 + FP_Data%MooringLine(ILine)%LFairVTe**2 )
   IF ( FairTe == 0.0 )  THEN ! .TRUE. if the effective tension at the fairlead is zero so that ATAN2() will be ill-conditioned; return zero instead
      FairTeAng = 0.0
   ELSE
      FairTeAng = ATAN2( FP_Data%MooringLine(ILine)%LFairVTe    , FP_Data%MooringLine(ILine)%LFairHTe    )
   ENDIF



   RETURN
   END SUBROUTINE FairleadTension
!=======================================================================
   SUBROUTINE FltngPtfmLd ( X, XD, ZTime, PtfmAM, PtfmFt, WaveDat, FP_Data, ErrStat )


      ! This routine implements the time domain hydrodynamic loading and
      ! mooring system dynamics equations for the floating platform.  The
      ! loads from the mooring system are obtained through an interface
      ! with the LINES module of SML and include contributions from
      ! inertia, restoring, and viscous separation damping effects,
      ! including the elastic response of multi-segment lines and their
      ! interaction with the seabed.  The hydrodynamic loads on the support
      ! platform include the restoring contributions of buoyancy and
      ! waterplane area from hydrostatics; the viscous drag contributions
      ! from Morison's equation; the added mass and damping contributions
      ! from radiation, including free surface memory effects; and the
      ! incident wave excitation from diffraction in regular or irregular
      ! seas.  The diffraction, radiation, and hydrostatics problems are
      ! implemented in their true linear form.  Not included in the model
      ! are the effects of vortex-induced vibration and loading from sea
      ! ice, as well as the nonlinear effects of slow-drift and sum-
      ! frequency excitation and high-order wave kinematics.

      ! The order of indices in all arrays passed to and from this routine
      !   is as follows:
      !      1 = Platform surge / xi-component of platform translation (internal DOF index = DOF_Sg)
      !      3 = Platform sway  / yi-component of platform translation (internal DOF index = DOF_Sw)
      !      3 = Platform heave / zi-component of platform translation (internal DOF index = DOF_Hv)
      !      4 = Platform roll  / xi-component of platform rotation    (internal DOF index = DOF_R )
      !      5 = Platform pitch / yi-component of platform rotation    (internal DOF index = DOF_P )
      !      6 = Platform yaw   / zi-component of platform rotation    (internal DOF index = DOF_Y )

      ! NOTE: In extreme wave conditions (i.e., hurricane sea states,
      !       etc.), radiation damping will not be important but Morison's
      !       equation with stretching is very important.  Therefore, these
      !       hydrodynamic loading equations, which are in true linear
      !       form, are limited in their accuracy in such sea conditions.


   USE                                  Waves  !RhoXg


   IMPLICIT                             NONE


      ! Passed Variables:

   REAL(ReKi),            INTENT(  OUT) :: PtfmAM   (6,6)                         ! Platform added mass matrix (kg, kg-m, kg-m^2)
   REAL(ReKi),            INTENT(  OUT) :: PtfmFt   (6)                           ! The 3 components of the portion of the platform force (in N  ) acting at the platform reference and the 3 components of the portion of the platform moment (in N-m  ) acting at the platform reference associated with everything but the added-mass effects; positive forces are in the direction of motion
   REAL(ReKi),            INTENT(IN   ) :: X        (6)                           ! The 3 components of the translational displacement    (in m  )        of the platform reference and the 3 components of the rotational displacement        (in rad  )        of the platform relative to the inertial frame
   REAL(ReKi),            INTENT(IN   ) :: XD       (6)                           ! The 3 components of the translational velocity        (in m/s)        of the platform reference and the 3 components of the rotational (angular) velocity  (in rad/s)        of the platform relative to the inertial frame
   REAL(ReKi),            INTENT(IN   ) :: ZTime                                  ! Current simulation time (sec)

!bjj start of proposed change v1.00.00a-bjj
   TYPE(Waves_DataType)  ,INTENT(INOUT) :: WaveDat                                ! the data structure containing the wave data
   TYPE(FltPtfm_DataType),INTENT(INOUT) :: FP_Data                                ! data structure containing the floating platform data
   INTEGER,               INTENT(  OUT) :: ErrStat                                ! A non-zero value indicates an error occurred
!bjj end

      ! Local Variables:

   REAL(ReKi)                           :: COSPhi                                  ! Cosine of the angle between the xi-axis of the inertia frame and the X-axis of the local coordinate system of the current mooring line (-)
   REAL(ReKi)                           :: F_HS     (6)                            ! Total load contribution from hydrostatics, including the effects of waterplane area and the center of buoyancy (N, N-m)
   REAL(ReKi)                           :: F_Lines  (6)                            ! Total load contribution from all mooring lines (N, N-m)
   REAL(ReKi)                           :: F_Rdtn   (6)                            ! Total load contribution from wave radiation damping (i.e., the diffraction problem) (N, N-m)
   REAL(ReKi)                           :: F_RdtnDT (6)                            ! The portion of the total load contribution from wave radiation damping associated with the convolution integral proportional to ( RdtnDT - RdtnRmndr ) (N, N-m)
   REAL(ReKi)                           :: F_RdtnRmndr(6)                          ! The portion of the total load contribution from wave radiation damping associated with the convolution integral proportional to (          RdtnRmndr ) (N, N-m)
   REAL(ReKi)                           :: F_Viscous(6)                            ! Total load contribution from viscous drag (N, N-m)
   REAL(ReKi)                           :: F_Waves  (6)                            ! Total load contribution from incident waves (i.e., the diffraction problem) (N, N-m)
   REAL(ReKi)                           :: IncrmntXD                               ! Incremental change in XD over a single radiation time step (m/s, rad/s)
!bjj start of proposed change   
!this must be stored in the data type so we can call multiple instances
!rm   REAL(ReKi), SAVE                     :: LastTime    = 0.0                       ! Last time the values in XDHistory where saved (sec)
!bjj end of proposed change   
!bjj start of proposed change v1.00.00a-bjj
   REAL(ReKi)                           :: Lamda0                                  ! Catenary parameter used to generate the initial guesses of the horizontal and vertical tensions at the fairlead for the Newton-Raphson iteration (-)
!bjj end of proposed change v1.00.00a-bjj   
   REAL(ReKi)                           :: LFairxi                                 ! xi-coordinate of the current fairlead in the inertial frame coordinate system (meters)
   REAL(ReKi)                           :: LFairxiRef                              ! xi-coordinate of the current fairlead relative to the platform reference point in the inertial frame coordinate system (meters)
   REAL(ReKi)                           :: LFairxiTe                               ! xi-component of the effective tension at the fairlead of the current mooring line (N)
   REAL(ReKi)                           :: LFairyi                                 ! yi-coordinate of the current fairlead in the inertial frame coordinate system (meters)
   REAL(ReKi)                           :: LFairyiRef                              ! yi-coordinate of the current fairlead relative to the platform reference point in the inertial frame coordinate system (meters)
   REAL(ReKi)                           :: LFairyiTe                               ! yi-component of the effective tension at the fairlead of the current mooring line (N)
   REAL(ReKi)                           :: LFairzi                                 ! zi-coordinate of the current fairlead in the inertial frame coordinate system (meters)
   REAL(ReKi)                           :: LFairziRef                              ! zi-coordinate of the current fairlead relative to the platform reference point in the inertial frame coordinate system (meters)
   REAL(ReKi)                           :: LFairziTe                               ! zi-component of the effective tension at the fairlead of the current mooring line (N)
   REAL(ReKi)                           :: MagVRel                                 ! The magnitude of the horizontal incident wave velocity relative to the current platform node at the current time (m/s)
!bjj START OF PROPOSED CHANGE V1.00.00A-BJJ: this parameter can be in the module data
!RM   REAL(ReKi), PARAMETER                :: OnePlusEps  = 1.0 + EPSILON(OnePlusEps) ! The number slighty greater than unity in the precision of ReKi.
!bjj END OF PROPOSED CHANGE: this parameter can be in the module data
   REAL(ReKi)                           :: PtfmVelocity     (2)                    ! Velocity of the current platform node in the xi- (1) and yi- (2) directions, respectively, at the current time (m/s)
   REAL(ReKi)                           :: RdtnRmndr                               ! ZTime - RdtnTime(IndRdtn)
   REAL(ReKi)                           :: SINPhi                                  ! Sine   of the angle between the xi-axis of the inertia frame and the X-axis of the local coordinate system of the current mooring line (-)
   REAL(ReKi)                           :: TransMat (3,3)                          ! Transformation matrix from the inertial frame to the tower base / platform coordinate system (-)
   REAL(ReKi)                           :: ViscousForce     (2)                    ! Viscous drag force in the xi- (1) and yi- (2) directions, respectively, on the current platform element at the current time (N)
   REAL(ReKi)                           :: WaveVelocity0    (2)                    ! Velocity of incident waves in the xi- (1) and yi- (2) directions, respectively, at the current platform node and time (m/s)
   REAL(ReKi)                           :: XF                                      ! Horizontal distance between anchor and fairlead of the current mooring line (meters)
   REAL(ReKi)                           :: XF2                                     ! = XF*XF
   REAL(ReKi)                           :: ZF                                      ! Vertical   distance between anchor and fairlead of the current mooring line (meters)
   REAL(ReKi)                           :: ZF2                                     ! = ZF*ZF



   INTEGER                              :: I                                       ! Generic index
   INTEGER                              :: IndRdtn                                 ! Generic index for the radiation problem
   INTEGER                              :: J                                       ! Generic index
   INTEGER                              :: JNode                                   ! The index of the current platform node / element (-) [1 to PtfmNodes]
   INTEGER                              :: K                                       ! Generic index
!bjj start of proposed change   
!this must be stored in the data type so we can call multiple instances
!rm   INTEGER,    SAVE                     :: LastIndRdtn                             ! Index into the radiation     arrays saved from the last call as a starting point for this call.
!rm   INTEGER,    SAVE                     :: LastIndRdtn2                            ! Index into the radiation     arrays saved from the last call as a starting point for this call.
!rm   INTEGER,    SAVE                     :: LastIndWave = 1                         ! Index into the incident wave arrays saved from the last call as a starting point for this call.
!bjj end of proposed change

!bjj start of proposed change v1.00.00a-bjj
      ! initialize the error status
   ErrStat = 0

!bjj end of proposed change

      ! Abort if the wave excitation loads have not been computed yet:

   IF ( .NOT. ALLOCATED ( FP_Data%WaveExctn ) )  THEN
      CALL ProgAbort ( ' Routine InitFltngPtfmLd() must be called before routine FltngPtfmLd().', TrapErrors = .TRUE. )
      ErrStat = 1
      RETURN
   END IF

      ! Compute the load contribution from incident waves (i.e., the diffraction
      !   problem):

   DO I = 1,6     ! Loop through all wave excitation forces and moments
      F_Waves(I) = InterpStp ( ZTime, WaveDat%WaveTime(:), FP_Data%WaveExctn(:,I), FP_Data%LastIndWave, WaveDat%NStepWave )
   ENDDO          ! I - All wave excitation forces and moments


!bjj start of proposed change v1.00.00a-bjj
   IF ( FP_Data%CalculateFirstGuess ) THEN
   
      ! Get the transformation matrix, TransMat0, from the inertial frame to the
      !   initial tower base / platform coordinate system:

      CALL SmllRotTrans ( 'platform displacement', X(4), X(5), X(6), TransMat )


      DO I = 1,FP_Data%NumLines ! Loop through all mooring lines


      ! Transform the fairlead location from the initial platform to the inertial
      !    frame coordinate system:
      ! NOTE: TransMat0^T = TransMat0^-1 where ^T = matrix transpose and ^-1 =
      !       matrix inverse.

         LFairxi = X(1) + TransMat(1,1)*FP_Data%MooringLine(I)%LFairxt &
                        + TransMat(2,1)*FP_Data%MooringLine(I)%LFairyt &
                        + TransMat(3,1)*FP_Data%MooringLine(I)%LFairzt
         LFairyi = X(2) + TransMat(1,2)*FP_Data%MooringLine(I)%LFairxt &
                        + TransMat(2,2)*FP_Data%MooringLine(I)%LFairyt &
                        + TransMat(3,2)*FP_Data%MooringLine(I)%LFairzt
         LFairzi = X(3) + TransMat(1,3)*FP_Data%MooringLine(I)%LFairxt &
                        + TransMat(2,3)*FP_Data%MooringLine(I)%LFairyt &
                        + TransMat(3,3)*FP_Data%MooringLine(I)%LFairzt


      ! Transform the fairlead location from the inertial frame coordinate system
      !   to the local coordinate system of the current line (this coordinate
      !   system lies at the current anchor, Z being vertical, and X directed from
      !   the current anchor to the current fairlead):

         XF      = SQRT( ( LFairxi - FP_Data%MooringLine(I)%LAnchxi )**2 + ( LFairyi - FP_Data%MooringLine(I)%LAnchyi )**2 )
         ZF      =         LFairzi - FP_Data%MooringLine(I)%LAnchzi

         XF2     = XF*XF
         ZF2     = ZF*ZF


      ! Generate the initial guess values for the horizontal and vertical tensions
      !   at the fairlead in the Newton-Raphson iteration for the catenary mooring
      !   line solution.  Use starting values documented in: Peyrot, Alain H. and
      !   Goulois, A. M., "Analysis Of Cable Structures," Computers & Structures,
      !   Vol. 10, 1979, pp. 805-813:

         IF     ( XF                               == 0.0               )  THEN ! .TRUE. if the current mooring line is exactly vertical
            Lamda0 = 1.0E+06
         ELSEIF ( FP_Data%MooringLine(I)%LUnstrLen <= SQRT( XF2 + ZF2 ) )  THEN ! .TRUE. if the current mooring line is taut
            Lamda0 = 0.2
         ELSE                                                                   ! The current mooring line must be slack and not vertical
            Lamda0 = SQRT( 3.0*( ( FP_Data%MooringLine(I)%LUnstrLen**2 - ZF2 )/XF2 - 1.0 ) )
         ENDIF

         FP_Data%MooringLine(I)%LFairHTe = ABS( 0.5*FP_Data%MooringLine(I)%LFldWght*  XF/     Lamda0    )
         FP_Data%MooringLine(I)%LFairVTe =      0.5*FP_Data%MooringLine(I)%LFldWght*( ZF/TANH(Lamda0) + &
                                                    FP_Data%MooringLine(I)%LUnstrLen                    )


      ENDDO             ! I - All mooring lines
   
   
      FP_Data%CalculateFirstGuess = .FALSE.
   END IF
!bjj end of proposed change v1.00.00a-bjj

      ! Compute the load contribution from hydrostatics:

   F_HS(:) = 0.0                      ! Initialize to zero...
   F_HS(3) = WaveDat%RhoXg*FP_Data%PtfmVol0   ! except for the hydrostatic buoyancy force from Archimede's Principle when the support platform is in its undisplaced position
   DO I = 1,6     ! Loop through all hydrostatic forces and moments
      DO J = 1,6  ! Loop through all platform DOFs
         F_HS(I) = F_HS(I) - FP_Data%HdroSttc(I,J)*X(J)
      ENDDO       ! J - All platform DOFs
   ENDDO          ! I - All hydrostatic forces and moments


      ! If necessary, compute the load contribution from wave radiation damping
      !   (i.e., the radiation problem):

   IF ( FP_Data%UseRdtn )  THEN ! .TRUE. when we will be modeling wave radiation damping.


      ! Find the index IndRdtn, where RdtnTime(IndRdtn) is the largest value in
      !   RdtnTime(:) that is less than or equal to ZTime and find the amount of
      !   time remaining from this calculation:
      ! NOTE: ZTime is scaled by OnePlusEps to ensure that IndRdtn increments in
      !       steps of 1 when RdtnDT = DT, even in the presence of numerical
      !       precision errors.

      IndRdtn   = FLOOR ( ( ZTime*OnePlusEps )/FP_Data%RdtnDT )
      RdtnRmndr = ZTime - ( IndRdtn*FP_Data%RdtnDT ) ! = ZTime - RdtnTime(IndRdtn); however, RdtnTime(:) has a maximum index of NStepRdtn-1


      ! Save the new values of XD in XDHistory if:
      !   (1) we are on the initialization pass where ZTime = 0.0,
      !   (2) we have increased in time by at least RdtnDT, or
      !   (3) the time has not changed since the last time we have increased in
      !       time by at least RdtnDT (i.e., on a call to the corrector)
      !   When saving the new values, interpolate to find all of the values
      !   between index LastIndRdtn and index IndRdtn.  Also, if the XDHistory
      !   array is full, use MOD(Index,NStepRdtn1) to replace the oldest values
      !   with the newest values:
      ! NOTE: When IndRdtn > LastIndRdtn, IndRdtn will equal           LastIndRdtn + 1 if DT <= RdtnDT;
      !       When IndRdtn > LastIndRdtn, IndRdtn will be greater than LastIndRdtn + 1 if DT >  RdtnDT.
!BJJ: this needs a better check so that it is ALWAYS done (MATLAB/Simulink could possibly avoid this step by starting at ZTime>0, OR there may be some numerical issues where this is NOT EXACTLY zero)
      IF ( ZTime == 0.0 )  THEN              ! (1) .TRUE. if we are on the initialization pass where ZTime = 0.0 (and IndRdtn = 0)

         DO J = 1,6  ! Loop through all platform DOFs
            FP_Data%XDHistory(IndRdtn,J) = XD(J)
         ENDDO       ! J - All platform DOFs

         FP_Data%LastIndRdtn  =     IndRdtn       ! Save the value of                           IndRdtn for the next call to this routine (in this case IndRdtn = 0)

      ELSEIF ( IndRdtn > FP_Data%LastIndRdtn )  THEN ! (2) .TRUE. if we have increased in time by at least RdtnDT

         DO J = 1,6                       ! Loop through all platform DOFs

            IncrmntXD = ( FP_Data%RdtnDT/( ZTime -     ( FP_Data%LastIndRdtn*FP_Data%RdtnDT ) ) ) * &
                        ( XD(J) - FP_Data%XDHistory(MOD(FP_Data%LastIndRdtn ,FP_Data%NStepRdtn1),J) )

            DO K = FP_Data%LastIndRdtn +1,IndRdtn ! Loop through all radiation time steps where the time history of XD has yet to be stored
               FP_Data%XDHistory(MOD(K,FP_Data%NStepRdtn1),J) = FP_Data%XDHistory(MOD(FP_Data%LastIndRdtn ,FP_Data%NStepRdtn1),J) &
                                                                              + ( K - FP_Data%LastIndRdtn  )*IncrmntXD
            ENDDO                         ! K - All radiation time steps where the time history of XD has yet to be stored

         ENDDO                            ! J - All platform DOFs

         FP_Data%LastIndRdtn2 = FP_Data%LastIndRdtn       ! Save the value of                       LastIndRdtn for the next call to this routine
         FP_Data%LastIndRdtn  =             IndRdtn       ! Save the value of                           IndRdtn for the next call to this routine
         FP_Data%LastTime     = ZTime                     ! Save the value of ZTime associated with LastIndRdtn for the next call to this routine

!BJJ: this needs a better check in case there may be some numerical issues where this is NOT EXACTLY the same...
      ELSEIF ( ZTime == FP_Data%LastTime )  THEN     ! (3). .TRUE. if the time has not changed since the last time we have increased in time by at least RdtnDt (i.e., on a call to the corrector)

         DO J = 1,6                       ! Loop through all platform DOFs

            IncrmntXD = ( FP_Data%RdtnDT/( ZTime -    ( FP_Data%LastIndRdtn2*FP_Data%RdtnDT ) ) ) * &
                        ( XD(J) - FP_Data%XDHistory(MOD(FP_Data%LastIndRdtn2,FP_Data%NStepRdtn1),J) )

            DO K = FP_Data%LastIndRdtn2+1,IndRdtn ! Loop through all radiation time steps where the time history of XD should be updated
               FP_Data%XDHistory(MOD(K,FP_Data%NStepRdtn1),J) = FP_Data%XDHistory(MOD(FP_Data%LastIndRdtn2,FP_Data%NStepRdtn1),J) &
                                                                              + ( K - FP_Data%LastIndRdtn2 )*IncrmntXD
            ENDDO                         ! K - All radiation time steps where the time history of XD should be updated

         ENDDO                            ! J - All platform DOFs

      ENDIF


      ! Perform numerical convolution to determine the load contribution from wave
      !   radiation damping:

      DO I = 1,6                 ! Loop through all wave radiation damping forces and moments

         F_RdtnDT   (I) = 0.0
         F_RdtnRmndr(I) = 0.0

         DO J = 1,6              ! Loop through all platform DOFs

            DO K = MAX(0,IndRdtn-FP_Data%NStepRdtn  ),IndRdtn-1  ! Loop through all NStepRdtn time steps in the radiation Kernel (less than NStepRdtn time steps are used when ZTime < RdtnTmax)
               F_RdtnDT   (I) = F_RdtnDT   (I) - FP_Data%RdtnKrnl(IndRdtn-1-K,I,J)*FP_Data%XDHistory(MOD(K,FP_Data%NStepRdtn1),J)
            ENDDO                                        ! K - All NStepRdtn time steps in the radiation Kernel (less than NStepRdtn time steps are used when ZTime < RdtnTmax)

            DO K = MAX(0,IndRdtn-FP_Data%NStepRdtn+1),IndRdtn    ! Loop through all NStepRdtn time steps in the radiation Kernel (less than NStepRdtn time steps are used when ZTime < RdtnTmax)
               F_RdtnRmndr(I) = F_RdtnRmndr(I) - FP_Data%RdtnKrnl(IndRdtn  -K,I,J)*FP_Data%XDHistory(MOD(K,FP_Data%NStepRdtn1),J)
            ENDDO                                        ! K - All NStepRdtn time steps in the radiation Kernel (less than NStepRdtn time steps are used when ZTime < RdtnTmax)

         ENDDO                   ! J - All platform DOFs

         F_Rdtn     (I) = ( FP_Data%RdtnDT - RdtnRmndr )*F_RdtnDT(I) + RdtnRmndr*F_RdtnRmndr(I)

      ENDDO                      ! I - All wave radiation damping forces and moments


   ELSE                 ! We must not be modeling wave radiation damping.


      ! Set the total load contribution from radiation damping to zero:

      F_Rdtn        (:) = 0.0


   ENDIF




      ! Compute the load contribution from viscous drag using the viscous drag
      !   term from Morison's equation:
      ! NOTE: It is inconsistent to use stretching (which is a nonlinear
      !       correction) for the viscous drag term in Morison's equation while
      !       not accounting for stretching in the diffraction and radiation
      !       problems (according to Paul Sclavounos, there are such corrections).
      !       Instead, the viscous drag term from Morison's equation is computed
      !       by integrating up to the MSL, regardless of the instantaneous free
      !       surface elevation.  Also, the undisplaced platform location and
      !       orientation are used in the calculations (when finding the incident
      !       wave velocity relative to the current platform node, for example)
      !       because this is the approach most consistent with the linearized
      !       implementation of the diffraction and radiation problems.  Finally,
      !       the viscous drag load only contains components in the surge, sway,
      !       pitch, and roll directions.  There is no viscous drag in the heave
      !       or yaw directions.

   F_Viscous = 0.0   ! First initialize this force to zero.

   DO JNode = 1,WaveDat%NWaveKin0  ! Loop through the platform nodes / elements


      ! Compute the velocity of the incident waves in the xi- (1) and yi- (2)
      !   directions, respectively, at the current platform node and time:

      DO K = 1,2     ! Loop through the xi- (1) and yi- (2) directions
         WaveVelocity0(K) = WaveVelocity ( JNode, K, ZTime, WaveDat, ErrStat )
         IF ( ErrStat /= 0 ) RETURN
      ENDDO          ! K - The xi- (1) and yi- (2) directions


      ! Compute the velocity of the current platform node in the xi- (1) and yi-
      !   (2) directions, respectively, at the current time:

      PtfmVelocity(1) = XD(1) + XD(5)*WaveDat%WaveKinzi0(JNode)
      PtfmVelocity(2) = XD(2) - XD(4)*WaveDat%WaveKinzi0(JNode)


      ! Compute the magnitude of the horizontal incident wave velocity relative to
      !   the current platform node at the current time:

      MagVRel = SQRT(   ( WaveVelocity0(1) - PtfmVelocity(1) )**2 &
                      + ( WaveVelocity0(2) - PtfmVelocity(2) )**2   )

!!bjj: stuff for OC3
!IF (     WaveDat%WaveKinzi0(JNode) > - 4.0 )  THEN
!	FP_Data%PtfmDiam = 6.5
!ELSEIF ( WaveDat%WaveKinzi0(JNode) < -12.0 )  THEN
!	FP_Data%PtfmDiam = 9.4
!ELSE
!	FP_Data%PtfmDiam = 9.4 + ( 6.5 - 9.4 )*( WaveDat%WaveKinzi0(JNode) + 12.0 )/( -4.0 + 12.0 )
!ENDIF
!!bjj end of stuff for OC3

      ! Compute the viscous drag force in the xi- (1) and yi- (2) directions,
      !   respectively, on the current platform element at the current time:

      DO K = 1,2     ! Loop through the xi- (1) and yi- (2) directions
         ViscousForce (K) = 0.5*FP_Data%PtfmCD *WaveDat%WtrDens *FP_Data%PtfmDiam*( WaveVelocity0(K) - PtfmVelocity(K) ) &
                               *MagVRel*WaveDat%DZNodes(JNode)
      ENDDO          ! K - The xi- (1) and yi- (2) directions


      ! Compute the portion of the viscous drag load on the platform associated
      !   with the current element:

      F_Viscous(1) = F_Viscous(1) + ViscousForce(1)                              ! surge component
      F_Viscous(2) = F_Viscous(2) + ViscousForce(2)                              ! sway  component
      F_Viscous(4) = F_Viscous(4) - ViscousForce(2)*WaveDat%WaveKinzi0(JNode)    ! roll  component
      F_Viscous(5) = F_Viscous(5) + ViscousForce(1)*WaveDat%WaveKinzi0(JNode)    ! pitch component


   ENDDO                   ! JNode - Platform nodes / elements




      ! Compute the load contribution from mooring lines:

   SELECT CASE ( FP_Data%LineMod ) ! Which incident wave kinematics model are we using?

   CASE ( 0 )              ! None.


      ! Set the total load contribution from all mooring lines to zero:

      F_Lines(:) = 0.0



   CASE ( 1 )              ! Standard quasi-static.


      ! First initialize the total load contribution from all mooring lines to
      !   zero:

      F_Lines(:) = 0.0


      ! Get the transformation matrix, TransMat, from the inertial frame to the
      !   tower base / platform coordinate system:

      CALL SmllRotTrans ( 'platform displacement', X(4), X(5), X(6), TransMat )


      DO I = 1,FP_Data%NumLines ! Loop through all mooring lines


      ! Transform the fairlead location from the platform to the inertial
      !   frame coordinate system:
      ! NOTE: TransMat^T = TransMat^-1 where ^T = matrix transpose and ^-1 =
      !       matrix inverse.

         LFairxiRef = TransMat(1,1)*FP_Data%MooringLine(I)%LFairxt + TransMat(2,1)*FP_Data%MooringLine(I)%LFairyt &
                                                                   + TransMat(3,1)*FP_Data%MooringLine(I)%LFairzt
         LFairyiRef = TransMat(1,2)*FP_Data%MooringLine(I)%LFairxt + TransMat(2,2)*FP_Data%MooringLine(I)%LFairyt &
                                                                   + TransMat(3,2)*FP_Data%MooringLine(I)%LFairzt
         LFairziRef = TransMat(1,3)*FP_Data%MooringLine(I)%LFairxt + TransMat(2,3)*FP_Data%MooringLine(I)%LFairyt &
                                                                   + TransMat(3,3)*FP_Data%MooringLine(I)%LFairzt

         LFairxi    = X(1) + LFairxiRef
         LFairyi    = X(2) + LFairyiRef
         LFairzi    = X(3) + LFairziRef


      ! Transform the fairlead location from the inertial frame coordinate system
      !   to the local coordinate system of the current line (this coordinate
      !   system lies at the current anchor, Z being vertical, and X directed from
      !   current anchor to the current fairlead).  Also, compute the orientation
      !   of this local coordinate system:

!bjj: note that i removed the .0 after the 2.0 on the exponents so that the compiler knows that these exponents are integers and it is therefore easier to compute!
         XF         = SQRT( ( LFairxi - FP_Data%MooringLine(I)%LAnchxi )**2 + ( LFairyi - FP_Data%MooringLine(I)%LAnchyi )**2 )
         ZF         =         LFairzi - FP_Data%MooringLine(I)%LAnchzi

!bjj: should this be compared with epsilon instead of 0.0 to avoid numerical instability?
         IF ( XF == 0.0 )  THEN  ! .TRUE. if the current mooring line is exactly vertical; thus, the solution below is ill-conditioned because the orientation is undefined; so set it such that the tensions and nodal positions are only vertical
            COSPhi  = 0.0
            SINPhi  = 0.0
         ELSE                    ! The current mooring line must not be vertical; use simple trigonometry
            COSPhi  =       ( LFairxi - FP_Data%MooringLine(I)%LAnchxi )/XF
            SINPhi  =       ( LFairyi - FP_Data%MooringLine(I)%LAnchyi )/XF
         ENDIF


      ! Solve the analytical, static equilibrium equations for a catenary (or
      !   taut) mooring line with seabed interaction in order to find the
      !   horizontal and vertical tensions at the fairlead in the local coordinate
      !   system of the current line:
      ! NOTE: The values for the horizontal and vertical tensions at the fairlead
      !       from the previous time step are used as the initial guess values at
      !       at this time step (because the LAnchHTe(:) and LAnchVTe(:) arrays
      !       are stored in a module and thus their values are saved from CALL to
      !       CALL).

         CALL Catenary ( XF                                , ZF                                 , &
                         FP_Data%MooringLine(I)%LUnstrLen  , FP_Data%MooringLine(I)%LEAStff     , &
                         FP_Data%MooringLine(I)%LFldWght   , FP_Data%MooringLine(I)%LSeabedCD   , &
                         FP_Data%MooringLine(I)%LTenTol    , FP_Data%MooringLine(I)%LFairHTe    , &
                         FP_Data%MooringLine(I)%LFairVTe   , FP_Data%MooringLine(I)%LAnchHTe    , &
                         FP_Data%MooringLine(I)%LAnchVTe   , FP_Data%MooringLine(I)%LineNodes   , &
                         FP_Data%MooringLine(I)%LSNodes (:), FP_Data%MooringLine(I)%LNodesX  (:), &
                         FP_Data%MooringLine(I)%LNodesZ (:), FP_Data%MooringLine(I)%LNodesTe (:)  )


      ! Transform the positions of each node on the current line from the local
      !   coordinate system of the current line to the inertial frame coordinate
      !   system:

         DO J = 1,FP_Data%MooringLine(I)%LineNodes ! Loop through all nodes per line where the line position and tension can be output
            FP_Data%MooringLine(I)%LNodesPi(J,1) = FP_Data%MooringLine(I)%LAnchxi + FP_Data%MooringLine(I)%LNodesX(J)*COSPhi
            FP_Data%MooringLine(I)%LNodesPi(J,2) = FP_Data%MooringLine(I)%LAnchyi + FP_Data%MooringLine(I)%LNodesX(J)*SINPhi
            FP_Data%MooringLine(I)%LNodesPi(J,3) = FP_Data%MooringLine(I)%LAnchzi + FP_Data%MooringLine(I)%LNodesZ(J)
         ENDDO              ! J - All nodes per line where the line position and tension can be output


      ! Compute the portion of the mooring system load on the platform associated
      !   with the current line:

         LFairxiTe  = FP_Data%MooringLine(I)%LFairHTe*COSPhi
         LFairyiTe  = FP_Data%MooringLine(I)%LFairHTe*SINPhi
         LFairziTe  = FP_Data%MooringLine(I)%LFairVTe

         F_Lines(1) = F_Lines(1) -            LFairxiTe
         F_Lines(2) = F_Lines(2) -            LFairyiTe
         F_Lines(3) = F_Lines(3) -            LFairziTe
         F_Lines(4) = F_Lines(4) - LFairyiRef*LFairziTe + LFairziRef*LFairyiTe
         F_Lines(5) = F_Lines(5) - LFairziRef*LFairxiTe + LFairxiRef*LFairziTe
         F_Lines(6) = F_Lines(6) - LFairxiRef*LFairyiTe + LFairyiRef*LFairxiTe


      ENDDO             ! I - All mooring lines



   CASE ( 2 )              ! User-defined mooring lines.


      ! CALL the user-defined platform loading model:
!bjj start of proposed change v1.00.00a-bjj
!UserLine is now a module
!rm      CALL UserLine ( X              , ZTime    , FP_Data%DirRoot        , F_Lines        , &
!rm                      FP_Data%NumLines, LineNodes, LFairHTe       , LFairVTe       , &
!rm                      LAnchHTe       , LAnchVTe , LNodesPi(:,:,1), LNodesPi(:,:,2), &
!rm                      LNodesPi(:,:,3), LNodesTe                         )
                      
      CALL UserLine ( X, ZTime, F_Lines, FP_Data%UserLine_Data, ErrStat )
      
      IF ( ErrStat /= 0 ) THEN
         CALL ProgAbort( ' Error in UserLine subroutine.', TrapErrors = .TRUE. ) 
         RETURN
      END IF
      
      
      DO I = 1,FP_Data%NumLines ! Loop through all mooring lines

         CALL UserLine_GetOutputs( I,                                    FP_Data%MooringLine(I)%LineNodes,     &
                                   FP_Data%MooringLine(I)%LFairHTe,      FP_Data%MooringLine(I)%LFairVTe,      &
                                   FP_Data%MooringLine(I)%LAnchHTe,      FP_Data%MooringLine(I)%LAnchVTe,      &
                                   FP_Data%MooringLine(I)%LNodesPi(:,1), FP_Data%MooringLine(I)%LNodesPi(:,2), &
                                   FP_Data%MooringLine(I)%LNodesPi(:,3), FP_Data%MooringLine(I)%LNodesTe,      &
                                   FP_Data%UserLine_Data,                ErrStat                               )         

         IF ( ErrStat /= 0 ) THEN
            CALL ProgAbort( ' Error in UserLine_GetOutputs subroutine.', TrapErrors = .TRUE. ) 
            RETURN
         END IF

!         FP_Data%MooringLine(I)%LFairHTe      = FP_Data%UserLine_Data%MooringLine(I)%FairHTen
!         FP_Data%MooringLine(I)%LFairVTe      = FP_Data%UserLine_Data%MooringLine(I)%FairVTen
!         FP_Data%MooringLine(I)%LAnchHTe      = FP_Data%UserLine_Data%MooringLine(I)%AnchHTen
!         FP_Data%MooringLine(I)%LAnchVTe      = FP_Data%UserLine_Data%MooringLine(I)%AnchVTen
!         FP_Data%MooringLine(I)%LNodesPi(:,1) = FP_Data%UserLine_Data%MooringLine(I)%Nodesxi               ! these arrays need to be allocated to the same size at initialization!!!
!         FP_Data%MooringLine(I)%LNodesPi(:,2) = FP_Data%UserLine_Data%MooringLine(I)%Nodesyi
!         FP_Data%MooringLine(I)%LNodesPi(:,3) = FP_Data%UserLine_Data%MooringLine(I)%Nodeszi
!         FP_Data%MooringLine(I)%LNodesTe      = FP_Data%UserLine_Data%MooringLine(I)%NodesTen
         
      END DO      
           
!bjj end of proposed change

   ENDSELECT




      ! Total up all of the loads that do not depend on platform acceleration:

   PtfmFt = F_Waves + F_HS + F_Rdtn + F_Viscous + F_Lines   ! This is an array operation.


!!bjj start of stuff for OC3
!PtfmFt(1) = PtfmFt(1) -   100000.0*XD(1)
!PtfmFt(2) = PtfmFt(2) -   100000.0*XD(2)
!PtfmFt(3) = PtfmFt(3) -   130000.0*XD(3)
!PtfmFt(6) = PtfmFt(6) - 13000000.0*XD(6)
!!bjj end of stuff for OC3

      ! Set the platform added mass matrix, PtfmAM, to be the infinite-frequency
      !   limit of the frequency-dependent hydrodynamic added mass matrix,
      !   HdroAdMsI:

   PtfmAM = FP_Data%HdroAdMsI   ! This is an array operation.




!JASON: USE THIS TO TEST RELATIVE MAGNITUDES:WRITE (*,*) ZTime, F_Waves(5), F_Rdtn(5), F_Viscous(5), F_Lines(5)   !JASON:USE THIS TO TEST RELATIVE MAGNITUDES:
   RETURN
!jmj Start of proposed change.  v6.02b-jmj  15-Nov-2006.
!jmj Replace the hard-coded mooring line restoring calculation with a general
!jmj   purpose, quasi-static solution based on the analytical catenary cable
!jmj   equations with seabed interaction:
   CONTAINS
!=======================================================================
!JASON: SHOULD THIS ROUTINE (Catenary) BE PLACED IN NWTC_Subs OR IN ITS OWN DLL?
      SUBROUTINE Catenary ( XF_In, ZF_In, L_In  , EA_In, &
                            W_In , CB_In, Tol_In, HF_In, &
                            VF_In, HA_In, VA_In , N    , &
                            s_In , X_In , Z_In  , Te_In    )


         ! This routine solves the analytical, static equilibrium equations
         ! for a catenary (or taut) mooring line with seabed interaction.
         ! Stretching of the line is accounted for, but bending stiffness
         ! is not.  Given the mooring line properties and the fairlead
         ! position relative to the anchor, this routine finds the line
         ! configuration and tensions.  Since the analytical solution
         ! involves two nonlinear equations (XF and  ZF) in two unknowns
         ! (HF and VF), a Newton-Raphson iteration scheme is implemented in
         ! order to solve for the solution.  The values of HF and VF that
         ! are passed into this routine are used as the initial guess in
         ! the iteration.  The Newton-Raphson iteration is only accurate in
         ! double precision, so all of the input/output arguments are
         ! converteds to/from double precision from/to default precision.


      IMPLICIT                        NONE


         ! Passed Variables:

      INTEGER,    INTENT(IN   )    :: N                                               ! Number of nodes where the line position and tension can be output (-)

      REAL(ReKi), INTENT(IN   )    :: CB_In                                           ! Coefficient of seabed static friction drag (a negative value indicates no seabed) (-)
      REAL(ReKi), INTENT(IN   )    :: EA_In                                           ! Extensional stiffness of line (N)
      REAL(ReKi), INTENT(  OUT)    :: HA_In                                           ! Effective horizontal tension in line at the anchor   (N)
      REAL(ReKi), INTENT(INOUT)    :: HF_In                                           ! Effective horizontal tension in line at the fairlead (N)
      REAL(ReKi), INTENT(IN   )    :: L_In                                            ! Unstretched length of line (meters)
      REAL(ReKi), INTENT(IN   )    :: s_In     (N)                                    ! Unstretched arc distance along line from anchor to each node where the line position and tension can be output (meters)
      REAL(ReKi), INTENT(  OUT)    :: Te_In    (N)                                    ! Effective line tensions at each node (N)
      REAL(ReKi), INTENT(IN   )    :: Tol_In                                          ! Convergence tolerance within Newton-Raphson iteration specified as a fraction of tension (-)
      REAL(ReKi), INTENT(  OUT)    :: VA_In                                           ! Effective vertical   tension in line at the anchor   (N)
      REAL(ReKi), INTENT(INOUT)    :: VF_In                                           ! Effective vertical   tension in line at the fairlead (N)
      REAL(ReKi), INTENT(IN   )    :: W_In                                            ! Weight of line in fluid per unit length (N/m)
      REAL(ReKi), INTENT(  OUT)    :: X_In     (N)                                    ! Horizontal locations of each line node relative to the anchor (meters)
      REAL(ReKi), INTENT(IN   )    :: XF_In                                           ! Horizontal distance between anchor and fairlead (meters)
      REAL(ReKi), INTENT(  OUT)    :: Z_In     (N)                                    ! Vertical   locations of each line node relative to the anchor (meters)
      REAL(ReKi), INTENT(IN   )    :: ZF_In                                           ! Vertical   distance between anchor and fairlead (meters)


         ! Local Variables:

      REAL(DbKi)                   :: CB                                              ! Coefficient of seabed static friction (a negative value indicates no seabed) (-)
      REAL(DbKi)                   :: CBOvrEA                                         ! = CB/EA
      REAL(DbKi)                   :: DET                                             ! Determinant of the Jacobian matrix (m^2/N^2)
      REAL(DbKi)                   :: dHF                                             ! Increment in HF predicted by Newton-Raphson (N)
      REAL(DbKi)                   :: dVF                                             ! Increment in VF predicted by Newton-Raphson (N)
      REAL(DbKi)                   :: dXFdHF                                          ! Partial derivative of the calculated horizontal distance with respect to the horizontal fairlead tension (m/N): dXF(HF,VF)/dHF
      REAL(DbKi)                   :: dXFdVF                                          ! Partial derivative of the calculated horizontal distance with respect to the vertical   fairlead tension (m/N): dXF(HF,VF)/dVF
      REAL(DbKi)                   :: dZFdHF                                          ! Partial derivative of the calculated vertical   distance with respect to the horizontal fairlead tension (m/N): dZF(HF,VF)/dHF
      REAL(DbKi)                   :: dZFdVF                                          ! Partial derivative of the calculated vertical   distance with respect to the vertical   fairlead tension (m/N): dZF(HF,VF)/dVF
      REAL(DbKi)                   :: EA                                              ! Extensional stiffness of line (N)
      REAL(DbKi)                   :: EXF                                             ! Error function between calculated and known horizontal distance (meters): XF(HF,VF) - XF
      REAL(DbKi)                   :: EZF                                             ! Error function between calculated and known vertical   distance (meters): ZF(HF,VF) - ZF
      REAL(DbKi)                   :: HA                                              ! Effective horizontal tension in line at the anchor   (N)
      REAL(DbKi)                   :: HF                                              ! Effective horizontal tension in line at the fairlead (N)
      REAL(DbKi)                   :: HFOvrW                                          ! = HF/W
      REAL(DbKi)                   :: HFOvrWEA                                        ! = HF/WEA
      REAL(DbKi)                   :: L                                               ! Unstretched length of line (meters)
      REAL(DbKi)                   :: Lamda0                                          ! Catenary parameter used to generate the initial guesses of the horizontal and vertical tensions at the fairlead for the Newton-Raphson iteration (-)
      REAL(DbKi)                   :: LMax                                            ! Maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead) (meters)
      REAL(DbKi)                   :: LMinVFOvrW                                      ! = L - VF/W
      REAL(DbKi)                   :: LOvrEA                                          ! = L/EA
      REAL(DbKi)                   :: s        (N)                                    ! Unstretched arc distance along line from anchor to each node where the line position and tension can be output (meters)
      REAL(DbKi)                   :: sOvrEA                                          ! = s(I)/EA
      REAL(DbKi)                   :: SQRT1VFOvrHF2                                   ! = SQRT( 1.0_DbKi + VFOvrHF2      )
      REAL(DbKi)                   :: SQRT1VFMinWLOvrHF2                              ! = SQRT( 1.0_DbKi + VFMinWLOvrHF2 )
      REAL(DbKi)                   :: SQRT1VFMinWLsOvrHF2                             ! = SQRT( 1.0_DbKi + VFMinWLsOvrHF*VFMinWLsOvrHF )
      REAL(DbKi)                   :: Te       (N)                                    ! Effective line tensions at each node (N)
      REAL(DbKi)                   :: Tol                                             ! Convergence tolerance within Newton-Raphson iteration specified as a fraction of tension (-)
      REAL(DbKi)                   :: VA                                              ! Effective vertical   tension in line at the anchor   (N)
      REAL(DbKi)                   :: VF                                              ! Effective vertical   tension in line at the fairlead (N)
      REAL(DbKi)                   :: VFMinWL                                         ! = VF - WL
      REAL(DbKi)                   :: VFMinWLOvrHF                                    ! = VFMinWL/HF
      REAL(DbKi)                   :: VFMinWLOvrHF2                                   ! = VFMinWLOvrHF*VFMinWLOvrHF
      REAL(DbKi)                   :: VFMinWLs                                        ! = VFMinWL + Ws
      REAL(DbKi)                   :: VFMinWLsOvrHF                                   ! = VFMinWLs/HF
      REAL(DbKi)                   :: VFOvrHF                                         ! = VF/HF
      REAL(DbKi)                   :: VFOvrHF2                                        ! = VFOvrHF*VFOvrHF
      REAL(DbKi)                   :: VFOvrWEA                                        ! = VF/WEA
      REAL(DbKi)                   :: W                                               ! Weight of line in fluid per unit length (N/m)
      REAL(DbKi)                   :: WEA                                             ! = W*EA
      REAL(DbKi)                   :: WL                                              ! Total weight of line in fluid (N): W*L
      REAL(DbKi)                   :: Ws                                              ! = W*s(I)
      REAL(DbKi)                   :: X        (N)                                    ! Horizontal locations of each line node relative to the anchor (meters)
      REAL(DbKi)                   :: XF                                              ! Horizontal distance between anchor and fairlead (meters)
      REAL(DbKi)                   :: XF2                                             ! = XF*XF
      REAL(DbKi)                   :: Z        (N)                                    ! Vertical   locations of each line node relative to the anchor (meters)
      REAL(DbKi)                   :: ZF                                              ! Vertical   distance between anchor and fairlead (meters)
      REAL(DbKi)                   :: ZF2                                             ! = ZF*ZF

      INTEGER                      :: I                                               ! Index for counting iterations or looping through line nodes (-)
      INTEGER                      :: MaxIter                                         ! Maximum number of Newton-Raphson iterations possible before giving up (-)

      LOGICAL                      :: FirstIter                                       ! Flag to determine whether or not this is the first time through the Newton-Raphson interation (flag)



         ! The Newton-Raphson iteration is only accurate in double precision, so
         !   convert the input arguments into double precision:

      CB     = REAL( CB_In    , DbKi )
      EA     = REAL( EA_In    , DbKi )
      HF     = REAL( HF_In    , DbKi )
      L      = REAL( L_In     , DbKi )
      s  (:) = REAL( s_In  (:), DbKi )
      Tol    = REAL( Tol_In   , DbKi )
      VF     = REAL( VF_In    , DbKi )
      W      = REAL( W_In     , DbKi )
      XF     = REAL( XF_In    , DbKi )
      ZF     = REAL( ZF_In    , DbKi )



         ! Abort when there is no solution or when the only possible solution is
         !   illogical:

      IF (    Tol <= 0.0_DbKi )  THEN   ! .TRUE. when the convergence tolerance is specified incorrectly

         CALL ProgAbort ( ' Convergence tolerance must be greater than zero in routine Catenary().', TrapErrors = .TRUE.)
         ErrStat = 1
         RETURN


      ELSEIF ( XF <  0.0_DbKi )  THEN   ! .TRUE. only when the local coordinate system is not computed correctly

         CALL ProgAbort ( ' The horizontal distance between an anchor and its'// &
                      ' fairlead must not be less than zero in routine Catenary().', TrapErrors = .TRUE.)
         ErrStat = 1
         RETURN


      ELSEIF ( ZF <  0.0_DbKi )  THEN   ! .TRUE. if the fairlead has passed below its anchor

         CALL ProgAbort ( ' A fairlead has passed below its anchor.', TrapErrors = .TRUE.)
         ErrStat = 1
         RETURN


      ELSEIF ( L  <= 0.0_DbKi )  THEN   ! .TRUE. when the unstretched line length is specified incorrectly

         CALL ProgAbort ( ' Unstretched length of line must be greater than zero in routine Catenary().', TrapErrors = .TRUE.)
         ErrStat = 1
         RETURN


      ELSEIF ( EA <= 0.0_DbKi )  THEN   ! .TRUE. when the unstretched line length is specified incorrectly

         CALL ProgAbort ( ' Extensional stiffness of line must be greater than zero in routine Catenary().', TrapErrors = .TRUE.)
         ErrStat = 1
         RETURN


      ELSEIF ( W  == 0.0_DbKi )  THEN   ! .TRUE. when the weight of the line in fluid is zero so that catenary solution is ill-conditioned

         CALL ProgAbort ( ' The weight of the line in fluid must not be zero. '// &
                      ' Routine Catenary() cannot solve quasi-static mooring line solution.', TrapErrors = .TRUE.)
         ErrStat = 1
         RETURN


      ELSEIF ( W  >  0.0_DbKi )  THEN   ! .TRUE. when the line will sink in fluid

         LMax      = XF - EA/W + SQRT( (EA/W)*(EA/W) + 2.0_DbKi*ZF*EA/W )  ! Compute the maximum stretched length of the line with seabed interaction beyond which the line would have to double-back on itself; here the line forms an "L" between the anchor and fairlead (i.e. it is horizontal along the seabed from the anchor, then vertical to the fairlead)

         IF ( ( L  >=  LMax   ) .AND. ( CB >= 0.0_DbKi ) )  THEN  ! .TRUE. if the line is as long or longer than its maximum possible value with seabed interaction
            CALL ProgAbort ( ' Unstretched mooring line length too large. '// &
                         ' Routine Catenary() cannot solve quasi-static mooring line solution.', TrapErrors = .TRUE.)
            ErrStat = 1
            RETURN
         END IF

      ENDIF



         ! Initialize some commonly used terms that don't depend on the iteration:

      WL      =          W  *L
      WEA     =          W  *EA
      LOvrEA  =          L  /EA
      CBOvrEA =          CB /EA
      MaxIter = INT(1.0_DbKi/Tol)   ! Smaller tolerances may take more iterations, so choose a maximum inversely proportional to the tolerance



         ! To avoid an ill-conditioned situation, ensure that the initial guess for
         !   HF is not less than or equal to zero.  Similarly, avoid the problems
         !   associated with having exactly vertical (so that HF is zero) or exactly
         !   horizontal (so that VF is zero) lines by setting the minimum values
         !   equal to the tolerance.  This prevents us from needing to implement
         !   the known limiting solutions for vertical or horizontal lines (and thus
         !   complicating this routine):

      HF = MAX( HF, Tol )
      XF = MAX( XF, Tol )
      ZF = MAX( ZF, TOl )



         ! Solve the analytical, static equilibrium equations for a catenary (or
         !   taut) mooring line with seabed interaction:

         ! Begin Newton-Raphson iteration:

      I         = 1        ! Initialize iteration counter
      FirstIter = .TRUE.   ! Initialize iteration flag

      DO


         ! Initialize some commonly used terms that depend on HF and VF:

         VFMinWL            = VF - WL
         LMinVFOvrW         = L  - VF/W
         HFOvrW             =      HF/W
         HFOvrWEA           =      HF/WEA
         VFOvrWEA           =      VF/WEA
         VFOvrHF            =      VF/HF
         VFMinWLOvrHF       = VFMinWL/HF
         VFOvrHF2           = VFOvrHF     *VFOvrHF
         VFMinWLOvrHF2      = VFMinWLOvrHF*VFMinWLOvrHF
         SQRT1VFOvrHF2      = SQRT( 1.0_DbKi + VFOvrHF2      )
         SQRT1VFMinWLOvrHF2 = SQRT( 1.0_DbKi + VFMinWLOvrHF2 )


         ! Compute the error functions (to be zeroed) and the Jacobian matrix
         !   (these depend on the anticipated configuration of the mooring line):

         IF ( ( CB <  0.0_DbKi ) .OR. ( W  <  0.0_DbKi ) .OR. ( VFMinWL >  0.0_DbKi ) )  THEN   ! .TRUE. when no portion of the line      rests on the seabed

            EXF    = (   LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                       &
                       - LOG( VFMinWLOvrHF +               SQRT1VFMinWLOvrHF2 )                                         )*HFOvrW &
                   + LOvrEA*  HF                         - XF
            EZF    = (                                     SQRT1VFOvrHF2                                              &
                       -                                   SQRT1VFMinWLOvrHF2                                           )*HFOvrW &
                   + LOvrEA*( VF - 0.5_DbKi*WL )         - ZF

            dXFdHF = (   LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                       &
                       - LOG( VFMinWLOvrHF +               SQRT1VFMinWLOvrHF2 )                                         )/     W &
                   - (      ( VFOvrHF      + VFOvrHF2     /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      ) &
                       -    ( VFMinWLOvrHF + VFMinWLOvrHF2/SQRT1VFMinWLOvrHF2 )/( VFMinWLOvrHF + SQRT1VFMinWLOvrHF2 )   )/     W &
                   + LOvrEA
            dXFdVF = (      ( 1.0_DbKi     + VFOvrHF      /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      ) &
                       -    ( 1.0_DbKi     + VFMinWLOvrHF /SQRT1VFMinWLOvrHF2 )/( VFMinWLOvrHF + SQRT1VFMinWLOvrHF2 )   )/     W
            dZFdHF = (                                     SQRT1VFOvrHF2                                              &
                       -                                   SQRT1VFMinWLOvrHF2                                           )/     W &
                   - (                       VFOvrHF2     /SQRT1VFOvrHF2                                              &
                       -                     VFMinWLOvrHF2/SQRT1VFMinWLOvrHF2                                           )/     W
            dZFdVF = (                       VFOvrHF      /SQRT1VFOvrHF2                                              &
                       -                     VFMinWLOvrHF /SQRT1VFMinWLOvrHF2                                           )/     W &
                   + LOvrEA


         ELSEIF (                                           -CB*VFMinWL <  HF         )  THEN   ! .TRUE. when a  portion of the line      rests on the seabed and the anchor tension is nonzero

            EXF    =     LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          *HFOvrW &
                   - 0.5_DbKi*CBOvrEA*W*  LMinVFOvrW*LMinVFOvrW                                                                  &
                   + LOvrEA*  HF           + LMinVFOvrW  - XF
            EZF    = (                                     SQRT1VFOvrHF2                                   - 1.0_DbKi   )*HFOvrW &
                   + 0.5_DbKi*VF*VFOvrWEA                - ZF

            dXFdHF =     LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          /     W &
                   - (      ( VFOvrHF      + VFOvrHF2     /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2        ) )/     W &
                   + LOvrEA
            dXFdVF = (      ( 1.0_DbKi     + VFOvrHF      /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2        ) )/     W &
                   + CBOvrEA*LMinVFOvrW - 1.0_DbKi/W
            dZFdHF = (                                     SQRT1VFOvrHF2                                   - 1.0_DbKi &
                       -                     VFOvrHF2     /SQRT1VFOvrHF2                                                )/     W
            dZFdVF = (                       VFOvrHF      /SQRT1VFOvrHF2                                                )/     W &
                   + VFOvrWEA


         ELSE                                                ! 0.0_DbKi <  HF  <= -CB*VFMinWL   !             A  portion of the line must rest  on the seabed and the anchor tension is    zero

            EXF    =     LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          *HFOvrW &
                   - 0.5_DbKi*CBOvrEA*W*( LMinVFOvrW*LMinVFOvrW - ( LMinVFOvrW - HFOvrW/CB )*( LMinVFOvrW - HFOvrW/CB ) )        &
                   + LOvrEA*  HF           + LMinVFOvrW  - XF
            EZF    = (                                     SQRT1VFOvrHF2                                   - 1.0_DbKi   )*HFOvrW &
                   + 0.5_DbKi*VF*VFOvrWEA                - ZF

            dXFdHF =     LOG( VFOvrHF      +               SQRT1VFOvrHF2      )                                          /     W &
                   - (      ( VFOvrHF      + VFOvrHF2     /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      )   )/     W &
                   + LOvrEA - ( LMinVFOvrW - HFOvrW/CB )/EA
            dXFdVF = (      ( 1.0_DbKi     + VFOvrHF      /SQRT1VFOvrHF2      )/( VFOvrHF      + SQRT1VFOvrHF2      )   )/     W &
                   + HFOvrWEA           - 1.0_DbKi/W
            dZFdHF = (                                     SQRT1VFOvrHF2                                   - 1.0_DbKi &
                       -                     VFOvrHF2     /SQRT1VFOvrHF2                                                )/     W
            dZFdVF = (                       VFOvrHF      /SQRT1VFOvrHF2                                                )/     W &
                   + VFOvrWEA


         ENDIF


         ! Compute the determinant of the Jacobian matrix and the incremental
         !   tensions predicted by Newton-Raphson:

         DET = dXFdHF*dZFdVF - dXFdVF*dZFdHF

         dHF = ( -dZFdVF*EXF + dXFdVF*EZF )/DET    ! This is the incremental change in horizontal tension at the fairlead as predicted by Newton-Raphson
         dVF = (  dZFdHF*EXF - dXFdHF*EZF )/DET    ! This is the incremental change in vertical   tension at the fairlead as predicted by Newton-Raphson

         dHF = dHF*( 1.0_DbKi - Tol*I )            ! Reduce dHF by factor (between 1 at I = 1 and 0 at I = MaxIter) that reduces linearly with iteration count to ensure that we converge on a solution even in the case were we obtain a nonconvergent cycle about the correct solution (this happens, for example, if we jump to quickly between a taut and slack catenary)
         dVF = dVF*( 1.0_DbKi - Tol*I )            ! Reduce dHF by factor (between 1 at I = 1 and 0 at I = MaxIter) that reduces linearly with iteration count to ensure that we converge on a solution even in the case were we obtain a nonconvergent cycle about the correct solution (this happens, for example, if we jump to quickly between a taut and slack catenary)

         dHF = MAX( dHF, ( Tol - 1.0_DbKi )*HF )   ! To avoid an ill-conditioned situation, make sure HF does not go less than or equal to zero by having a lower limit of Tol*HF [NOTE: the value of dHF = ( Tol - 1.0_DbKi )*HF comes from: HF = HF + dHF = Tol*HF when dHF = ( Tol - 1.0_DbKi )*HF]


         ! Check if we have converged on a solution, or restart the iteration, or
         !   Abort if we cannot find a solution:

         IF ( ( ABS(dHF) <= ABS(Tol*HF) ) .AND. ( ABS(dVF) <= ABS(Tol*VF) ) )  THEN ! .TRUE. if we have converged; stop iterating! [The converge tolerance, Tol, is a fraction of tension]

            EXIT


         ELSEIF ( ( I == MaxIter )        .AND. (       FirstIter         ) )  THEN ! .TRUE. if we've iterated MaxIter-times for the first time;

         ! Perhaps we failed to converge because our initial guess was too far off.
         !   (This could happen, for example, while linearizing a model via large
         !   pertubations in the DOFs.)  Instead, use starting values documented in:
         !   Peyrot, Alain H. and Goulois, A. M., "Analysis Of Cable Structures,"
         !   Computers & Structures, Vol. 10, 1979, pp. 805-813:
         ! NOTE: We don't need to check if the current mooring line is exactly
         !       vertical (i.e., we don't need to check if XF == 0.0), because XF is
         !       limited by the tolerance above.

            XF2 = XF*XF
            ZF2 = ZF*ZF

            IF ( L <= SQRT( XF2 + ZF2 ) )  THEN ! .TRUE. if the current mooring line is taut
               Lamda0 = 0.2_DbKi
            ELSE                                ! The current mooring line must be slack and not vertical
               Lamda0 = SQRT( 3.0_DbKi*( ( L*L - ZF2 )/XF2 - 1.0_DbKi ) )
            ENDIF

            HF  = MAX( ABS( 0.5_DbKi*W*  XF/     Lamda0      ), Tol )   ! As above, set the lower limit of the guess value of HF to the tolerance
            VF  =           0.5_DbKi*W*( ZF/TANH(Lamda0) + L )


         ! Restart Newton-Raphson iteration:

            I         = 0
            FirstIter = .FALSE.
            dHF       = 0.0_DbKi
            dVF       = 0.0_DbKi


         ELSEIF ( ( I == MaxIter )        .AND. ( .NOT. FirstIter         ) )  THEN ! .TRUE. if we've iterated as much as we can take without finding a solution; Abort

            CALL ProgAbort ( ' Iteration not convergent. '// &
                         ' Routine Catenary() cannot solve quasi-static mooring line solution.', TrapErrors = .TRUE.)
            ErrStat = 1
            RETURN


         ENDIF


         ! Increment fairlead tensions and iteration counter so we can try again:

         HF = HF + dHF
         VF = VF + dVF

         I  = I  + 1


      ENDDO



         ! We have found a solution for the tensions at the fairlead!

         ! Now compute the tensions at the anchor and the line position and tension
         !   at each node (again, these depend on the configuration of the mooring
         !   line):

      IF ( ( CB <  0.0_DbKi ) .OR. ( W  <  0.0_DbKi ) .OR. ( VFMinWL >  0.0_DbKi ) )  THEN   ! .TRUE. when no portion of the line      rests on the seabed

         ! Anchor tensions:

         HA = HF
         VA = VFMinWL


         ! Line position and tension at each node:

         DO I = 1,N  ! Loop through all nodes where the line position and tension are to be computed

            IF ( ( s(I) <  0.0_DbKi ) .OR. ( s(I) >  L ) )  THEN
               CALL ProgAbort ( ' All line nodes must be located between the anchor ' &
                              //'and fairlead (inclusive) in routine Catenary().' )
            END IF

            Ws                  = W       *s(I)                                  ! Initialize
            VFMinWLs            = VFMinWL + Ws                                   ! some commonly
            VFMinWLsOvrHF       = VFMinWLs/HF                                    ! used terms
            sOvrEA              = s(I)    /EA                                    ! that depend
            SQRT1VFMinWLsOvrHF2 = SQRT( 1.0_DbKi + VFMinWLsOvrHF*VFMinWLsOvrHF ) ! on s(I)

            X (I)    = (   LOG( VFMinWLsOvrHF + SQRT1VFMinWLsOvrHF2 ) &
                         - LOG( VFMinWLOvrHF  + SQRT1VFMinWLOvrHF2  )   )*HFOvrW                     &
                     + sOvrEA*  HF
            Z (I)    = (                        SQRT1VFMinWLsOvrHF2   &
                         -                      SQRT1VFMinWLOvrHF2      )*HFOvrW                     &
                     + sOvrEA*(         VFMinWL + 0.5_DbKi*Ws    )
            Te(I)    = SQRT( HF*HF +    VFMinWLs*VFMinWLs )

         ENDDO       ! I - All nodes where the line position and tension are to be computed


      ELSEIF (                                           -CB*VFMinWL <  HF         )  THEN   ! .TRUE. when a  portion of the line      rests on the seabed and the anchor tension is nonzero

         ! Anchor tensions:

         HA = HF + CB*VFMinWL
         VA = 0.0_DbKi


         ! Line position and tension at each node:

         DO I = 1,N  ! Loop through all nodes where the line position and tension are to be computed

            IF ( ( s(I) <  0.0_DbKi ) .OR. ( s(I) >  L ) )  THEN
               CALL ProgAbort ( ' All line nodes must be located between the anchor ' &
                              //'and fairlead (inclusive) in routine Catenary().', TrapErrors = .TRUE.)
               ErrStat = 1
               RETURN
            END IF

            Ws                  = W       *s(I)                                  ! Initialize
            VFMinWLs            = VFMinWL + Ws                                   ! some commonly
            VFMinWLsOvrHF       = VFMinWLs/HF                                    ! used terms
            sOvrEA              = s(I)    /EA                                    ! that depend
            SQRT1VFMinWLsOvrHF2 = SQRT( 1.0_DbKi + VFMinWLsOvrHF*VFMinWLsOvrHF ) ! on s(I)

            IF (     s(I) <= LMinVFOvrW             )  THEN ! .TRUE. if this node rests on the seabed and the tension is nonzero

               X (I) = s(I)                                                                          &
                     + sOvrEA*( HF + CB*VFMinWL + 0.5_DbKi*Ws*CB )
               Z (I) = 0.0_DbKi
               Te(I) =       HF    + CB*VFMinWLs

            ELSE           ! LMinVFOvrW < s <= L            !           This node must be above the seabed

               X (I) =     LOG( VFMinWLsOvrHF + SQRT1VFMinWLsOvrHF2 )    *HFOvrW                     &
                     + sOvrEA*  HF + LMinVFOvrW                    - 0.5_DbKi*CB*VFMinWL*VFMinWL/WEA
               Z (I) = ( - 1.0_DbKi           + SQRT1VFMinWLsOvrHF2     )*HFOvrW                     &
                     + sOvrEA*(         VFMinWL + 0.5_DbKi*Ws    ) + 0.5_DbKi*   VFMinWL*VFMinWL/WEA
               Te(I) = SQRT( HF*HF +    VFMinWLs*VFMinWLs )

            ENDIF

         ENDDO       ! I - All nodes where the line position and tension are to be computed


      ELSE                                                ! 0.0_DbKi <  HF  <= -CB*VFMinWL   !             A  portion of the line must rest  on the seabed and the anchor tension is    zero

         ! Anchor tensions:

         HA = 0.0_DbKi
         VA = 0.0_DbKi


         ! Line position and tension at each node:

         DO I = 1,N  ! Loop through all nodes where the line position and tension are to be computed

            IF ( ( s(I) <  0.0_DbKi ) .OR. ( s(I) >  L ) )  THEN
               CALL ProgAbort ( ' All line nodes must be located between the anchor ' &
                              //'and fairlead (inclusive) in routine Catenary().', TrapErrors = .TRUE.)
               ErrStat = 1
               RETURN
            END IF

            Ws                  = W       *s(I)                                  ! Initialize
            VFMinWLs            = VFMinWL + Ws                                   ! some commonly
            VFMinWLsOvrHF       = VFMinWLs/HF                                    ! used terms
            sOvrEA              = s(I)    /EA                                    ! that depend
            SQRT1VFMinWLsOvrHF2 = SQRT( 1.0_DbKi + VFMinWLsOvrHF*VFMinWLsOvrHF ) ! on s(I)

            IF (     s(I) <= LMinVFOvrW - HFOvrW/CB )  THEN ! .TRUE. if this node rests on the seabed and the tension is    zero

               X (I) = s(I)
               Z (I) = 0.0_DbKi
               Te(I) = 0.0_DbKi

            ELSEIF ( s(I) <= LMinVFOvrW             )  THEN ! .TRUE. if this node rests on the seabed and the tension is nonzero

               X (I) = s(I)                     - ( LMinVFOvrW - 0.5_DbKi*HFOvrW/CB )*HF/EA          &
                     + sOvrEA*( HF + CB*VFMinWL + 0.5_DbKi*Ws*CB ) + 0.5_DbKi*CB*VFMinWL*VFMinWL/WEA
               Z (I) = 0.0_DbKi
               Te(I) =       HF    + CB*VFMinWLs

            ELSE           ! LMinVFOvrW < s <= L            !           This node must be above the seabed

               X (I) =     LOG( VFMinWLsOvrHF + SQRT1VFMinWLsOvrHF2 )    *HFOvrW                     &
                     + sOvrEA*  HF + LMinVFOvrW - ( LMinVFOvrW - 0.5_DbKi*HFOvrW/CB )*HF/EA
               Z (I) = ( - 1.0_DbKi           + SQRT1VFMinWLsOvrHF2     )*HFOvrW                     &
                     + sOvrEA*(         VFMinWL + 0.5_DbKi*Ws    ) + 0.5_DbKi*   VFMinWL*VFMinWL/WEA
               Te(I) = SQRT( HF*HF +    VFMinWLs*VFMinWLs )

            ENDIF

         ENDDO       ! I - All nodes where the line position and tension are to be computed


      ENDIF



         ! The Newton-Raphson iteration is only accurate in double precision, so
         !   convert the output arguments back into the default precision for real
         !   numbers:

      HA_In    = REAL( HA   , ReKi )
      HF_In    = REAL( HF   , ReKi )
      Te_In(:) = REAL( Te(:), ReKi )
      VA_In    = REAL( VA   , ReKi )
      VF_In    = REAL( VF   , ReKi )
      X_In (:) = REAL( X (:), ReKi )
      Z_In (:) = REAL( Z (:), ReKi )



      RETURN
      END SUBROUTINE Catenary
!=======================================================================
!bjj start of proposed change v1.00.00a-bjj
!bjj: this is now stored in UserLine.f90
!!JASON: MOVE THIS USER-DEFINED ROUTINE (UserLine) TO THE UserSubs.f90 OF HydroDyn WHEN THE PLATFORM LOADING FUNCTIONALITY HAS BEEN DOCUMENTED!!!!!
!      SUBROUTINE UserLine ( X       , ZTime    , DirRoot , F       , &
!                            NumLines, LineNodes, FairHTen, FairVTen, &
!                            AnchHTen, AnchVTen , Nodesxi , Nodesyi , &
!                            Nodeszi , NodesTen                         )
!
!
!         ! This is a dummy routine for holding the place of a user-specified
!         ! mooring system.  Modify this code to create your own model of an
!         ! array of mooring lines.  The local variables and associated
!         ! calculations below provide a template for making this
!         ! user-specified mooring system model include a linear 6x6 restoring
!         ! matrix with offset.  These are provided as an example only and can
!         ! be modified or deleted as desired by the user without detriment to
!         ! the interface (i.e., they are not necessary for the interface).
!
!         ! The primary output of this routine is array F(:), which must
!         ! contain the 3 components of the total force from all mooring lines
!         ! (in N) acting at the platform reference and the 3 components of the
!         ! total moment from all mooring lines (in N-m) acting at the platform
!         ! reference; positive forces are in the direction of positive
!         ! platform displacement.  This primary output effects the overall
!         ! dynamic response of the system.  However, this routine must also
!         ! compute:
!         !   Array FairHTen(:)   - Effective horizontal tension at the fairlead of each mooring line
!         !   Array FairVTen(:)   - Effective vertical   tension at the fairlead of each mooring line
!         !   Array AnchHTen(:)   - Effective horizontal tension at the anchor   of each mooring line
!         !   Array AnchVTen(:)   - Effective vertical   tension at the anchor   of each mooring line
!         !   Array NodesTen(:,:) - Effective line tensions              at each node of each line
!         !   Array Nodesxi (:,:) - xi-coordinates in the inertial frame of each node of each line
!         !   Array Nodesyi (:,:) - yi-coordinates in the inertial frame of each node of each line
!         !   Array Nodeszi (:,:) - zi-coordinates in the inertial frame of each node of each line
!         ! These secondary outputs are only used to determine what to output
!         ! for the associated parameters placed in the OutList from the
!         ! primary input file.  The number of mooring lines where the fairlead
!         ! and anchor tensions can be output and the number of nodes per line
!         ! where the mooring line position and tension can be output, NumLines
!         ! and LineNodes, respectively, are additional inputs to this routine.
!
!
!      IMPLICIT                        NONE
!
!
!         ! Passed Variables:
!
!      INTEGER,    INTENT(IN )      :: LineNodes                                       ! Number of nodes per line where the mooring line position and tension can be output, (-).
!      INTEGER,    INTENT(IN )      :: NumLines                                        ! Number of mooring lines where the fairlead and anchor tensions can be output, (-).
!
!!bjj: these outputs must be in a derived type such as
!! TYPE :: LineType
!!      INTEGER                      :: LineNodes
!!      REAL(ReKi)                   :: AnchHTen                                        ! Effective horizontal tension at the anchor   of this mooring line, N.
!!      REAL(ReKi)                   :: AnchVTen                                        ! Effective vertical   tension at the anchor   of this mooring line, N.
!!      REAL(ReKi)                   :: FairHTen                                        ! Effective horizontal tension at the fairlead of this mooring line, N.
!!      REAL(ReKi)                   :: FairVTen                                        ! Effective vertical   tension at the fairlead of this mooring line, N.
!!      REAL(ReKi), ALLOCATABLE      :: NodesTen (:) !(LineNodes)                       ! Effective line tensions              at each node of this line, N.
!!      REAL(ReKi), ALLOCATABLE      :: Nodesxi  (:) !(LineNodes)                       ! xi-coordinates in the inertial frame of each node of this line, meters.
!!      REAL(ReKi), ALLOCATABLE      :: Nodesyi  (:) !(LineNodes)                       ! yi-coordinates in the inertial frame of each node of this line, meters.
!!      REAL(ReKi), ALLOCATABLE      :: Nodeszi  (:) !(LineNodes)                       ! zi-coordinates in the inertial frame of each node of this line, meters.
!!END TYPE LineType
!!
!!TYPE(LineType), intent(out)        :: LineOutputData(NumLines)
!!bjj: so that LineNodes can vary by Line
!
!      REAL(ReKi), INTENT(OUT)      :: F        (6)                                    ! The 3 components of the total force from all mooring lines (in N) acting at the platform reference and the 3 components of the total moment from all mooring lines (in N-m) acting at the platform reference; positive forces are in the direction of positive platform displacement.
!      REAL(ReKi), INTENT(OUT)      :: AnchHTen (NumLines          )                   ! Effective horizontal tension at the anchor   of each mooring line, N.
!      REAL(ReKi), INTENT(OUT)      :: AnchVTen (NumLines          )                   ! Effective vertical   tension at the anchor   of each mooring line, N.
!      REAL(ReKi), INTENT(OUT)      :: FairHTen (NumLines          )                   ! Effective horizontal tension at the fairlead of each mooring line, N.
!      REAL(ReKi), INTENT(OUT)      :: FairVTen (NumLines          )                   ! Effective vertical   tension at the fairlead of each mooring line, N.
!      REAL(ReKi), INTENT(OUT)      :: NodesTen (NumLines,LineNodes)                   ! Effective line tensions              at each node of each line, N.
!      REAL(ReKi), INTENT(OUT)      :: Nodesxi  (NumLines,LineNodes)                   ! xi-coordinates in the inertial frame of each node of each line, meters.
!      REAL(ReKi), INTENT(OUT)      :: Nodesyi  (NumLines,LineNodes)                   ! yi-coordinates in the inertial frame of each node of each line, meters.
!      REAL(ReKi), INTENT(OUT)      :: Nodeszi  (NumLines,LineNodes)                   ! zi-coordinates in the inertial frame of each node of each line, meters.
!      REAL(ReKi), INTENT(IN )      :: X        (6)                                    ! The 3 components of the translational displacement         (in m) of        the platform reference and the 3 components of the rotational displacement             (in rad) of        the platform relative to the inertial frame.
!      REAL(ReKi), INTENT(IN )      :: ZTime                                           ! Current simulation time, sec.
!
!      CHARACTER(1024), INTENT(IN ) :: DirRoot                                         ! The name of the root file including the full path to the current working directory.  This may be useful if you want this routine to write a permanent record of what it does to be stored with the simulation results: the results should be stored in a file whose name (including path) is generated by appending any suitable extension to DirRoot.
!
!
!         ! Local Variables:
!
!      REAL(ReKi)                   :: F0       (6)                                    ! Total mooring line load acting on the support platform in its undisplaced position (N, N-m)
!      REAL(ReKi)                   :: Stff     (6,6)                                  ! Linear restoring matrix from all mooring lines (kg/s^2, kg-m/s^2, kg-m^2/s^2)
!
!      INTEGER                      :: I                                               ! Generic index.
!      INTEGER                      :: J                                               ! Generic index.
!
!
!
!      F0  (1  ) = 0.0
!      F0  (2  ) = 0.0
!      F0  (3  ) = 0.0
!      F0  (4  ) = 0.0
!      F0  (5  ) = 0.0
!      F0  (6  ) = 0.0
!
!      Stff(1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!      Stff(2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!      Stff(3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!      Stff(4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!      Stff(5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!      Stff(6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
!
!!JASON: VALUES FOR OldTLP; REMOVE THIS!!!F0  (1  ) =         0.0 !JASON: VALUES FOR OldTLP; REMOVE THIS!!!
!!JASON: VALUES FOR OldTLP; REMOVE THIS!!!F0  (2  ) =         0.0 !JASON: VALUES FOR OldTLP; REMOVE THIS!!!
!!JASON: VALUES FOR OldTLP; REMOVE THIS!!!F0  (3  ) = -41050000.0 !JASON: VALUES FOR OldTLP; REMOVE THIS!!!
!!JASON: VALUES FOR OldTLP; REMOVE THIS!!!F0  (4  ) =         0.0 !JASON: VALUES FOR OldTLP; REMOVE THIS!!!
!!JASON: VALUES FOR OldTLP; REMOVE THIS!!!F0  (5  ) =         0.0 !JASON: VALUES FOR OldTLP; REMOVE THIS!!!
!!JASON: VALUES FOR OldTLP; REMOVE THIS!!!F0  (6  ) =         0.0 !JASON: VALUES FOR OldTLP; REMOVE THIS!!!
!!JASON: VALUES FOR OldTLP; REMOVE THIS!!!Stff(1,:) = (/    907000.0,        0.0,         0.0,           0.0,   -16100000.0,        0.0 /)   !JASON: VALUES FOR OldTLP; REMOVE THIS!!!
!!JASON: VALUES FOR OldTLP; REMOVE THIS!!!Stff(2,:) = (/         0.0,   907000.0,         0.0,    16100000.0,           0.0,        0.0 /)   !JASON: VALUES FOR OldTLP; REMOVE THIS!!!
!!JASON: VALUES FOR OldTLP; REMOVE THIS!!!Stff(3,:) = (/         0.0,        0.0, 213000000.0,           0.0,           0.0,        0.0 /)   !JASON: VALUES FOR OldTLP; REMOVE THIS!!!
!!JASON: VALUES FOR OldTLP; REMOVE THIS!!!Stff(4,:) = (/         0.0, 15600000.0,         0.0, 10600000000.0,           0.0,        0.0 /)   !JASON: VALUES FOR OldTLP; REMOVE THIS!!!
!!JASON: VALUES FOR OldTLP; REMOVE THIS!!!Stff(5,:) = (/ -15600000.0,        0.0,         0.0,           0.0, 10600000000.0,        0.0 /)   !JASON: VALUES FOR OldTLP; REMOVE THIS!!!
!!JASON: VALUES FOR OldTLP; REMOVE THIS!!!Stff(6,:) = (/         0.0,        0.0,         0.0,           0.0,           0.0, 82900000.0 /)   !JASON: VALUES FOR OldTLP; REMOVE THIS!!!
!F0  (1  ) = 0.0   !JASON: VALUES FOR MIT/NREL SDB; REMOVE THIS!!!
!F0  (2  ) = 0.0   !JASON: VALUES FOR MIT/NREL SDB; REMOVE THIS!!!
!F0  (3  ) = 0.0   !JASON: VALUES FOR MIT/NREL SDB; REMOVE THIS!!!
!F0  (4  ) = 0.0   !JASON: VALUES FOR MIT/NREL SDB; REMOVE THIS!!!
!F0  (5  ) = 0.0   !JASON: VALUES FOR MIT/NREL SDB; REMOVE THIS!!!
!F0  (6  ) = 0.0   !JASON: VALUES FOR MIT/NREL SDB; REMOVE THIS!!!
!Stff(1,:) = (/ 4000000.0,       0.0, 0.0, 0.0, 0.0, 0.0 /)  !JASON: VALUES FOR MIT/NREL SDB; REMOVE THIS!!!
!Stff(2,:) = (/       0.0, 4000000.0, 0.0, 0.0, 0.0, 0.0 /)  !JASON: VALUES FOR MIT/NREL SDB; REMOVE THIS!!!
!Stff(3,:) = (/       0.0,       0.0, 0.0, 0.0, 0.0, 0.0 /)  !JASON: VALUES FOR MIT/NREL SDB; REMOVE THIS!!!
!Stff(4,:) = (/       0.0,       0.0, 0.0, 0.0, 0.0, 0.0 /)  !JASON: VALUES FOR MIT/NREL SDB; REMOVE THIS!!!
!Stff(5,:) = (/       0.0,       0.0, 0.0, 0.0, 0.0, 0.0 /)  !JASON: VALUES FOR MIT/NREL SDB; REMOVE THIS!!!
!Stff(6,:) = (/       0.0,       0.0, 0.0, 0.0, 0.0, 0.0 /)  !JASON: VALUES FOR MIT/NREL SDB; REMOVE THIS!!!
!!JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!F0  (1  ) =         0.0 !JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!
!!JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!F0  (2  ) =         0.0 !JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!
!!JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!F0  (3  ) =  -6150000.0 !JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!
!!JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!F0  (4  ) =         0.0 !JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!
!!JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!F0  (5  ) =         0.0 !JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!
!!JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!F0  (6  ) =         0.0 !JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!
!!JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!Stff(1,:) = (/  15920.0,       0.0,     0.0,        0.0,   144700.0,       0.0 /)   !JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!
!!JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!Stff(2,:) = (/      0.0,   15920.0,     0.0,  -144600.0,        0.0,       0.0 /)   !JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!
!!JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!Stff(3,:) = (/      0.0,       0.0, 24930.0,        0.0,        0.0,       0.0 /)   !JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!
!!JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!Stff(4,:) = (/      0.0, -144500.0,     0.0, 38740000.0,        0.0,       0.0 /)   !JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!
!!JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!Stff(5,:) = (/ 144500.0,       0.0,     0.0,        0.0, 38740000.0,       0.0 /)   !JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!
!!JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!Stff(6,:) = (/      0.0,       0.0,     0.0,        0.0,        0.0, 2797000.0 /)   !JASON: VALUES FOR ITI BARGE; REMOVE THIS!!!
!
!      DO I = 1,6     ! Loop through all mooring line forces and moments
!            F(I) = F0(I)
!         DO J = 1,6  ! Loop through all platform DOFs
!            F(I) = F (I) - Stff(I,J)*X(J)
!         ENDDO       ! J - All platform DOFs
!      ENDDO          ! I - All mooring line forces and moments
!
!
!
!      DO I = 1,NumLines       ! Loop through all mooring lines where the fairlead and anchor tensions can be output
!
!         FairHTen   (I  ) = 0.0
!         FairVTen   (I  ) = 0.0
!         AnchHTen   (I  ) = 0.0
!         AnchVTen   (I  ) = 0.0
!
!         DO J = 1,LineNodes   ! Loop through all nodes per line where the line position and tension can be output
!
!            Nodesxi (I,J) = 0.0
!            Nodesyi (I,J) = 0.0
!            Nodeszi (I,J) = 0.0
!            NodesTen(I,J) = 0.0
!
!         ENDDO                ! J - All nodes per line where the line position and tension can be output
!
!      ENDDO                   ! I - All mooring lines where the fairlead and anchor tensions can be output
!
!
!
!      RETURN
!      END SUBROUTINE UserLine
!!=======================================================================
!!jmj End of proposed change.  v6.02b-jmj  15-Nov-2006.
!bjj end of propsoed change v1.00.00a-bjj
   END SUBROUTINE FltngPtfmLd
!=======================================================================
!bjj start of proposed change v1.00.00a-bjj
!rm   SUBROUTINE InitFltngPtfmLd ( WAMITFile  , PtfmVol0In, PtfmDiamIn , PtfmCDIn , &
!rm                                RdtnTMaxIn , RdtnDTIn  , NumLinesIn , LineModIn, &
!rm                                LAnchxiIn  , LAnchyiIn , LAnchziIn  , LFairxtIn, &
!rm                                LFairytIn  , LFairztIn , LUnstrLenIn, LDiam    , &
!rm                                LMassDen   , LEAStffIn , LSeabedCDIn, LTenTolIn, &
!rm                                LineNodesIn, LSNodesIn , X0                        )

!   SUBROUTINE InitFltngPtfmLd ( WAMITFile  , PtfmVol0In, PtfmDiamIn , PtfmCDIn , &
!                                RdtnTMaxIn , RdtnDTIn  , NumLinesIn , LineModIn, &
!                                LAnchxiIn  , LAnchyiIn , LAnchziIn  , LFairxtIn, &
!                                LFairytIn  , LFairztIn , LUnstrLenIn, LDiam    , &
!                                LMassDen   , LEAStffIn , LSeabedCDIn, LTenTolIn, &
!                                LineNodesIn, LSNodesIn , X0, DirRootIn,      )
   SUBROUTINE InitFltngPtfmLd ( FltPtfm_InitData , WaveDat, FP_Data , ErrStat )
!bjj end of proposed change

      ! This routine is used to initialize the variables used in the time
      ! domain hydrodynamic loading and mooring system dynamics routines
      ! for various floating platform concepts.


   USE                                    FFT_Module
   USE                                    Waves


   IMPLICIT                               NONE


      ! Passed Variables:
!BJJ START OF PROPOSED CHANGE V1.00.00A-BJJ
!   INTEGER,    INTENT(IN )      :: LineNodesIn                                     ! Number of nodes per line where the mooring line position and tension can be output (-)
!   INTEGER,    INTENT(IN )      :: NumLinesIn                                      ! Number of mooring lines (-)
!
!   REAL(ReKi), INTENT(IN )      :: LAnchxiIn  (NumLinesIn)                         ! xi-coordinate of each anchor   in the inertial frame        coordinate system (meters)
!   REAL(ReKi), INTENT(IN )      :: LAnchyiIn  (NumLinesIn)                         ! yi-coordinate of each anchor   in the inertial frame        coordinate system (meters)
!   REAL(ReKi), INTENT(IN )      :: LAnchziIn  (NumLinesIn)                         ! zi-coordinate of each anchor   in the inertial frame        coordinate system (meters)
!   REAL(ReKi), INTENT(IN )      :: LDiam      (NumLinesIn)                         ! Effective diameter of each mooring line for calculation of the line buoyancy (meters)
!   REAL(ReKi), INTENT(IN )      :: LEAStffIn  (NumLinesIn)                         ! Extensional stiffness of each mooring line (N)
!   REAL(ReKi), INTENT(IN )      :: LFairxtIn  (NumLinesIn)                         ! xt-coordinate of each fairlead in the tower base / platform coordinate system (meters)
!   REAL(ReKi), INTENT(IN )      :: LFairytIn  (NumLinesIn)                         ! yt-coordinate of each fairlead in the tower base / platform coordinate system (meters)
!   REAL(ReKi), INTENT(IN )      :: LFairztIn  (NumLinesIn)                         ! zt-coordinate of each fairlead in the tower base / platform coordinate system (meters)
!   REAL(ReKi), INTENT(IN )      :: LMassDen   (NumLinesIn)                         ! Mass density of each mooring line (kg/m)
!   REAL(ReKi), INTENT(IN )      :: LSeabedCDIn(NumLinesIn)                         ! Coefficient of seabed static friction drag of each mooring line (a negative value indicates no seabed) (-)
!   REAL(ReKi), INTENT(IN )      :: LSNodesIn  (NumLinesIn,LineNodesIn)             ! Unstretched arc distance along mooring line from anchor to each node where the line position and tension can be output (meters)
!   REAL(ReKi), INTENT(IN )      :: LTenTolIn  (NumLinesIn)                         ! Convergence tolerance within Newton-Raphson iteration of each mooring line specified as a fraction of tension (-)
!   REAL(ReKi), INTENT(IN )      :: LUnstrLenIn(NumLinesIn)                         ! Unstretched length of each mooring line (meters)
!   REAL(ReKi), INTENT(IN )      :: PtfmCDIn                                        ! Effective platform normalized hydrodynamic viscous drag coefficient in calculation of viscous drag term from Morison's equation (-)
!   REAL(ReKi), INTENT(IN )      :: PtfmDiamIn                                      ! Effective platform diameter in calculation of viscous drag term from Morison's equation (meters)
!   REAL(ReKi), INTENT(IN )      :: PtfmVol0In                                      ! Displaced volume of water when the platform is in its undisplaced position (m^3)
!   REAL(ReKi), INTENT(IN )      :: RdtnDTIn                                        ! Time step for wave radiation kernel calculations (sec)
!   REAL(ReKi), INTENT(IN )      :: RdtnTMaxIn                                      ! Analysis time for wave radiation kernel calculations; the actual analysis time may be larger than this value in order for the maintain an effecient (co)sine transform (sec)
!   REAL(ReKi), INTENT(IN )      :: X0         (6)                                  ! The 3 components of the initial translational displacement (in m) of the platform reference and the 3 components of the initial rotational displacement (in rad) of the platform relative to the inertial frame
!
!   INTEGER,    INTENT(IN )      :: LineModIn                                       ! Mooring line model switch {0: none, 1: standard quasi-static, 2: user-defined from routine UserLine} (switch)
!
!   CHARACTER(1024), INTENT(IN )   :: WAMITFile                                       ! Root name of WAMIT output files containing the linear, nondimensionalized, hydrostatic restoring matrix (.hst extension), frequency-dependent hydrodynamic added mass matrix and damping matrix (.1 extension), and frequency- and direction-dependent wave excitation force vector per unit wave amplitude (.3 extension)

   TYPE(FltPtfm_InitDataType), INTENT(IN) :: FltPtfm_InitData                      ! data to initialize the floating platform module
   TYPE(Waves_DataType),       INTENT(IN) :: WaveDat
   TYPE(FltPtfm_DataType),     INTENT(OUT):: FP_Data                               ! data for the Floating platform module
   INTEGER,                    INTENT(OUT):: ErrStat                               ! a non-zero value indicates an error has occurred
!BJJ END OF PROPOSED CHANGE

      ! Local Variables:

   COMPLEX(ReKi), ALLOCATABLE             :: HdroExctn (:,:,:)                    ! Frequency- and direction-dependent complex hydrodynamic wave excitation force per unit wave amplitude vector (kg/s^2, kg-m/s^2)
!jbj: start of proposed change v1.00.00b-jbj
!rm   COMPLEX(ReKi), ALLOCATABLE             :: WaveExctnC(:,:)                      ! Fourier transform of the instantaneous value of the total excitation force on the support platfrom from incident waves (N-s, N-m-s)
   COMPLEX(ReKi), ALLOCATABLE             :: WaveExctnC(:,:)                      ! Discrete Fourier transform of the instantaneous value of the total excitation force on the support platfrom from incident waves (N, N-m)
!jbj: end of proposed change v1.00.00b-jbj
   COMPLEX(ReKi), ALLOCATABLE             :: X_Diffrctn(:,:)                      ! Frequency-dependent complex hydrodynamic wave excitation force per unit wave amplitude vector at the chosen wave heading direction, WaveDir (kg/s^2, kg-m/s^2)

   REAL(ReKi)                             :: DffrctDim (6)                        ! Matrix used to redimensionalize WAMIT hydrodynamic wave excitation force  output (kg/s^2, kg-m/s^2            )
   REAL(ReKi), ALLOCATABLE                :: HdroAddMs (:,:)                      ! The upper-triangular portion (diagonal and above) of the frequency-dependent hydrodynamic added mass matrix from the radiation problem (kg  , kg-m  , kg-m^2  )
   REAL(ReKi), ALLOCATABLE                :: HdroDmpng (:,:)                      ! The upper-triangular portion (diagonal and above) of the frequency-dependent hydrodynamic damping    matrix from the radiation problem (kg/s, kg-m/s, kg-m^2/s)
   REAL(ReKi), ALLOCATABLE                :: HdroFreq  (:)                        ! Frequency components inherent in the hydrodynamic added mass matrix, hydrodynamic daming matrix, and complex wave excitation force per unit wave amplitude vector (rad/s)
   REAL(ReKi), ALLOCATABLE                :: HdroWvDir (:)                        ! Incident wave propagation heading direction components inherent in the complex wave excitation force per unit wave amplitude vector (degrees)
   REAL(ReKi)                             :: HighFreq    = 0.0                    ! The highest frequency component in the WAMIT file, not counting infinity.
   REAL(ReKi)                             :: Krnl_Fact                            ! Factor used to scale the magnitude of the RdtnKnrl  as required by the discrete time (co)sine transform (-)
!bjj rm v1.00.00a-bjj:   REAL(ReKi)                             :: Lamda0                               ! Catenary parameter used to generate the initial guesses of the horizontal and vertical tensions at the fairlead for the Newton-Raphson iteration (-)
!bjj rm v1.00.00a-bjj:   REAL(ReKi)                             :: LFairxi                              ! xi-coordinate of the current fairlead in the inertial frame coordinate system (meters)
!bjj rm v1.00.00a-bjj:   REAL(ReKi)                             :: LFairyi                              ! yi-coordinate of the current fairlead in the inertial frame coordinate system (meters)
!bjj rm v1.00.00a-bjj:   REAL(ReKi)                             :: LFairzi                              ! zi-coordinate of the current fairlead in the inertial frame coordinate system (meters)
   REAL(ReKi)                             :: Omega                                ! Wave frequency (rad/s)
   REAL(ReKi)                             :: PrvDir                               ! The value of TmpDir from the previous line (degrees)
   REAL(ReKi)                             :: PrvPer                               ! The value of TmpPer from the previous line (sec    )
!bjj start of proposed change v1.00.00a-bjj
!we'll rename this WAMITULEN and make it an input parameter
!rm   REAL(ReKi), PARAMETER                  :: PtfmULEN    = 1.0                    ! Characteristic body length scale used to redimensionalize WAMIT output (meters) !JASON: BECAUSE I HAVE FIXED THIS TO UNITY, THE WAMIT .GDF FILE MUST HAVE ULEN SET TO 1.0.  SHOULD WE INSTEAD MAKE PtfmULEN AN ACTUAL INPUT TO THE PROGRAM WITHIN PtfmFile???
!bjj end of proposed change v1.00.00a-bjj
   REAL(ReKi)                             :: RdtnDim   (6,6)                      ! Matrix used to redimensionalize WAMIT hydrodynamic added mass and damping output (kg    , kg-m    , kg-m^2    )
   REAL(ReKi)                             :: RdtnDOmega                           ! Frequency step for wave radiation kernel calculations (rad/s)
   REAL(ReKi)                             :: RdtnOmegaMax                         ! Maximum frequency used in the (co)sine transform to fine the radiation impulse response functions (rad/s)
   REAL(ReKi), ALLOCATABLE                :: RdtnTime  (:)                        ! Simulation times at which the instantaneous values of the wave radiation kernel are determined (sec)
   REAL(ReKi)                             :: SttcDim   (6,6)                      ! Matrix used to redimensionalize WAMIT hydrostatic  restoring              output (kg/s^2, kg-m/s^2, kg-m^2/s^2)
   REAL(ReKi)                             :: TmpData1                             ! A temporary           value  read in from a WAMIT file (-      )
   REAL(ReKi)                             :: TmpData2                             ! A temporary           value  read in from a WAMIT file (-      )
   REAL(ReKi)                             :: TmpDir                               ! A temporary direction        read in from a WAMIT file (degrees)
   REAL(ReKi)                             :: TmpIm                                ! A temporary imaginary value  read in from a WAMIT file (-      ) - stored as a REAL value
   REAL(ReKi)                             :: TmpPer                               ! A temporary period           read in from a WAMIT file (sec    )
   REAL(ReKi)                             :: TmpRe                                ! A temporary real      value  read in from a WAMIT file (-      )
!bjj rm v1.00.00a-bjj:   REAL(ReKi)                             :: TransMat0 (3,3)                      ! Transformation matrix from the inertial frame to the initial tower base / platform coordinate system (-)
   REAL(ReKi), ALLOCATABLE                :: WAMITFreq (:)                        ! Frequency      components as ordered in the WAMIT output files (rad/s  )
   REAL(ReKi), ALLOCATABLE                :: WAMITPer  (:)                        ! Period         components as ordered in the WAMIT output files (sec    )
   REAL(ReKi), ALLOCATABLE                :: WAMITWvDir(:)                        ! Wave direction components as ordered in the WAMIT output files (degrees)
!bjj rm v1.00.00a-bjj:   REAL(ReKi)                             :: XF                                   ! Horizontal distance between anchor and fairlead of the current mooring line (meters)
!bjj rm v1.00.00a-bjj:   REAL(ReKi)                             :: XF2                                  ! = XF*XF
!bjj rm v1.00.00a-bjj:   REAL(ReKi)                             :: ZF                                   ! Vertical   distance between anchor and fairlead of the current mooring line (meters)
!bjj rm v1.00.00a-bjj:   REAL(ReKi)                             :: ZF2                                  ! = ZF*ZF

   INTEGER                                :: I                                    ! Generic index
   INTEGER                                :: Indx                                 ! Cycles through the upper-triangular portion (diagonal and above) of the frequency-dependent hydrodynamic added mass and damping matrices from the radiation problem
   INTEGER                                :: InsertInd                            ! The lowest sorted index whose associated frequency component is higher than the current frequency component -- this is to sort the frequency components from lowest to highest
   INTEGER                                :: J                                    ! Generic index
   INTEGER                                :: K                                    ! Generic index
!bjj start of proposed change
!b/c this gets SAVEd by default when initialized here, I'm going to initialize in the subroutine body instead   
!rm   INTEGER                                :: LastInd     = 1                      ! Index into the arrays saved from the last call as a starting point for this call
   INTEGER                                :: LastInd                              ! Index into the arrays saved from the last call as a starting point for this call
!bjj end of proposed change   
   INTEGER                                :: NInpFreq                             ! Number of input frequency components inherent in the hydrodynamic added mass matrix, hydrodynamic daming matrix, and complex wave excitation force per unit wave amplitude vector (-)
   INTEGER                                :: NInpWvDir                            ! Number of input incident wave propagation heading direction components inherent in the complex wave excitation force per unit wave amplitude vector (-)
   INTEGER                                :: NStepRdtn2                           ! ( NStepRdtn-1 )/2
   INTEGER,    ALLOCATABLE                :: SortFreqInd (:)                      ! The array of indices such that WAMITFreq (SortFreqInd (:)) is sorted from lowest to highest frequency (-)
   INTEGER,    ALLOCATABLE                :: SortWvDirInd(:)                      ! The array of indices such that WAMITWvDir(SortWvDirInd(:)) is sorted from lowest to highest agnle     (-)
   INTEGER                                :: Sttus                                ! Status returned by an attempted allocation or READ.
!bjj start of proposed change v1.00.00a-bjj   
   INTEGER,    ALLOCATABLE                :: User_LineNodes (:)                   ! Temp array to hold line nodes on each line
!bjj end of proposed change v1.00.00a-bjj   
   INTEGER                                :: UnW1        = 31                     ! I/O unit number for the WAMIT output file with the .1   extension; this file contains the linear, nondimensionalized, frequency-dependent solution to the radiation   problem.
   INTEGER                                :: UnW3        = 32                     ! I/O unit number for the WAMIT output file with the .3   extension; this file contains the linear, nondimensionalized, frequency-dependent solution to the diffraction problem.
   INTEGER                                :: UnWh        = 33                     ! I/O unit number for the WAMIT output file with the .hst extension; this file contains the linear, nondimensionalized hydrostatic restoring matrix.

!bjj: these with the initialization are SAVEd by default???
   LOGICAL                                :: FirstFreq                            ! When .TRUE., indicates we're still looping through the first frequency component.
   LOGICAL                                :: FirstPass                            ! When .TRUE., indicates we're on the first pass through a loop.
!bjj start of proposed change
!b/c this gets SAVEd by default when initialized here, I'm going to initialize in the subroutine body instead   
!rm   LOGICAL                                :: InfFreq     = .FALSE.                ! When .TRUE., indicates that the infinite-frequency limit of added mass is contained within the WAMIT output files.
   LOGICAL                                :: InfFreq                              ! When .TRUE., indicates that the infinite-frequency limit of added mass is contained within the WAMIT output files.
!bjj end of proposed change   
   LOGICAL                                :: NewPer                               ! When .TRUE., indicates that the period has just changed.
!bjj start of proposed change
!b/c this gets SAVEd by default when initialized here, I'm going to initialize in the subroutine body instead   
!rm   LOGICAL                                :: RdtnFrmAM   = .FALSE.                ! Determine the wave radiation kernel from the frequency-dependent hydrodynamic added mass matrix? (.TRUE = yes, .FALSE. = determine the wave radiation kernel from the frequency-dependent hydrodynamic damping matrix) !JASON: SHOULD YOU MAKE THIS AN INPUT???<--JASON: IT IS NOT WISE TO COMPUTE THE RADIATION KERNEL FROM THE FREQUENCY-DEPENDENT ADDED MASS MATRIX, UNLESS A CORRECTION IS APPLIED.  THIS IS DESCRIBED IN THE WAMIT USER'S GUIDE!!!!
!rm   LOGICAL                                :: ZeroFreq    = .FALSE.                ! When .TRUE., indicates that the zero    -frequency limit of added mass is contained within the WAMIT output files.
   LOGICAL                                :: RdtnFrmAM                            ! Determine the wave radiation kernel from the frequency-dependent hydrodynamic added mass matrix? (.TRUE = yes, .FALSE. = determine the wave radiation kernel from the frequency-dependent hydrodynamic damping matrix) !JASON: SHOULD YOU MAKE THIS AN INPUT???<--JASON: IT IS NOT WISE TO COMPUTE THE RADIATION KERNEL FROM THE FREQUENCY-DEPENDENT ADDED MASS MATRIX, UNLESS A CORRECTION IS APPLIED.  THIS IS DESCRIBED IN THE WAMIT USER'S GUIDE!!!!
   LOGICAL                                :: ZeroFreq                             ! When .TRUE., indicates that the zero    -frequency limit of added mass is contained within the WAMIT output files.
!bjj end of proposed change   

   CHARACTER(1024)                        :: Line                                 ! String to temporarily hold the value of a line within a WAMIT output file.

!bjj start of proposed change v1.00.00a-bjj
!   CHARACTER(*), INTENT(IN)     :: DirRootIn
   TYPE(FFT_DataType)                     :: FFT_Data                             ! the instance of the FFT module we're using

      ! Initialize data
      
   ErrStat   = 0
   LastInd   = 1
   InfFreq   = .FALSE.
   RdtnFrmAM = .FALSE.
   ZeroFreq  = .FALSE.
!bjj end of proposed change v1.00.00a-bjj


      ! Save these values for future use:

!bjj start of proposed change v1.00.00a-bjj
!   PtfmVol0     = PtfmVol0In
!   PtfmDiam     = PtfmDiamIn
!   PtfmCD       = PtfmCDIn
!   RdtnDT       = RdtnDTIn
!   IF ( RdtnTMaxIn == 0.0 )  THEN   ! .TRUE. when we don't want to model wave radiation damping; set RdtnTMax to some minimum value greater than zero to avoid an error in the calculations below.
!      RdtnTMax  = RdtnDTIn
!      UseRdtn   = .FALSE.
!   ELSE                             ! We will be modeling wave radiation damping.
!      RdtnTMax  = RdtnTMaxIn
!      UseRdtn   = .TRUE.
!   ENDIF

   FP_Data%PtfmVol0     = FltPtfm_InitData%PtfmVol0
   FP_Data%PtfmDiam     = FltPtfm_InitData%PtfmDiam
   FP_Data%PtfmCD       = FltPtfm_InitData%PtfmCD
   FP_Data%RdtnDT       = FltPtfm_InitData%RdtnDT
   IF ( FltPtfm_InitData%RdtnTMax == 0.0 )  THEN   ! .TRUE. when we don't want to model wave radiation damping; set RdtnTMax to some minimum value greater than zero to avoid an error in the calculations below.
      FP_Data%RdtnTMax  = FltPtfm_InitData%RdtnDT
      FP_Data%UseRdtn   = .FALSE.
   ELSE                             ! We will be modeling wave radiation damping.
      FP_Data%RdtnTMax  = FltPtfm_InitData%RdtnTMax
      FP_Data%UseRdtn   = .TRUE.
   ENDIF

!bjj end of proposed change v1.00.00a-bjj


   RdtnOmegaMax = Pi/FP_Data%RdtnDT




      ! ProgAbort if the wave elevation has not been computed yet:
!bjj start of proposed change v1.00.00a-bjj
!rm   IF ( .NOT. ALLOCATED ( WaveElev ) )  THEN
   IF ( .NOT. Waves_IsAllocated( WaveDat, 'WaveElev  ', ErrStat ) ) THEN
!bjj end of proposed change v1.00.00a-bjj
      CALL ProgAbort ( ' Routine InitWaves() must be called before routine InitFltngPtfmLd().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   END IF

      ! Initialize the variables associated with the mooring system:

   FP_Data%LineMod      = FltPtfm_InitData%LineMod
!bjj start of proposed change V1.00.00a-bjj
   
   IF ( FP_Data%LineMod == 2 )  THEN  ! .TRUE if we have user-defined mooring lines.
   
         ! call the UserLineModule's initialization routine to get NumLines and LineNodes
         
!      CALL UserLine_Init( FltPtfm_InitData%DirRoot, FP_Data%UserLine_Data, ErrStat )
      CALL UserLine_Init( FltPtfm_InitData%DirRoot, FP_Data%NumLines, User_LineNodes, &
                                                    FP_Data%UserLine_Data, ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgAbort ( ' Error in call to UserLine_Init().', TrapErrors = .TRUE.)
         RETURN
      END IF
      
!      FP_Data%NumLines = FP_Data%UserLine_Data%NumLines
      
      
      ALLOCATE ( FP_Data%MooringLine( FP_Data%NumLines ), STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the MooringLine array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF
           
      
      DO I = 1,FP_Data%NumLines

!         FP_Data%MooringLine(I)%LineNodes = FP_Data%UserLine_Data%MooringLine(I)%LineNodes
         FP_Data%MooringLine(I)%LineNodes = User_LineNodes(I)

         ALLOCATE( FP_Data%MooringLine(I)%LNodesPi( FP_Data%MooringLine(I)%LineNodes, 3), STAT = ErrStat )
         IF ( ErrStat /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the LNodesPi array.', TrapErrors = .TRUE.)
            RETURN
         ENDIF

         FP_Data%MooringLine(I)%LNodesPi = 0.0_ReKi


         ALLOCATE( FP_Data%MooringLine(I)%LNodesTe( FP_Data%MooringLine(I)%LineNodes   ), STAT = ErrStat )
         IF ( ErrStat /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the LNodesTe array.', TrapErrors = .TRUE.)
            RETURN
         ENDIF

         FP_Data%MooringLine(I)%LNodesTe = 0.0_ReKi

      END DO
      
      
   ELSE  ! .TRUE if we have standard quasi-static mooring lines

      FP_Data%NumLines = FltPtfm_InitData%NumLines
      
      ALLOCATE ( FP_Data%MooringLine( FP_Data%NumLines ), STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort ( ' Error allocating memory for the MooringLine array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

!   END IF 
!rm   LineNodes    = FltPtfm_InitData%LineNodes
!RM   FP_Data%NumLines     = FltPtfm_InitData%NumLines
!REMOVE THESE NOW, TOO:
!RM   ALLOCATE ( LAnchHTe (NumLines            ) , STAT=ErrStat )
!RM   IF ( ErrStat /= 0 )  THEN
!RM      CALL ProgAbort ( ' Error allocating memory for the LAnchHTe array.', TrapErrors = .TRUE.)
!RM      RETURN
!RM   ENDIF
!RM
!RM   ALLOCATE ( LAnchVTe (NumLines            ) , STAT=ErrStat )
!RM   IF ( ErrStat /= 0 )  THEN
!RM      CALL ProgAbort ( ' Error allocating memory for the LAnchVTe array.', TrapErrors = .TRUE.)
!RM      RETURN
!RM   ENDIF
!RM
!RM   ALLOCATE ( LFairHTe (NumLines            ) , STAT=ErrStat )
!RM   IF ( ErrStat /= 0 )  THEN
!RM      CALL ProgAbort ( ' Error allocating memory for the LFairHTe array.', TrapErrors = .TRUE.)
!RM      RETURN
!RM   ENDIF
!RM
!RM   ALLOCATE ( LFairVTe (NumLines            ) , STAT=ErrStat )
!RM   IF ( ErrStat /= 0 )  THEN
!RM      CALL ProgAbort ( ' Error allocating memory for the LFairVTe array.', TrapErrors = .TRUE.)
!RM      RETURN
!RM   ENDIF
!RM
!RM   ALLOCATE ( LNodesPi (NumLines,LineNodes,3) , STAT=ErrStat )
!RM   IF ( ErrStat /= 0 )  THEN
!RM      CALL ProgAbort ( ' Error allocating memory for the LNodesPi array.', TrapErrors = .TRUE.)
!RM      RETURN
!RM   ENDIF
!RM
!RM   ALLOCATE ( LNodesTe (NumLines,LineNodes  ) , STAT=ErrStat )
!RM   IF ( ErrStat /= 0 )  THEN
!RM      CALL ProgAbort ( ' Error allocating memory for the LNodesTe array.', TrapErrors = .TRUE.)
!RM      RETURN
!RM   ENDIF
!RM
!RM   LAnchHTe(:    ) = 0.0
!RM   LAnchVTe(:    ) = 0.0
!RM   LFairHTe(:    ) = 0.0
!RM   LFairVTe(:    ) = 0.0
!RM   LNodesPi(:,:,:) = 0.0
!RM   LNodesTe(:,:  ) = 0.0
!
!   ALLOCATE ( FP_Data%MooringLine( FP_Data%NumLines ), STAT=ErrStat )
!   IF ( ErrStat /= 0 )  THEN
!      CALL ProgAbort ( ' Error allocating memory for the MooringLine array.', TrapErrors = .TRUE.)
!      RETURN
!   ENDIF

!rm   IF ( FP_Data%LineMod == 1 )  THEN  ! .TRUE if we have standard quasi-static mooring lines.
!rm
!rm      ALLOCATE ( LAnchxi  (NumLines          ) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LAnchxi array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LAnchyi  (NumLines          ) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LAnchyi array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LAnchzi  (NumLines          ) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LAnchzi array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LFairxt  (NumLines          ) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LFairxt array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LFairyt  (NumLines          ) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LFairyt array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LFairzt  (NumLines          ) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LFairzt array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LUnstrLen(NumLines          ) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LUnstrLen array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LEAStff  (NumLines          ) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LEAStff array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LSeabedCD(NumLines          ) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LSeabedCD array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LTenTol  (NumLines          ) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LTenTol array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LFldWght (NumLines          ) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LFldWght array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LSNodes  (NumLines,LineNodes) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LSNodes array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LNodesX  (         LineNodes) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LNodesX array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm      ALLOCATE ( LNodesZ  (         LineNodes) , STAT=ErrStat )
!rm      IF ( ErrStat /= 0 )  THEN
!rm         CALL ProgAbort ( ' Error allocating memory for the LNodesZ array.', TrapErrors = .TRUE.)
!rm         RETURN
!rm      ENDIF
!rm
!rm
!rm
!rm      LAnchxi  (:  ) = FltPtfm_InitData%LAnchxi  (:  )
!rm      LAnchyi  (:  ) = FltPtfm_InitData%LAnchyi  (:  )
!rm      LAnchzi  (:  ) = FltPtfm_InitData%LAnchzi  (:  )
!rm      LFairxt  (:  ) = FltPtfm_InitData%LFairxt  (:  )
!rm      LFairyt  (:  ) = FltPtfm_InitData%LFairyt  (:  )
!rm      LFairzt  (:  ) = FltPtfm_InitData%LFairzt  (:  )
!rm      LUnstrLen(:  ) = FltPtfm_InitData%LUnstrLen(:  )
!rm      LEAStff  (:  ) = FltPtfm_InitData%LEAStff  (:  )
!rm      LSeabedCD(:  ) = FltPtfm_InitData%LSeabedCD(:  )
!rm      LTenTol  (:  ) = FltPtfm_InitData%LTenTol  (:  )
!rm      LSNodes  (:,:) = FltPtfm_InitData%LSNodes  (:,:)


      DO I = 1,FP_Data%NumLines

         FP_Data%MooringLine(I)%LAnchxi   = FltPtfm_InitData%MooringLine(I)%LAnchxi
         FP_Data%MooringLine(I)%LAnchyi   = FltPtfm_InitData%MooringLine(I)%LAnchyi
         FP_Data%MooringLine(I)%LAnchzi   = FltPtfm_InitData%MooringLine(I)%LAnchzi
         FP_Data%MooringLine(I)%LFairxt   = FltPtfm_InitData%MooringLine(I)%LFairxt
         FP_Data%MooringLine(I)%LFairyt   = FltPtfm_InitData%MooringLine(I)%LFairyt
         FP_Data%MooringLine(I)%LFairzt   = FltPtfm_InitData%MooringLine(I)%LFairzt
         FP_Data%MooringLine(I)%LUnstrLen = FltPtfm_InitData%MooringLine(I)%LUnstrLen
         FP_Data%MooringLine(I)%LEAStff   = FltPtfm_InitData%MooringLine(I)%LEAStff
         FP_Data%MooringLine(I)%LSeabedCD = FltPtfm_InitData%MooringLine(I)%LSeabedCD
         FP_Data%MooringLine(I)%LTenTol   = FltPtfm_InitData%MooringLine(I)%LTenTol

         FP_Data%MooringLine(I)%LineNodes = FltPtfm_InitData%MooringLine(I)%LineNodes

         ALLOCATE( FP_Data%MooringLine(I)%LSNodes(  FP_Data%MooringLine(I)%LineNodes   ), STAT = ErrStat )
         IF ( ErrStat /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the LSNodes array.', TrapErrors = .TRUE.)
            RETURN
         ENDIF

         IF ( FP_Data%MooringLine(I)%LineNodes > 0 ) THEN
            FP_Data%MooringLine(I)%LSNodes = FltPtfm_InitData%MooringLine(I)%LSNodes
         END IF


         ALLOCATE( FP_Data%MooringLine(I)%LNodesX(  FP_Data%MooringLine(I)%LineNodes   ), STAT = ErrStat )
         IF ( ErrStat /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the LNodesX array.', TrapErrors = .TRUE.)
            RETURN
         ENDIF


         ALLOCATE( FP_Data%MooringLine(I)%LNodesZ(  FP_Data%MooringLine(I)%LineNodes   ), STAT = ErrStat )
         IF ( ErrStat /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the LNodesZ array.', TrapErrors = .TRUE.)
            RETURN
         ENDIF


         ALLOCATE( FP_Data%MooringLine(I)%LNodesPi( FP_Data%MooringLine(I)%LineNodes, 3), STAT = ErrStat )
         IF ( ErrStat /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the LNodesPi array.', TrapErrors = .TRUE.)
            RETURN
         ENDIF

         FP_Data%MooringLine(I)%LNodesPi = 0.0_ReKi

         ALLOCATE( FP_Data%MooringLine(I)%LNodesTe( FP_Data%MooringLine(I)%LineNodes   ), STAT = ErrStat )
         IF ( ErrStat /= 0 )  THEN
            CALL ProgAbort ( ' Error allocating memory for the LNodesTe array.', TrapErrors = .TRUE.)
            RETURN
         ENDIF

         FP_Data%MooringLine(I)%LNodesTe = 0.0_ReKi

      END DO
!bjj END of proposed change v1.00.00a-bjj


!bjj start of proposed change v1.00.00a-bjj
!bjj: move this initialization to the load calculation so that we don't need an initial displacement for our guess:
!rm      ! Get the transformation matrix, TransMat0, from the inertial frame to the
!rm      !   initial tower base / platform coordinate system:
!rm
!rm      CALL SmllRotTrans ( 'platform displacement', FltPtfm_InitData%X0(4), FltPtfm_InitData%X0(5), &
!rm                                                   FltPtfm_InitData%X0(6), TransMat0 )
!bjj end of proposed change v1.00.00a-bjj


      DO I = 1,FP_Data%NumLines ! Loop through all mooring lines


      ! Compute the weight of each mooring line in fluid per unit length based on
      !   their mass density and effective diameter, water density, and gravity:
      ! NOTE: The buoyancy is calculated assuming that the entire length of the
      !       mooring line is submerged in the water.

!bjj start of proposed change v1.00.00a-bjj
!rm         LFldWght(I) = ( LMassDen(I) - WtrDens*PiOvr4*LDiam(I)*LDiam(I) )*Gravity
         FP_Data%MooringLine(I)%LFldWght = ( FltPtfm_InitData%MooringLine(I)%LMassDen - &
                                             WaveDat%WtrDens*PiOvr4*(FltPtfm_InitData%MooringLine(I)%LDiam**2) )*WaveDat%Gravity
!bjj end of proposed change

!bjj start of proposed change v1.00.00a-bjj
!rm      ! Transform the fairlead location from the initial platform to the inertial
!rm      !    frame coordinate system:
!rm      ! NOTE: TransMat0^T = TransMat0^-1 where ^T = matrix transpose and ^-1 =
!rm      !       matrix inverse.
!rm
!rm         LFairxi = FltPtfm_InitData%X0(1) + TransMat0(1,1)*FP_Data%MooringLine(I)%LFairxt &
!rm                                          + TransMat0(2,1)*FP_Data%MooringLine(I)%LFairyt &
!rm                                          + TransMat0(3,1)*FP_Data%MooringLine(I)%LFairzt
!rm         LFairyi = FltPtfm_InitData%X0(2) + TransMat0(1,2)*FP_Data%MooringLine(I)%LFairxt &
!rm                                          + TransMat0(2,2)*FP_Data%MooringLine(I)%LFairyt &
!rm                                          + TransMat0(3,2)*FP_Data%MooringLine(I)%LFairzt
!rm         LFairzi = FltPtfm_InitData%X0(3) + TransMat0(1,3)*FP_Data%MooringLine(I)%LFairxt &
!rm                                          + TransMat0(2,3)*FP_Data%MooringLine(I)%LFairyt &
!rm                                          + TransMat0(3,3)*FP_Data%MooringLine(I)%LFairzt
!rm
!rm
!rm      ! Transform the fairlead location from the inertial frame coordinate system
!rm      !   to the local coordinate system of the current line (this coordinate
!rm      !   system lies at the current anchor, Z being vertical, and X directed from
!rm      !   the current anchor to the current fairlead):
!rm
!rm!bjj: note that i changed the exponents from 2.0 to 2 so that the compiler knows that the exponent is an integer and can perform the computation easier!
!rm         XF      = SQRT( ( LFairxi - FP_Data%MooringLine(I)%LAnchxi )**2 + ( LFairyi - FP_Data%MooringLine(I)%LAnchyi )**2 )
!rm         ZF      =         LFairzi - FP_Data%MooringLine(I)%LAnchzi
!rm
!rm         XF2     = XF*XF
!rm         ZF2     = ZF*ZF
!rm
!rm
!rm      ! Generate the initial guess values for the horizontal and vertical tensions
!rm      !   at the fairlead in the Newton-Raphson iteration for the catenary mooring
!rm      !   line solution.  Use starting values documented in: Peyrot, Alain H. and
!rm      !   Goulois, A. M., "Analysis Of Cable Structures," Computers & Structures,
!rm      !   Vol. 10, 1979, pp. 805-813:
!rm
!rm         IF     ( XF                               == 0.0               )  THEN ! .TRUE. if the current mooring line is exactly vertical
!rm            Lamda0 = 1.0E+06
!rm         ELSEIF ( FP_Data%MooringLine(I)%LUnstrLen <= SQRT( XF2 + ZF2 ) )  THEN ! .TRUE. if the current mooring line is taut
!rm            Lamda0 = 0.2
!rm         ELSE                                                                   ! The current mooring line must be slack and not vertical
!rm            Lamda0 = SQRT( 3.0*( ( FP_Data%MooringLine(I)%LUnstrLen**2 - ZF2 )/XF2 - 1.0 ) )
!rm         ENDIF
!rm
!rm         FP_Data%MooringLine(I)%LFairHTe = ABS( 0.5*FP_Data%MooringLine(I)%LFldWght*  XF/     Lamda0    )
!rm         FP_Data%MooringLine(I)%LFairVTe =      0.5*FP_Data%MooringLine(I)%LFldWght*( ZF/TANH(Lamda0) + &
!rm                                                    FP_Data%MooringLine(I)%LUnstrLen                    )
!rm
!rm
!bjj end of proposed change v1.00.00a-bjj
      ENDDO             ! I - All mooring lines

   ENDIF




      ! Tell our nice users what is about to happen that may take a while:

   CALL WrScr ( ' Reading in WAMIT output with root name "'//TRIM(FltPtfm_InitData%WAMITFile)//'".' )



      ! Let's set up the matrices used to redimensionalize the hydrodynamic data
      !   from WAMIT; all these matrices are symmetric and need to be used with
      !   element-by-element multiplication, instead of matrix-by-matrix
      !   multiplication:

!BJJ: I replaced PtfmULEN with FltPtfm_InitData%WAMITULEN in the following 8 equations:
   SttcDim(1,1) = WaveDat%RhoXg  *FltPtfm_InitData%WAMITULEN**2  ! Force-translation
   SttcDim(1,4) = WaveDat%RhoXg  *FltPtfm_InitData%WAMITULEN**3  ! Force-rotation/Moment-translation - Hydrostatic restoring
   SttcDim(4,4) = WaveDat%RhoXg  *FltPtfm_InitData%WAMITULEN**4  ! Moment-rotation

   RdtnDim(1,1) = WaveDat%WtrDens*FltPtfm_InitData%WAMITULEN**3  ! Force-translation
   RdtnDim(1,4) = WaveDat%WtrDens*FltPtfm_InitData%WAMITULEN**4  ! Force-rotation/Moment-translation - Hydrodynamic added mass and damping
   RdtnDim(4,4) = WaveDat%WtrDens*FltPtfm_InitData%WAMITULEN**5  ! Moment-rotation

   DffrctDim(1) = WaveDat%RhoXg  *FltPtfm_InitData%WAMITULEN**2  ! Force-translation - Hydrodynamic wave excitation force
   DffrctDim(4) = WaveDat%RhoXg  *FltPtfm_InitData%WAMITULEN**3  ! Moment-rotation

   DO I = 1,3     ! Loop through all force-translation elements (rows)

      DO J = 1,3  ! Loop through all force-translation elements (columns)

         SttcDim(I,J) = SttcDim(1,1)

         RdtnDim(I,J) = RdtnDim(1,1)

      ENDDO       ! J - All force-translation elements (columns)

      DffrctDim (I  ) = DffrctDim(1)

   ENDDO          ! I - All force-translation elements (rows)

   DO I = 1,3     ! Loop through all force-rotation/moment-translation elements (rows/columns)

      DO J = 4,6  ! Loop through all force-rotation/moment-translation elements (columns/rows)

         SttcDim(I,J) = SttcDim(1,4)
         SttcDim(J,I) = SttcDim(1,4)

         RdtnDim(I,J) = RdtnDim(1,4)
         RdtnDim(J,I) = RdtnDim(1,4)

      ENDDO       ! J - All force-rotation/moment-translation elements (rows/columns)

   ENDDO          ! I - All force-rotation/moment-translation elements (columns/rows)

   DO I = 4,6     ! Loop through all moment-rotation elements (rows)

      DO J = 4,6  ! Loop through all moment-rotation elements (columns)

         SttcDim(I,J) = SttcDim(4,4)

         RdtnDim(I,J) = RdtnDim(4,4)

      ENDDO       ! J - All moment-rotation elements (columns)

      DffrctDim (I  ) = DffrctDim(4)

   ENDDO          ! I - All moment-rotation elements (rows)




      ! Let's read in and redimensionalize the hydrodynamic data from the WAMIT
      !   output files:



      ! Linear restoring from the hydrostatics problem:

   CALL OpenFInpFile ( UnWh, TRIM(FltPtfm_InitData%WAMITFile)//'.hst', ErrStat )  ! Open file.
   IF ( ErrStat /= 0 ) RETURN

   FP_Data%HdroSttc (:,:) = 0.0 ! Initialize to zero

   DO    ! Loop through all rows in the file


      READ (UnWh,*,IOSTAT=Sttus)  I, J, TmpData1   ! Read in the row index, column index, and nondimensional data from the WAMIT file

      IF ( Sttus == 0 )  THEN                ! .TRUE. when data is read in successfully

         FP_Data%HdroSttc (I,J) = TmpData1*SttcDim(I,J)    ! Redimensionalize the data and place it at the appropriate location within the array

      ELSE                                   ! We must have reached the end of the file, so stop reading in data

         EXIT

      ENDIF


   ENDDO ! End loop through all rows in the file

   CLOSE ( UnWh ) ! Close file.



      ! Linear, frequency-dependent hydrodynamic added mass and damping from the
      !   radiation problem:

   CALL OpenFInpFile ( UnW1, TRIM(FltPtfm_InitData%WAMITFile)//'.1', ErrStat   )  ! Open file.
   IF ( ErrStat /= 0 ) RETURN


      ! First find the number of input frequency components inherent in the
      !   hydrodynamic added mass matrix, hydrodynamic daming matrix, and complex
      !   wave excitation force per unit wave amplitude vector:

   NInpFreq  = 0        ! Initialize to zero
   PrvPer    = 0.0      ! Initialize to a don't care
   FirstPass = .TRUE.   ! Initialize to .TRUE. for the first pass

   DO    ! Loop through all rows in the file


      READ (UnW1,*,IOSTAT=Sttus)  TmpPer  ! Read in only the period from the WAMIT file

      IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully

         IF ( FirstPass .OR. ( TmpPer /= PrvPer ) )  THEN   ! .TRUE. if we are on the first pass or if the period currently read in is different than the previous period read in; thus we found a new frequency in the WAMIT file!
            NInpFreq  = NInpFreq + 1      ! Since we found a new frequency, count it in the total
            PrvPer    = TmpPer            ! Store the current period as the previous period for the next pass
            FirstPass = .FALSE.           ! Sorry, you can only have one first pass
         ENDIF

      ELSE                    ! We must have reached the end of the file, so stop reading in data  !bjj -- thiw isn't necessarially true....

         EXIT

      ENDIF


   ENDDO ! End loop through all rows in the file


   REWIND (UNIT=UnW1)   ! REWIND the file so we can read it in a second time.


   ! Now that we know how many frequencies there are, we can ALLOCATE the arrays
   !   to store the frequencies and frequency-dependent hydrodynamic added mass
   !   and damping matrices:

   ALLOCATE ( WAMITFreq  (NInpFreq   ) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the WAMITFreq array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

   ALLOCATE ( WAMITPer   (NInpFreq   ) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the WAMITPer array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

   ALLOCATE ( SortFreqInd(NInpFreq   ) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the SortFreqInd array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

   ALLOCATE ( HdroFreq   (NInpFreq   ) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the HdroFreq array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

   ALLOCATE ( HdroAddMs  (NInpFreq,21) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the HdroAddMs array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

   ALLOCATE ( HdroDmpng  (NInpFreq,21) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the HdroDmpng array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF


      ! Now find out how the frequencies are ordered in the file.  When we read in
      !   the added mass and damping matrices, we need to have them sorted by
      !   increasing frequency.  Thus, find the array of indices, SortFreqInd(),
      !   such that WAMITFreq(SortFreqInd(:)) is sorted from lowest to highest
      !   frequency:

   K         = 0        ! Initialize to zero
   PrvPer    = 0.0      ! Initialize to a don't care
   FirstPass = .TRUE.   ! Initialize to .TRUE. for the first pass

   DO    ! Loop through all rows in the file


      READ (UnW1,*,IOSTAT=Sttus)  TmpPer  ! Read in only the period from the WAMIT file

      IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully

         IF ( FirstPass .OR. ( TmpPer /= PrvPer ) )  THEN   ! .TRUE. if we are on the first pass or if the period currently read in is different than the previous period read in; thus we found a new frequency in the WAMIT file!

            K               = K + 1       ! This is current count of which frequency component we are on
            PrvPer          = TmpPer      ! Store the current period as the previous period for the next pass
            FirstPass       = .FALSE.     ! Sorry, you can only have one first pass

            WAMITPer    (K) = TmpPer         ! Store the periods                         in the order they appear in the WAMIT file
            IF (     TmpPer <  0.0 )  THEN   ! Periods less than zero in WAMIT represent infinite period = zero frequency
               WAMITFreq(K) = 0.0
               ZeroFreq     = .TRUE.
            ELSEIF ( TmpPer == 0.0 )  THEN   ! Periods equal to  zero in WAMIT represent infinite frequency
               WAMITFreq(K) = HUGE(TmpPer)   ! Use HUGE() to approximate infinity in the precision of ReKi
               InfFreq      = .TRUE.
            ELSE                             ! We must have positive, non-infinite frequency
               WAMITFreq(K) = TwoPi/TmpPer   ! Store the periods as frequencies in rad/s in the order they appear in the WAMIT file
               HighFreq     = MAX( HighFreq, WAMITFreq(K) ) ! Find the highest frequency (HighFreq) in the WAMIT output file, not counting infinity (even if the infinite frequency limit is in the file).
            ENDIF

            InsertInd       = K           ! Initialize as the K'th component
            DO I = 1,K-1   ! Loop throuh all previous frequencies
               IF ( ( WAMITFreq(I) > WAMITFreq(K) ) )  THEN ! .TRUE. if a previous frequency component is higher than the current frequency component
                  InsertInd      = MIN( InsertInd, SortFreqInd(I) )  ! Store the lowest sorted index whose associated frequency component is higher than the current frequency component
                  SortFreqInd(I) = SortFreqInd(I) + 1                ! Shift all of the sorted indices up by 1 whose associated frequency component is higher than the current frequency component
               ENDIF
            ENDDO          ! I - All previous frequencies
            SortFreqInd(K)  = InsertInd   ! Store the index such that WAMITFreq(SortFreqInd(:)) is sorted from lowest to highest frequency

         ENDIF

      ELSE                    ! We must have reached the end of the file, so stop reading in data

         EXIT

      ENDIF


   ENDDO ! End loop through all rows in the file


   REWIND (UNIT=UnW1)   ! REWIND the file so we can read it in a third time.  (This is getting ridiculous!)


      ! Now we can finally read in the frequency-dependent added mass and damping
      !   matrices; only store the upper-triangular portions (diagonal and above)
      !   of these matrices:

   K              = 0      ! Initialize to zero
   PrvPer         = 0.0    ! Initialize to a don't care
   FirstPass      = .TRUE. ! Initialize to .TRUE. for the first pass

   HdroAddMs(:,:) = 0.0    ! Initialize to zero
   HdroDmpng(:,:) = 0.0    ! Initialize to zero

   DO    ! Loop through all rows in the file


      READ (UnW1,'(A)',IOSTAT=Sttus)  Line   ! Read in the entire line

      IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully


         READ (Line,*)  TmpPer               ! Read in only the period from the WAMIT file


         IF ( FirstPass .OR. ( TmpPer /= PrvPer ) )  THEN   ! .TRUE. if we are on the first pass or if the period currently read in is different than the previous period read in; thus we found a new frequency in the WAMIT file!

            K              = K + 1           ! This is current count of which frequency component we are on
            PrvPer         = TmpPer          ! Store the current period as the previous period for the next pass
            FirstPass      = .FALSE.         ! Sorry, you can only have one first pass

            IF (     TmpPer <  0.0 )  THEN   ! Periods less than zero in WAMIT represent infinite period = zero frequency
               HdroFreq (SortFreqInd(K)) = 0.0
            ELSEIF ( TmpPer == 0.0 )  THEN   ! Periods equal to  zero in WAMIT represent infinite frequency; a value slightly larger than HighFreq is returned to approximate infinity while still maintaining an effective interpolation later on.
               HdroFreq (SortFreqInd(K)) = HighFreq*OnePlusEps ! Set the infinite frequency to a value slightly larger than HighFreq
            ELSE                             ! We must have positive, non-infinite frequency
               HdroFreq (SortFreqInd(K)) = TwoPi/TmpPer  ! Convert the period in seconds to a frequency in rad/s and store them sorted from lowest to highest
            ENDIF

         ENDIF


         IF ( TmpPer <= 0.0 )  THEN          ! .TRUE. if the current period is less than or equal to zero, which in WAMIT represents the zero and infinite frequency limits, respectively; in these cases, only the added mass matrix is computed and output by WAMIT (and based on hydrodynamic theory, the damping matrix is zero as initialized above)

            READ (Line,*)  TmpPer, I, J, TmpData1           ! Read in the period, row index, column index, and nondimensional data from the WAMIT file

            IF ( J >= I )  THEN  ! .TRUE. if we are on or above the diagonal
               Indx = 6*( I - 1 ) + J - ( I*( I - 1 ) )/2                                       ! Convert from row/column indices to an index in the format used to save only the upper-triangular portion of the matrix.  NOTE: ( I*( I - 1 ) )/2 = SUM(I,START=1,END=I-1).

               HdroAddMs(SortFreqInd(K),Indx) = TmpData1*RdtnDim(I,J)                           ! Redimensionalize the data and place it at the appropriate location within the array
            ENDIF

         ELSE                                ! We must have a positive, non-infinite frequency.

            READ (Line,*)  TmpPer, I, J, TmpData1, TmpData2 ! Read in the period, row index, column index, and nondimensional data from the WAMIT file

            IF ( J >= I )  THEN  ! .TRUE. if we are on or above the diagonal
               Indx = 6*( I - 1 ) + J - ( I*( I - 1 ) )/2                                       ! Convert from row/column indices to an index in the format used to save only the upper-triangular portion of the matrix.  NOTE: ( I*( I - 1 ) )/2 = SUM(I,START=1,END=I-1).

               HdroAddMs(SortFreqInd(K),Indx) = TmpData1*RdtnDim(I,J)                           ! Redimensionalize the data and place it at the appropriate location within the array
               HdroDmpng(SortFreqInd(K),Indx) = TmpData2*RdtnDim(I,J)*HdroFreq(SortFreqInd(K))  ! Redimensionalize the data and place it at the appropriate location within the array
            ENDIF

         ENDIF


      ELSE                    ! We must have reached the end of the file, so stop reading in data


         EXIT


      ENDIF


   ENDDO ! End loop through all rows in the file


   CLOSE ( UnW1 ) ! Close file.



      ! Linear, frequency- and direction-dependent complex hydrodynamic wave
      !   excitation force per unit wave amplitude vector from the diffraction
      !   problem:

   CALL OpenFInpFile ( UnW3, TRIM(FltPtfm_InitData%WAMITFile)//'.3', ErrStat   )  ! Open file.
   IF ( ErrStat /= 0 ) RETURN


      ! First find the number of input incident wave propagation heading direction
      !   components inherent in the complex wave excitation force per unit wave
      !   amplitude vector:

   NInpWvDir = 0        ! Initialize to zero
   PrvDir    = 0.0      ! Initialize to a don't care
   FirstPass = .TRUE.   ! Initialize to .TRUE. for the first pass

   DO    ! Loop through all rows in the file


      READ (UnW3,'(A)',IOSTAT=Sttus)  Line   ! Read in the entire line

      IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully


         READ (Line,*)  TmpPer, TmpDir ! Read in only the period and direction from the WAMIT file !bjj why don't we check IOSTAT here, too????


         IF ( FirstPass                           )  THEN   ! .TRUE. if we are on the first pass
            PrvPer = TmpPer            ! Store the current period    as the previous period    for the next pass
         ENDIF


         IF (                  TmpPer /= PrvPer   )  THEN   ! .TRUE.                                if the period    currently read in is different than the previous period    read in; thus we found a new period    in the WAMIT file, so stop reading in data
            EXIT
         ENDIF


         IF ( FirstPass .OR. ( TmpDir /= PrvDir ) )  THEN   ! .TRUE. if we are on the first pass or if the direction currently read in is different than the previous direction read in; thus we found a new direction in the WAMIT file!
            NInpWvDir = NInpWvDir + 1  ! Since we found a new direction, count it in the total
            PrvDir    = TmpDir         ! Store the current direction as the previous direction for the next pass
            FirstPass = .FALSE.        ! Sorry, you can only have one first pass
         ENDIF


      ELSE                    ! We must have reached the end of the file, so stop reading in data


         EXIT


      ENDIF


   ENDDO ! End loop through all rows in the file


   REWIND (UNIT=UnW3)   ! REWIND the file so we can read it in a second time.


   ! Now that we know how many directions there are, we can ALLOCATE the arrays to
   !   to store the directions and frequency- and direction-dependent complex wave
   !   excitation force per unit wave amplitude vector:

   ALLOCATE ( WAMITWvDir  (NInpWvDir           ) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the WAMITWvDir array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

   ALLOCATE ( SortWvDirInd(NInpWvDir           ) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the SortWvDirInd array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

   ALLOCATE ( HdroWvDir   (NInpWvDir           ) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the HdroWvDir array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

   ALLOCATE ( HdroExctn   (NInpFreq,NInpWvDir,6) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the HdroExctn array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF


      ! Now find out how the directions are ordered in the file.  When we read in
      !   the wave excitation force vector, we need to have them sorted by
      !   increasing angle.  Thus, find the array of indices, SortWvDirInd(),
      !   such that WAMITWvDir(SortWvDirInd(:)) is sorted from lowest to highest
      !   angle.  At the same time, make sure that the frequencies in the .3 file are
      !   ordered in the same way they are in the .1 file and make sure that the
      !   directions are the same for each frequency component:

   K         = 0        ! Initialize to zero
   PrvPer    = 0.0      ! Initialize to a don't care
   PrvDir    = 0.0      ! Initialize to a don't care
   FirstPass = .TRUE.   ! Initialize to .TRUE. for the first pass

   DO    ! Loop through all rows in the file


      READ (UnW3,'(A)',IOSTAT=Sttus)  Line   ! Read in the entire line

      IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully


         READ (Line,*)  TmpPer, TmpDir ! Read in only the period and direction from the WAMIT file


         IF ( FirstPass .OR. ( TmpPer /= PrvPer ) )  THEN   ! .TRUE. if we are on the first pass or if the period    currently read in is different than the previous period    read in; thus we found a new period    in the WAMIT file!

            J         = 0           ! Reset the count of directions to zero
            K         = K + 1       ! This is current count of which frequency component we are on
            PrvPer    = TmpPer      ! Store the current period    as the previous period    for the next pass
            FirstFreq = FirstPass   ! Sorry, you can only loop through the first frequency once
            NewPer    = .TRUE.      ! Reset the new period flag

            DO WHILE ( WAMITPer(K) <= 0.0 )  ! Periods less than or equal to zero in WAMIT represent infinite period = zero frequency and infinite frequency, respectively.  However, only the added mass is output by WAMIT at these limits.  The damping and wave excitation are left blank, so skip them!
               K = K + 1
            ENDDO

            IF ( TmpPer /= WAMITPer(K) )  THEN  ! Abort if the .3 and .1 files do not contain the same frequency components (not counting zero and infinity)
               CALL ProgAbort ( ' Other than zero and infinite frequencies, "'   //TRIM(FltPtfm_InitData%WAMITFile)//'.3",' // &
                            ' contains different frequency components than "'//TRIM(FltPtfm_InitData%WAMITFile)//'.1". '// &
                            ' Both WAMIT output files must be generated from the same run.', TrapErrors = .TRUE.)
               ErrStat = 1
               RETURN
            ENDIF

         ENDIF


         IF ( FirstPass .OR. ( TmpDir /= PrvDir ) .OR. NewPer )  THEN   ! .TRUE. if we are on the first pass, or if this is new period, or if the direction currently read in is different than the previous direction read in; thus we found a new direction in the WAMIT file!

            J         = J + 1       ! This is current count of which direction component we are on
            PrvDir    = TmpDir      ! Store the current direction as the previous direction for the next pass
            FirstPass = .FALSE.     ! Sorry, you can only have one first pass
            NewPer    = .FALSE.     ! Disable the new period flag

            IF ( FirstFreq )  THEN                    ! .TRUE. while we are still looping through all directions for the first frequency component
               WAMITWvDir(J)   = TmpDir      ! Store the directions in the order they appear in the WAMIT file

               InsertInd       = J           ! Initialize as the J'th component
               DO I = 1,J-1   ! Loop throuh all previous directions
                  IF ( ( WAMITWvDir(I) > WAMITWvDir(J) ) )  THEN  ! .TRUE. if a previous direction component is higher than the current direction component
                     InsertInd       = MIN( InsertInd, SortWvDirInd(I) )   ! Store the lowest sorted index whose associated direction component is higher than the current direction component
                     SortWvDirInd(I) = SortWvDirInd(I) + 1                 ! Shift all of the sorted indices up by 1 whose associated direction component is higher than the current direction component
                  ENDIF
               ENDDO          ! I - All previous directions
               SortWvDirInd(J) = InsertInd   ! Store the index such that WAMITWvDir(SortWvDirInd(:)) is sorted from lowest to highest direction
            ELSEIF ( TmpDir /= WAMITWvDir(J) )  THEN  ! We must have looped through all directions at least once; so check to make sure all subsequent directions are consistent with the directions from the first frequency component, otherwise Abort
               CALL ProgAbort ( ' Not every frequency component in "'//TRIM(FltPtfm_InitData%WAMITFile)//'.3"'// &
                            ' contains the same listing of direction angles.  Check for' // &
                            ' errors in the WAMIT output file.', TrapErrors = .TRUE.)
               ErrStat = 1
               RETURN
            ENDIF

         ENDIF


      ELSE                    ! We must have reached the end of the file, so stop reading in data


         EXIT


      ENDIF


   ENDDO ! End loop through all rows in the file


   REWIND (UNIT=UnW3)   ! REWIND the file so we can read it in a third time.  (This is getting ridiculous!)


      ! Now we can finally read in the frequency- and direction-dependent complex
      !   wave excitation force per unit wave amplitude vector:

   K                = 0       ! Initialize to zero
   PrvPer           = 0.0     ! Initialize to a don't care
   PrvDir           = 0.0     ! Initialize to a don't care
   FirstPass        = .TRUE.  ! Initialize to .TRUE. for the first pass

   HdroExctn(:,:,:) = 0.0     ! Initialize to zero

   DO    ! Loop through all rows in the file


      READ (UnW3,'(A)',IOSTAT=Sttus)  Line   ! Read in the entire line

      IF ( Sttus == 0 )  THEN ! .TRUE. when data is read in successfully


         READ (Line,*)  TmpPer, TmpDir, I, TmpData1, TmpData2, TmpRe, TmpIm   ! Read in the period, direction, row index, and nondimensional data from the WAMIT file


         IF ( FirstPass .OR. ( TmpPer /= PrvPer ) )  THEN   ! .TRUE. if we are on the first pass or if the period    currently read in is different than the previous period    read in; thus we found a new period    in the WAMIT file!

            J            = 0           ! Reset the count of directions to zero
            K            = K + 1       ! This is current count of which frequency component we are on
            PrvPer       = TmpPer      ! Store the current period    as the previous period    for the next pass
            FirstFreq    = FirstPass   ! Sorry, you can only loop through the first frequency once
            NewPer       = .TRUE.      ! Reset the new period flag

            DO WHILE ( WAMITPer(K) <= 0.0 )  ! Periods less than or equal to zero in WAMIT represent infinite period = zero frequency and infinite frequency, respectively.  However, only the added mass is output by WAMIT at these limits.  The damping and wave excitation are left blank, so skip them!
               K = K + 1
            ENDDO

         ENDIF


         IF ( FirstPass .OR. ( TmpDir /= PrvDir ) .OR. NewPer )  THEN   ! .TRUE. if we are on the first pass, or if this is new period, or if the direction currently read in is different than the previous direction read in; thus we found a new direction in the WAMIT file!

            J            = J + 1       ! This is current count of which direction component we are on
            PrvDir       = TmpDir      ! Store the current direction as the previous direction for the next pass
            FirstPass    = .FALSE.     ! Sorry, you can only have one first pass
            NewPer       = .FALSE.     ! Disable the new period flag

            IF ( FirstFreq )  THEN  ! .TRUE. while we are still looping through all directions for the first frequency component
               HdroWvDir(SortWvDirInd(J)) = TmpDir ! Store the directions sorted from lowest to highest
            ENDIF

         ENDIF


         HdroExctn(SortFreqInd(K),SortWvDirInd(J),I) = CMPLX( TmpRe, TmpIm )*DffrctDim(I) ! Redimensionalize the data and place it at the appropriate location within the array


      ELSE                    ! We must have reached the end of the file, so stop reading in data


         EXIT


      ENDIF


   ENDDO ! End loop through all rows in the file


   CLOSE ( UnW3 ) ! Close file.


   ! For some reason, WAMIT computes the zero- and infinite- frequency limits for
   !   only the added mass.  Based on hydrodynamic theory, the damping is zero at
   !   these limits (as initialized).  Hydrodynamic theory also says that the
   !   infinite-frequency limit of the diffraction force is zero (as initialized);
   !   however, the zero-frequency limit need not be zero.  Thus, if necessary
   !   (i.e., if we have read in a WAMIT output file that contains the
   !   zero-frequency limit of the added mass), compute the zero-frequency limit
   !   of the diffraction problem using the known values at the lowest
   !   nonzero-valued frequency available:

   DO I = 1,NInpFreq       ! Loop through all input frequency components

      IF ( HdroFreq(I) > 0.0 )  THEN ! .TRUE. at the lowest nonzero-valued frequency component

         DO J = I-1,1,-1   ! Loop through all zero-valued frequency components
            HdroExctn(J,:,:) = HdroExctn(I,:,:) ! Set the zero-frequency limits to equal the known values at the lowest nonzero-valued frequency available
         ENDDO             ! J - All zero-valued frequency components

         EXIT  ! Since HdroFreq(:) is sorted from lowest to highest frequency, there is no reason to continue on once we have found the lowest nonzero-valued frequency component

      ENDIF

   ENDDO                   ! I - All input frequency components




      ! Tell our nice users what is about to happen that may take a while:

   CALL WrScr ( ' Computing radiation impulse response functions and wave diffraction forces.' )



      ! Abort if the WAMIT files do not contain both the zero- and and infinite-
      !   frequency limits of added mass.  Also, if HighFreq is greater than
      !   RdtnOmegaMax, Abort because RdtnDT must be reduced in order to have
      !   sufficient accuracy in the computation of the radiation impulse response
      !   functions:

   IF ( .NOT. ( ZeroFreq .AND. InfFreq ) )  THEN   ! .TRUE. if both the zero- and infinite-frequency limits of added mass are contained within the WAMIT file
      CALL ProgAbort ( ' "'//TRIM(FltPtfm_InitData%WAMITFile)// &
                       '.1" must contain both the zero- and infinite-frequency limits of added mass.', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( HighFreq > RdtnOmegaMax      )  THEN   ! .TRUE. if the highest frequency component (not counting infinity) in the WAMIT file is greater than RdtnOmegaMax
      CALL ProgAbort ( ' Based on the frequency range found in "'//TRIM(FltPtfm_InitData%WAMITFile)//'.1",'       // &
                   ' RdtnDT must be set smaller than '//TRIM(Flt2LStr( Pi/HighFreq ))//' sec'// &
                   ' in order to accurately compute the radiation impulse response functions.', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ENDIF



      ! Set the infinite-frequency limit of the frequency-dependent hydrodynamic
      !   added mass matrix, HdroAdMsI, based on the highest frequency available:

   Indx = 0
   DO J = 1,6        ! Loop through all rows    of HdroAdMsI
      DO K = J,6     ! Loop through all columns of HdroAdMsI above and including the diagonal
         Indx = Indx + 1
         FP_Data%HdroAdMsI(J,K) = HdroAddMs(NInpFreq,Indx)
      ENDDO          ! K - All columns of HdroAdMsI above and including the diagonal
      DO K = J+1,6   ! Loop through all rows    of HdroAdMsI below the diagonal
         FP_Data%HdroAdMsI(K,J) = FP_Data%HdroAdMsI(J,K)
      ENDDO          ! K - All rows    of HdroAdMsI below the diagonal
   ENDDO             ! J - All rows    of HdroAdMsI



      ! Perform some initialization computations including calculating the total
      !   number of frequency components = total number of time steps in the wave,
      !   radiation kernel, calculating the frequency step, and ALLOCATing the
      !   arrays:
      ! NOTE: RdtnDOmega = Pi/RdtnTMax since, in the (co)sine transforms:
      !          Omega = (K-1)*RdtnDOmega
      !          Time  = (J-1)*RdtnDT
      !       and therefore:
      !          Omega*Time = (K-1)*(J-1)*RdtnDOmega*RdtnDT
      !                     = (K-1)*(J-1)*Pi/(NStepRdtn-1) [see FFT_Module]
      !       or:
      !          RdtnDOmega = Pi/((NStepRdtn-1)*RdtnDT)
      !                     = Pi/RdtnTMax

   FP_Data%NStepRdtn  = CEILING ( FP_Data%RdtnTMax/FP_Data%RdtnDT )                 ! Set NStepRdtn to an odd integer
   IF ( MOD(FP_Data%NStepRdtn,2) == 0 )  FP_Data%NStepRdtn = FP_Data%NStepRdtn + 1  !   larger or equal to RdtnTMax/RdtnDT.
   NStepRdtn2 = MAX( ( FP_Data%NStepRdtn-1 )/2, 1 )                                 ! Make sure that NStepRdtn-1 is an even product of small factors (PSF) that is greater
   FP_Data%NStepRdtn  = 2*PSF ( NStepRdtn2, 9 ) + 1                                 !   or equal to RdtnTMax/RdtnDT to ensure that the (co)sine transform is efficient.

   FP_Data%NStepRdtn1 = FP_Data%NStepRdtn + 1                                       ! Save the value of NStepRdtn + 1 for future use.
   NStepRdtn2 = ( FP_Data%NStepRdtn-1 )/2                                           ! Update the value of NStepRdtn2 based on the value needed for NStepRdtn.
   FP_Data%RdtnTMax   = ( FP_Data%NStepRdtn-1 )*FP_Data%RdtnDT                      ! Update the value of RdtnTMax   based on the value needed for NStepRdtn.
   RdtnDOmega = Pi/FP_Data%RdtnTMax                                                 ! Compute the frequency step for wave radiation kernel calculations.

   ALLOCATE ( RdtnTime (0:FP_Data%NStepRdtn-1    ) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the RdtnTime array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

   ALLOCATE ( FP_Data%RdtnKrnl (0:FP_Data%NStepRdtn-1,6,6) , STAT=ErrStat )
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the RdtnKrnl array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF

   ALLOCATE ( FP_Data%XDHistory(0:FP_Data%NStepRdtn  ,6  ) , STAT=ErrStat )   ! In the numerical convolution we must have NStepRdtn1 elements within the XDHistory array, which is one more than the NStepRdtn elements that are in the RdtnKrnl array
   IF ( ErrStat /= 0 )  THEN
      CALL ProgAbort(' Error allocating memory for the XDHistory array.', TrapErrors = .TRUE.)
      RETURN
   ENDIF



   IF ( RdtnFrmAM )  THEN  ! .TRUE. if we will determine the wave radiation kernel from the frequency-dependent hydrodynamic added mass matrix



      ! Calculate the factor needed by the discrete sine transform in the
      !   calculation of the wave radiation kernel:

      Krnl_Fact = -1.0/FP_Data%RdtnDT ! This factor is needed by the discrete time sine transform



      ! Compute all frequency components (including zero) of the sine transform
      !   of the wave radiation kernel:

      DO I = 0,FP_Data%NStepRdtn-1 ! Loop through all frequency components (including zero) of the sine transform


      ! Calculate the array of simulation times at which the instantaneous values
      !   of the wave radiation kernel are to be determined:

         RdtnTime(I) = I*FP_Data%RdtnDT


      ! Compute the frequency of this component:

         Omega = I*RdtnDOmega


      ! Compute the upper-triangular portion (diagonal and above) of the sine
      !   transform of the wave radiation kernel:

         Indx = 0
         DO J = 1,6        ! Loop through all rows    of RdtnKrnl
            DO K = J,6     ! Loop through all columns of RdtnKrnl above and including the diagonal
               Indx = Indx + 1
               FP_Data%RdtnKrnl(I,J,K) = Krnl_Fact*Omega*( InterpStp( Omega, HdroFreq(:), &
                                                                             HdroAddMs(:       ,Indx), LastInd, NInpFreq ) &
                                                   -                         HdroAddMs(NInpFreq,Indx)                      )
            ENDDO          ! K - All columns of RdtnKrnl above and including the diagonal
         ENDDO             ! J - All rows    of RdtnKrnl


      ENDDO                ! I - All frequency components (including zero) of the sine transform



      ! Compute the sine transforms to find the time-domain representation of
      !   the wave radiation kernel:

!bjj start of proposed change v1.00.00a-bjj
!i added the FFT_Data and ErrStat parameters to the FFT_Module subroutine calls
!rm      CALL InitSINT ( NStepRdtn, .TRUE. )
      CALL InitSINT ( FP_Data%NStepRdtn, FFT_Data, .TRUE., ErrStat )
      IF ( ErrStat /= 0 ) RETURN

      DO J = 1,6                 ! Loop through all rows    of RdtnKrnl
         DO K = J,6              ! Loop through all columns of RdtnKrnl above and including the diagonal
            CALL ApplySINT( FP_Data%RdtnKrnl(:,J,K), FFT_Data, ErrStat )
            IF ( ErrStat /= 0 ) RETURN
         ENDDO                   ! K - All columns of RdtnKrnl above and including the diagonal
         DO K = J+1,6            ! Loop through all rows    of RdtnKrnl below the diagonal
            DO I = 0,FP_Data%NStepRdtn-1 ! Loop through all frequency components (including zero) of the sine transform
               FP_Data%RdtnKrnl(I,K,J) = FP_Data%RdtnKrnl(I,J,K)
            ENDDO                ! I - All frequency components (including zero) of the sine transform
         ENDDO                   ! K - All rows    of RdtnKrnl below the diagonal
      ENDDO                      ! J - All rows    of RdtnKrnl

      CALL ExitSINT(FFT_Data, ErrStat)
      IF ( ErrStat /= 0 ) RETURN
!bjj end of proposed change


   ELSE                    ! We must be determining the wave radiation kernel from the frequency-dependent hydrodynamic damping matrix



      ! Calculate the factor needed by the discrete cosine transform in the
      !   calculation of the wave radiation kernel:

      Krnl_Fact = 1.0/FP_Data%RdtnDT  ! This factor is needed by the discrete time cosine transform



      ! Compute all frequency components (including zero) of the cosine transform
      !   of the wave radiation kernel:

      DO I = 0,FP_Data%NStepRdtn-1 ! Loop through all frequency components (including zero) of the cosine transform


      ! Calculate the array of simulation times at which the instantaneous values
      !   of the wave radiation kernel are to be determined:

         RdtnTime(I) = I*FP_Data%RdtnDT


      ! Compute the frequency of this component:

         Omega = I*RdtnDOmega


      ! Compute the upper-triangular portion (diagonal and above) of the cosine
      !   transform of the wave radiation kernel:

         Indx = 0
         DO J = 1,6        ! Loop through all rows    of RdtnKrnl
            DO K = J,6     ! Loop through all columns of RdtnKrnl above and including the diagonal
               Indx = Indx + 1
               FP_Data%RdtnKrnl(I,J,K) = Krnl_Fact*InterpStp ( Omega, HdroFreq(:), HdroDmpng(:,Indx), LastInd, NInpFreq )
            ENDDO          ! K - All columns of RdtnKrnl above and including the diagonal
         ENDDO             ! J - All rows    of RdtnKrnl


      ENDDO                ! I - All frequency components (including zero) of the cosine transform



      ! Compute the cosine transforms to find the time-domain representation of
      !   the wave radiation kernel:

!bjj start of proposed change v1.00.00a-bjj
!i added the FFT_Data and ErrStat parameters to the FFT_Module subroutine calls
      CALL InitCOST ( FP_Data%NStepRdtn, FFT_Data, .TRUE., ErrStat )
      IF ( ErrStat /= 0 ) RETURN

      DO J = 1,6                          ! Loop through all rows    of RdtnKrnl
         DO K = J,6                       ! Loop through all columns of RdtnKrnl above and including the diagonal
            CALL ApplyCOST( FP_Data%RdtnKrnl(:,J,K), FFT_Data, ErrStat )
            IF ( ErrStat /= 0 ) RETURN
         ENDDO                            ! K - All columns of RdtnKrnl above and including the diagonal
         DO K = J+1,6                     ! Loop through all rows    of RdtnKrnl below the diagonal
            DO I = 0,FP_Data%NStepRdtn-1  ! Loop through all radiation time steps
               FP_Data%RdtnKrnl(I,K,J) = FP_Data%RdtnKrnl(I,J,K)
            ENDDO                         ! I - All radiation time steps
         ENDDO                            ! K - All rows    of RdtnKrnl below the diagonal
      ENDDO                               ! J - All rows    of RdtnKrnl

      CALL ExitCOST(FFT_Data, ErrStat)
      IF ( ErrStat /= 0 ) RETURN
!bjj end of proposed change


   ENDIF
!JASON:USE THIS TO TEST ADDED MASS:DO I = 1,6  !JASON:USE THIS TO TEST ADDED MASS:
!JASON:USE THIS TO TEST ADDED MASS:DO J = 1,6  !JASON:USE THIS TO TEST ADDED MASS:
!JASON:USE THIS TO TEST ADDED MASS:   WRITE (*,*) I, J, HdroAdMsI(I,J) !JASON:USE THIS TO TEST ADDED MASS:
!JASON:USE THIS TO TEST ADDED MASS:ENDDO  !JASON:USE THIS TO TEST ADDED MASS:
!JASON:USE THIS TO TEST ADDED MASS:ENDDO  !JASON:USE THIS TO TEST ADDED MASS:
!JASON:USE THIS TO TEST IMPULSE RESPONSE FUNCTIONS:DO I = 0,NStepRdtn-1 ! Loop through all radiation time steps   !JASON:USE THIS TO TEST IMPULSE RESPONSE FUNCTIONS:
!JASON:USE THIS TO TEST IMPULSE RESPONSE FUNCTIONS:   WRITE (*,*) RdtnTime(I), RdtnKrnl(I,1,1), RdtnKrnl(I,3,3)   !JASON:USE THIS TO TEST IMPULSE RESPONSE FUNCTIONS:
!JASON:USE THIS TO TEST IMPULSE RESPONSE FUNCTIONS:   WRITE (*,*) RdtnTime(I), RdtnKrnl(I,4,4), RdtnKrnl(I,6,6)   !JASON:USE THIS TO TEST IMPULSE RESPONSE FUNCTIONS:
!JASON:USE THIS TO TEST IMPULSE RESPONSE FUNCTIONS:   WRITE (*,*) RdtnTime(I), RdtnKrnl(I,1,5), RdtnKrnl(I,2,4)   !JASON:USE THIS TO TEST IMPULSE RESPONSE FUNCTIONS:
!JASON:USE THIS TO TEST IMPULSE RESPONSE FUNCTIONS:ENDDO                ! I - All radiation time steps   !JASON:USE THIS TO TEST IMPULSE RESPONSE FUNCTIONS:




      ! Initialize the variables associated with the incident wave:

   SELECT CASE ( WaveDat%WaveMod ) ! Which incident wave kinematics model are we using?

   CASE ( 0 )              ! None=still water.



      ! Initialize everything to zero:

      ALLOCATE ( FP_Data%WaveExctn (0:WaveDat%NStepWave-1,6) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveExctn array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      FP_Data%WaveExctn = 0.0




!jbj: start of proposed change v1.00.00b-jbj
!rm   CASE ( 1, 2, 3 )        ! Plane progressive (regular) wave, JONSWAP/Pierson-Moskowitz spectrum (irregular) wave, or user-defined spectrum (irregular) wave.
   CASE ( 1, 2, 3, 10 )    ! Plane progressive (regular) wave, JONSWAP/Pierson-Moskowitz spectrum (irregular) wave, or user-defined spectrum (irregular) wave.
!jbj: end of proposed change v1.00.00b-jbj



      ! Abort if we have chosen a wave heading direction that is outside the range
      !   of directions where the complex wave excitation force per unit wave
      !   amplitude vector has been defined, else interpolate to find the complex
      !   wave excitation force per unit wave amplitude vector at the chosen wave
      !   heading direction:

      IF ( ( WaveDat%WaveDir < HdroWvDir(1) ) .OR. ( WaveDat%WaveDir > HdroWvDir(NInpWvDir) ) )  THEN
         CALL ProgAbort ( ' WaveDir must be within the wave heading angle range available in "' &
                           //TRIM(FltPtfm_InitData%WAMITFile)//'.3" (inclusive).', TrapErrors = .TRUE.)
         ErrStat = 1
         RETURN
      END IF

      ALLOCATE ( X_Diffrctn(NInpFreq,6) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the X_Diffrctn array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      DO J = 1,6           ! Loop through all wave excitation forces and moments
         DO I = 1,NInpFreq ! Loop through all input frequency components inherent in the complex wave excitation force per unit wave amplitude vector
            X_Diffrctn(I,J) = InterpStp( WaveDat%WaveDir, HdroWvDir(:), HdroExctn(I,:,J), LastInd, NInpWvDir )
         ENDDO             ! I - All input frequency components inherent in the complex wave excitation force per unit wave amplitude vector
      ENDDO                ! J - All wave excitation forces and moments



      ! ALLOCATE the arrays:

      ALLOCATE (         WaveExctnC(0:WaveDat%NStepWave2 ,6) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveExctnC array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF

      ALLOCATE ( FP_Data%WaveExctn (0:WaveDat%NStepWave-1,6) , STAT=ErrStat )
      IF ( ErrStat /= 0 )  THEN
         CALL ProgAbort(' Error allocating memory for the WaveExctn array.', TrapErrors = .TRUE.)
         RETURN
      ENDIF


!jbj: start of proposed change v1.00.00b-jbj
!rm
!rm      ! Compute the positive-frequency components (including zero) of the Fourier
!rm      !  transforms of the wave excitation force:
!rm
!rm      DO I = 0,WaveDat%NStepWave2  ! Loop through the positive frequency components (including zero) of the Fourier transforms
      ! Compute the positive-frequency components (including zero) of the discrete
      !   Fourier transform of the wave excitation force:

      DO I = 0,WaveDat%NStepWave2  ! Loop through the positive frequency components (including zero) of the discrete Fourier transform
!jbj: end of proposed change v1.00.00b-jbj


      ! Compute the frequency of this component:

         Omega = I*WaveDat%WaveDOmega

!jbj: start of proposed change v1.00.00b-jbj
!rm      ! Compute the Fourier transform of the instantaneous value of the total
!rm      !   excitation force on the support platfrom from incident waves:
      ! Compute the discrete Fourier transform of the instantaneous value of the
      !   total excitation force on the support platfrom from incident waves:
!jbj: end of proposed change v1.00.00b-jbj

         DO J = 1,6           ! Loop through all wave excitation forces and moments
            WaveExctnC(I,J) = WaveDat%WaveElevC0(I)*InterpStp ( Omega, HdroFreq(:), X_Diffrctn(:,J), LastInd, NInpFreq )
         ENDDO                ! J - All wave excitation forces and moments


!jbj: start of proposed change v1.00.00b-jbj
!rm      ENDDO                ! I - The positive frequency components (including zero) of the Fourier transforms
!rm
!rm
!rm      ! Compute the inverse Fourier transforms to find the time-domain
!rm      !   representations of the wave excitation forces:
      ENDDO                ! I - The positive frequency components (including zero) of the discrete Fourier transform



      ! Compute the inverse discrete Fourier transform to find the time-domain
      !   representation of the wave excitation force:
!jmj End of proposed change.  v6.10d-jmj  13-Aug-2009.

!jbj: end of proposed change v1.00.00b-jbj

!bjj start of proposed change v1.00.00a-bjj
!bjj i added FFT_Data and ErrStat arguments to the subroutine calls in the FFT_Module
      CALL InitFFT ( WaveDat%NStepWave, FFT_Data, .TRUE., ErrStat )
      IF ( ErrStat /= 0 ) RETURN

      DO J = 1,6           ! Loop through all wave excitation forces and moments
         CALL ApplyFFT_cx ( FP_Data%WaveExctn(:,J), WaveExctnC(:,J), FFT_Data, ErrStat )
         IF ( ErrStat /= 0 ) RETURN
      ENDDO                ! J - All wave excitation forces and moments

      CALL ExitFFT(FFT_Data, ErrStat)
      IF ( ErrStat /= 0 ) RETURN
!bjj end of proposed change



   CASE ( 4 )              ! GH Bladed wave data.



      CALL ProgAbort ( ' GH Bladed wave data not applicable for floating platforms. ', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN




   ENDSELECT


      ! deallocate arrays

   IF ( ALLOCATED( HdroExctn    ) ) DEALLOCATE( HdroExctn    )
   IF ( ALLOCATED( WaveExctnC   ) ) DEALLOCATE( WaveExctnC   )
   IF ( ALLOCATED( X_Diffrctn   ) ) DEALLOCATE( X_Diffrctn   )
   IF ( ALLOCATED( HdroAddMs    ) ) DEALLOCATE( HdroAddMs    )
   IF ( ALLOCATED( HdroDmpng    ) ) DEALLOCATE( HdroDmpng    )
   IF ( ALLOCATED( HdroFreq     ) ) DEALLOCATE( HdroFreq     )
   IF ( ALLOCATED( HdroWvDir    ) ) DEALLOCATE( HdroWvDir    )
   IF ( ALLOCATED( RdtnTime     ) ) DEALLOCATE( RdtnTime     )
   IF ( ALLOCATED( WAMITFreq    ) ) DEALLOCATE( WAMITFreq    )
   IF ( ALLOCATED( WAMITPer     ) ) DEALLOCATE( WAMITPer     )
   IF ( ALLOCATED( WAMITWvDir   ) ) DEALLOCATE( WAMITWvDir   )
   IF ( ALLOCATED( SortFreqInd  ) ) DEALLOCATE( SortFreqInd  )
   IF ( ALLOCATED( SortWvDirInd ) ) DEALLOCATE( SortWvDirInd )

!bjj start of proposed change v1.00.00a-bjj   
   IF ( ALLOCATED( User_LineNodes ) ) DEALLOCATE( User_LineNodes )
!bjj end of proposed change   

   RETURN
   END SUBROUTINE InitFltngPtfmLd
!=======================================================================
   FUNCTION LinePosition ( ILine, JNode, KDirection, FP_Data, ErrStat )


      ! This FUNCTION is used to return the instantaneous line position at
      ! node JNode of mooring line ILine in the xi- (KDirection=1), yi-
      ! (KDirection=2), or zi- (KDirection=3) direction, respectively, to
      ! the calling program.


   IMPLICIT                        NONE


      ! Passed Variables:

   REAL(ReKi)                        :: LinePosition                             ! This function = instantaneous line position at node JNode of mooring line ILine in the inertia frame (meters)

   INTEGER,    INTENT(IN )           :: ILine                                    ! Mooring line number (-)
   INTEGER,    INTENT(IN )           :: JNode                                    ! The index of the current mooring line node (-)
   INTEGER,    INTENT(IN )           :: KDirection                               ! 1, 2, or 3, for the xi-, yi-, or zi-directions, respectively (-)

!bjj start of proposed change v1.00.00a-bjj
   TYPE(FltPtfm_DataType),INTENT(IN) :: FP_Data
   INTEGER,    INTENT(OUT)           :: ErrStat                                  ! a non-zero value indicates an error has occurred

   ErrStat = 0
!bjj end of proposed change


      ! Abort if the mooring line parameters have not been computed yet, if ILine
      !   is not one of the existing mooring lines, if JNode is not one of the
      !   existing line nodes, or if KDirection is not specified properly:

!bjj start of propsoed change v1.00.00a-bjj
!bjj: this must be reorganized due to new data types
!   IF ( .NOT. ALLOCATED ( LNodesPi )                   )  THEN
!      CALL ProgAbort ( ' Routine InitFltngPtfmLd() must be called before routine LinePosition().', TrapErrors = .TRUE.)
!      ErrStat = 1
!      RETURN
!   ELSEIF ( ( ILine < 1 ) .OR. ( ILine > FP_Data%NumLines  )   )  THEN
!      CALL ProgAbort ( ' Mooring line '//TRIM( Int2LStr( ILine ) )//' has not been analyzed.', TrapErrors = .TRUE.)
!      ErrStat = 1
!      RETURN
!   ELSEIF ( ( JNode < 1 ) .OR. ( JNode > LineNodes )   )  THEN
!      CALL ProgAbort ( ' Line node '   //TRIM( Int2LStr( Jnode ) )//' has not been analyzed.', TrapErrors = .TRUE.)
!      ErrStat = 1
!      RETURN
!   ELSEIF ( ( KDirection < 1 ) .OR. ( KDirection > 3 ) )  THEN
!      CALL ProgAbort ( ' KDirection must be 1, 2, or 3 in routine LinePosition().'          , TrapErrors = .TRUE.)
!      ErrStat = 1
!      RETURN
!   ENDIF

   IF ( ( ILine < 1 ) .OR. ( ILine > FP_Data%NumLines  )                       )  THEN
      CALL ProgAbort ( ' Mooring line '//TRIM( Int2LStr( ILine ) )//' has not been analyzed.'    , TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( .NOT. ALLOCATED ( FP_Data%MooringLine )                             )  THEN
      CALL ProgAbort ( ' Routine InitFltngPtfmLd() must be called before routine LinePosition().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( ( JNode < 1 ) .OR. ( JNode > FP_Data%MooringLine(ILine)%LineNodes ) )  THEN
      CALL ProgAbort ( ' Line node '   //TRIM( Int2LStr( Jnode ) )//' has not been analyzed.'    , TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( .NOT. ALLOCATED ( FP_Data%MooringLine(ILine)%LNodesPi )             )  THEN
      CALL ProgAbort ( ' Routine InitFltngPtfmLd() must be called before routine LinePosition().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( ( KDirection < 1 ) .OR. ( KDirection > 3 ) )  THEN
      CALL ProgAbort ( ' KDirection must be 1, 2, or 3 in routine LinePosition().'               , TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ENDIF
!bjj end of proposed change

      ! Return the instantaneous line position:

   LinePosition = FP_Data%MooringLine(ILine)%LNodesPi(JNode,KDirection)



   RETURN
   END FUNCTION LinePosition
!=======================================================================
   FUNCTION LineTension ( ILine, JNode, FP_Data, ErrStat )


      ! This FUNCTION is used to return the instantaneous effective line
      ! tension at node JNode of mooring line ILine to the calling program.



   IMPLICIT                        NONE


      ! Passed Variables:

   REAL(ReKi)                        :: LineTension                              ! This function = instantaneous effective line tension at node JNode of mooring line ILine (N)

   INTEGER,    INTENT(IN )           :: ILine                                    ! Mooring line number (-)
   INTEGER,    INTENT(IN )           :: JNode                                    ! The index of the current mooring line node (-)

!bjj start of proposed change v1.00.00a-bjj
   TYPE(FltPtfm_DataType),INTENT(IN) :: FP_Data                                  ! data for this instance of the floating platform module
   INTEGER,    INTENT(OUT)           :: ErrStat                                  ! a non-zero value indicates an error has occurred

   ErrStat = 0
!bjj end of proposed change

      ! Abort if the mooring line parameters have not been computed yet, if ILine
      !   is not one of the existing mooring lines, or if JNode is not one of the
      !   existing line nodes:

!bjj start of proposed change v1.00.00a-bjj
!bjj this needs to be reordered due to new data types
!   IF ( .NOT. ALLOCATED ( LNodesTe )                 )  THEN
!      CALL ProgAbort ( ' Routine InitFltngPtfmLd() must be called before routine LineTension().', TrapErrors = .TRUE.)
!      ErrStat = 1
!      RETURN
!   ELSEIF ( ( ILine < 1 ) .OR. ( ILine > FP_Data%NumLines  ) )  THEN
!      CALL ProgAbort ( ' Mooring line '//TRIM( Int2LStr( ILine ) )//' has not been analyzed.', TrapErrors = .TRUE.)
!      ErrStat = 1
!      RETURN
!   ELSEIF ( ( JNode < 1 ) .OR. ( JNode > LineNodes ) )  THEN
!      CALL ProgAbort ( ' Line node '   //TRIM( Int2LStr( Jnode ) )//' has not been analyzed.', TrapErrors = .TRUE.)
!      ErrStat = 1
!      RETURN
!   ENDIF

   IF ( ( ILine < 1 ) .OR. ( ILine > FP_Data%NumLines  )                        )  THEN
      CALL ProgAbort ( ' Mooring line '//TRIM( Int2LStr( ILine ) )//' has not been analyzed.'   , TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( .NOT. ALLOCATED ( FP_Data%MooringLine )                             )  THEN
      CALL ProgAbort ( ' Routine InitFltngPtfmLd() must be called before routine LineTension().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( ( JNode < 1 ) .OR. ( JNode > FP_Data%MooringLine(ILine)%LineNodes ) )  THEN
      CALL ProgAbort ( ' Line node '   //TRIM( Int2LStr( JNode ) )//' has not been analyzed.'   , TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ELSEIF ( .NOT. ALLOCATED ( FP_Data%MooringLine(ILine)%LNodesTe )             )  THEN
      CALL ProgAbort ( ' Routine InitFltngPtfmLd() must be called before routine LineTension().', TrapErrors = .TRUE.)
      ErrStat = 1
      RETURN
   ENDIF
!bjj end of proposed change v1.00.00a-bjj



      ! Return the instantaneous effective tension:

   LineTension = FP_Data%MooringLine(ILine)%LNodesTe(JNode)



   RETURN
   END FUNCTION LineTension
!bjj start of proposed change v1.00.00a-bjj
!======================================================================= 
FUNCTION FP_GetNumLines ( FP_Data )
   
         ! Passed variables

   TYPE(FltPtfm_DataType),INTENT(IN)    :: FP_Data                                         ! data for this instance of the floating platform module
   REAL(ReKi)                           :: FP_GetNumLines
   
   
   FP_GetNumLines = FP_Data%NumLines
   
   
END FUNCTION FP_GetNumLines  
!======================================================================= 
FUNCTION FP_GetNumLineNodes ( FP_Data, ILine, ErrStat )
   
         ! Passed variables

   TYPE(FltPtfm_DataType),INTENT(IN)    :: FP_Data                                         ! data for this instance of the floating platform module
   INTEGER,               INTENT(IN)    :: ILine                                           ! the line for which to get number of nodes
   INTEGER,               INTENT(OUT)   :: ErrStat
   REAL(ReKi)                           :: FP_GetNumLineNodes
   
   ErrStat            = 0
   FP_GetNumLineNodes = 0
   
   IF ( FP_Data%NumLines < ILine ) THEN
      ErrStat = 1
      RETURN
   ELSE IF ( ALLOCATED( FP_Data%MooringLine ) ) THEN
      FP_GetNumLineNodes = FP_Data%MooringLine(ILine)%LineNodes
   END IF   
   
END FUNCTION FP_GetNumLineNodes  
!bjj end of proposed change v1.00.00a-bjj
!=======================================================================
   SUBROUTINE FP_Terminate(FP_Data, ErrStat )

         ! Passed variables

   TYPE(FltPtfm_DataType),INTENT(INOUT) :: FP_Data                                         ! data for this instance of the floating platform module
   INTEGER,               INTENT(OUT)   :: ErrStat

      ! Internal variables
   INTEGER                              :: I
   LOGICAL                              :: Err


      ! Initialize error status code

   ErrStat = 0
   Err     = .FALSE.


      ! Deallocate arrays

   IF ( ALLOCATED( FP_Data%MooringLine ) ) THEN

      DO I = 1,SIZE( FP_Data%MooringLine )

         IF ( ALLOCATED( FP_Data%MooringLine(I)%LNodesPi ) ) DEALLOCATE( FP_Data%MooringLine(I)%LNodesPi, STAT = ErrStat )
         IF ( ErrStat /= 0 ) Err = .TRUE.
         IF ( ALLOCATED( FP_Data%MooringLine(I)%LNodesTe ) ) DEALLOCATE( FP_Data%MooringLine(I)%LNodesTe, STAT = ErrStat )
         IF ( ErrStat /= 0 ) Err = .TRUE.
         IF ( ALLOCATED( FP_Data%MooringLine(I)%LNodesX  ) ) DEALLOCATE( FP_Data%MooringLine(I)%LNodesX , STAT = ErrStat )
         IF ( ErrStat /= 0 ) Err = .TRUE.
         IF ( ALLOCATED( FP_Data%MooringLine(I)%LNodesZ  ) ) DEALLOCATE( FP_Data%MooringLine(I)%LNodesZ , STAT = ErrStat )
         IF ( ErrStat /= 0 ) Err = .TRUE.
         IF ( ALLOCATED( FP_Data%MooringLine(I)%LSNodes  ) ) DEALLOCATE( FP_Data%MooringLine(I)%LSNodes , STAT = ErrStat )
         IF ( ErrStat /= 0 ) Err = .TRUE.

      END DO ! I

      DEALLOCATE ( FP_Data%MooringLine, STAT = ErrStat )
      IF ( ErrStat /= 0 ) Err = .TRUE.

   END IF

   FP_Data%NumLines = 0
   
!   IF ( ALLOCATED( LAnchHTe  ) ) DEALLOCATE(LAnchHTe , STAT = ErrStat )
!   IF ( ALLOCATED( LAnchVTe  ) ) DEALLOCATE(LAnchVTe , STAT = ErrStat )
!   IF ( ALLOCATED( LAnchxi   ) ) DEALLOCATE(LAnchxi  , STAT = ErrStat )
!   IF ( ALLOCATED( LAnchyi   ) ) DEALLOCATE(LAnchyi  , STAT = ErrStat )
!   IF ( ALLOCATED( LAnchzi   ) ) DEALLOCATE(LAnchzi  , STAT = ErrStat )
!   IF ( ALLOCATED( LEAStff   ) ) DEALLOCATE(LEAStff  , STAT = ErrStat )
!   IF ( ALLOCATED( LFairHTe  ) ) DEALLOCATE(LFairHTe , STAT = ErrStat )
!   IF ( ALLOCATED( LFairVTe  ) ) DEALLOCATE(LFairVTe , STAT = ErrStat )
!   IF ( ALLOCATED( LFairxt   ) ) DEALLOCATE(LFairxt  , STAT = ErrStat )
!   IF ( ALLOCATED( LFairyt   ) ) DEALLOCATE(LFairyt  , STAT = ErrStat )
!   IF ( ALLOCATED( LFairzt   ) ) DEALLOCATE(LFairzt  , STAT = ErrStat )
!   IF ( ALLOCATED( LFldWght  ) ) DEALLOCATE(LFldWght , STAT = ErrStat )
!   IF ( ALLOCATED( LNodesPi  ) ) DEALLOCATE(LNodesPi , STAT = ErrStat )
!   IF ( ALLOCATED( LNodesTe  ) ) DEALLOCATE(LNodesTe , STAT = ErrStat )
!   IF ( ALLOCATED( LNodesX   ) ) DEALLOCATE(LNodesX  , STAT = ErrStat )
!   IF ( ALLOCATED( LNodesZ   ) ) DEALLOCATE(LNodesZ  , STAT = ErrStat )
!   IF ( ALLOCATED( LSeabedCD ) ) DEALLOCATE(LSeabedCD, STAT = ErrStat )
!   IF ( ALLOCATED( LSNodes   ) ) DEALLOCATE(LSNodes  , STAT = ErrStat )
!   IF ( ALLOCATED( LTenTol   ) ) DEALLOCATE(LTenTol  , STAT = ErrStat )
!   IF ( ALLOCATED( LUnstrLen ) ) DEALLOCATE(LUnstrLen, STAT = ErrStat )

   IF ( ALLOCATED( FP_Data%RdtnKrnl  ) ) DEALLOCATE(FP_Data%RdtnKrnl , STAT = ErrStat)
   IF ( ErrStat /= 0 ) Err = .TRUE.
   IF ( ALLOCATED( FP_Data%WaveExctn ) ) DEALLOCATE(FP_Data%WaveExctn, STAT = ErrStat)
   IF ( ErrStat /= 0 ) Err = .TRUE.
   IF ( ALLOCATED( FP_Data%XDHistory ) ) DEALLOCATE(FP_Data%XDHistory, STAT = ErrStat)
   IF ( ErrStat /= 0 ) Err = .TRUE.


      ! Terminate UserLineModule
   CALL UserLine_Terminate( FP_Data%UserLine_Data, ErrStat )
   IF ( ErrStat /= 0 ) Err = .TRUE.
   

      ! Close files


      ! Reset variables
   FP_Data%CalculateFirstGuess    = .TRUE.      

      ! Return ErrStat = 1 if any error occurred while cleaning up

   IF ( Err ) ErrStat = 1



   END SUBROUTINE FP_Terminate

!=======================================================================

END MODULE FloatingPlatform
!=======================================================================
