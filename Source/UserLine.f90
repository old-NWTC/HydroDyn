MODULE UserLineModule
!
! The module that defines a user-specified mooring system.
! 
! The UserLine_Init() subroutine is called once at the beginning of the simulation to get the 
!     number of lines and line nodes.
! UserLine() is called to calculate the forces from the mooring lines at each time.
! UserLine_GetOutput() is called to get output values for the HydroDyn output file at each 
!     time and at each mooring line.
! UserLine_Terminate() is called once at the end of the simulation (and potentially called 
!     if an error is encountered)
!
! Please keep the following "rules" in mind:
! - The public types and subroutine interfaces defined in this module must not be changed;
!   however the contents inside the types and subroutines may be modified to achieve the user's
!   desired mooring line model.
! - NO global variables--including any variables with the "SAVE" attribute--
!   should be used in this module (so that multiple instances of HydroDyn can be called at the 
!   same time). Parameters are okay, as long as they are used only in this module. 
! - To store data that is accessed between subroutine calls, add the desired variable definitions 
!   to the UserLine_DataType. The calling program will keep track of this data (but will not be
!   able to read/write it).
! - Do not abort the code if an error occurs.  Always return a non-zero error code if there is an
!   error (set ErrStat = 0 if no error occurred) instead of actually aborting the code.
! - If any allocatable arrays have been declared in TYPE(UserLine_DataType) you must deallocate them
!   in UserLine_Terminate(); likewise, deallocate any local allocated arrays in subroutines/functions
! - Use caution in using a fixed unit number for file I/O. If multiple instances of the module are
!   in use, the instances could be reading from/writing to the same file. Consider using NWTC_Library's
!   GetNewUnit() routine instead.
!----------------------------------------------------------------------------------------------------

   USE                                             NWTC_Library
   
   IMPLICIT                                        NONE
   PRIVATE   
   
   
      !----------------------------------------------------------------------------------------------      
      ! Data types that are required to pass information to/from this module.     
      !----------------------------------------------------------------------------------------------      
   
   TYPE, PUBLIC :: UserLine_DataType
      PRIVATE                                                                       ! Data is private, though the data type definition is public
         ! Add your own variables here. (These variables take the place of any global or SAVEd variables.)
      CHARACTER(1024)   :: DirRoot                                                  ! a dummy variable to store DirRoot between calls. It is not used in this module, but is stored as an example.      
   END TYPE UserLine_DataType
   
   
      !----------------------------------------------------------------------------------------------      
      ! Public routines (all others are PRIVATE):
      !----------------------------------------------------------------------------------------------      
   PUBLIC                                       :: UserLine                         ! Do NOT change this line
   PUBLIC                                       :: UserLine_Init                    ! Do NOT change this line
   PUBLIC                                       :: UserLine_GetOutputs              ! Do NOT change this line
   PUBLIC                                       :: UserLine_Terminate               ! Do NOT change this line
         

CONTAINS
!----------------------------------------------------------------------------------------------------
   SUBROUTINE UserLine_Init( DirRoot, NumLines, LineNodes, UserLine_Data, ErrStat )
   ! This initialization subroutine is used to determine the number of mooring lines and number
   ! of nodes on each line. Users should change this so that it sets up the correct "NumLines" 
   ! (number of mooring lines) and "LineNodes" (analysis nodes) for each line. Note that "LineNodes"
   ! can vary by mooring line.
   !-------------------------------------------------------------------------------------------------
   
         ! Passed variables
         
      CHARACTER(1024),         INTENT( IN  )    :: DirRoot                       ! The name of the directory where the main input/output files are stored.
      INTEGER,                 INTENT( OUT )    :: ErrStat                       ! A non-zero number if an error occurs
      INTEGER, ALLOCATABLE,    INTENT( OUT )    :: LineNodes (:)                 ! Number of nodes on each line where the mooring line position and tension can be output (-)
      INTEGER,                 INTENT( OUT )    :: NumLines                      ! Number of mooring lines to be used

      TYPE(UserLine_DataType), INTENT( OUT )    :: UserLine_Data                 ! Internal data storage, which is initialized here
      
      
         ! Internal variables
         
      INTEGER                                   :: ILine                         ! Loop counter for mooring lines
      
      
         !-------------------------------------------------------------------------------------------         
         ! Initialize data
         !-------------------------------------------------------------------------------------------         
      ErrStat = 0                                                                ! No errors encountered, yet.   
      UserLine_Data%DirRoot = DirRoot                             
   
   
         !-------------------------------------------------------------------------------------------         
         ! Set the number of mooring lines
         !-------------------------------------------------------------------------------------------                  
      NumLines = 8                                                               ! USERS SHOULD CHANGE THIS NUMBER
         
         
         !-------------------------------------------------------------------------------------------         
         ! Allocate space for the LineNodes array
         !-------------------------------------------------------------------------------------------         
      ALLOCATE ( LineNodes( NumLines ), STAT=ErrStat )
      IF ( ErrStat /= 0 ) THEN
         CALL ProgAbort( ' Error allocating array for UserLine LineNodes.', TrapErrors = .TRUE. )
         RETURN
      END IF         
              
                     
         !-------------------------------------------------------------------------------------------         
         ! Set the number of analysis nodes on each line                         ! USERS SHOULD CHANGE THIS
         !-------------------------------------------------------------------------------------------                  
      DO ILine = 1,NumLines
      
         LineNodes(ILine) = 2                                                    ! USERS SHOULD CHANGE THIS NUMBER (note that this number can vary by mooring line, which may require taking this out of the DO LOOP)
                        
      END DO   ! ILine          
   
   
      RETURN
   
   END SUBROUTINE UserLine_Init
!----------------------------------------------------------------------------------------------------
   SUBROUTINE UserLine ( X, ZTime, F, UserLine_Data, ErrStat)
   !
   ! This is a dummy routine for holding the place of a user-specified mooring system. Modify
   ! this code to create your own model of an array of mooring lines. The local variables 
   ! and associated calculations below provide a template for making this user-specified
   ! mooring system model include a linear 6x6 restoring matrix with offset. These are 
   ! provided as an example only and can be modified or deleted as desired by the user 
   ! without detriment to the interface (i.e., they are not necessary for the interface).
   !
   ! The output of this routine is array F(:), which must contain the 3 components of
   ! the total force from all mooring lines (in N) acting at the platform reference and the 3
   ! components of the total moment from all mooring lines (in N-m) acting at the platform
   ! reference; positive forces are in the direction of positive platform displacement. This
   ! primary output effects the overall dynamic response of the system.  However, this routine 
   ! must also compute:
   !
   !-------------------------------------------------------------------------------------------------


         ! Passed Variables:
         
      REAL(ReKi),              INTENT( IN    )  :: X        (6)                  ! The 3 components of the translational displacement         (in m) of        the platform reference and the 3 components of the rotational displacement             (in rad) of        the platform relative to the inertial frame.
      REAL(ReKi),              INTENT( IN    )  :: ZTime                         ! Current simulation time, sec.

      REAL(ReKi),              INTENT(   OUT )  :: F        (6)                  ! The 3 components of the total force from all mooring lines (in N) acting at the platform reference and the 3 components of the total moment from all mooring lines (in N-m) acting at the platform reference; positive forces are in the direction of positive platform displacement.

      TYPE(UserLine_DataType), INTENT( INOUT )  :: UserLine_Data                 ! Internal data storage
      INTEGER,                 INTENT(   OUT )  :: ErrStat                       ! A non-zero number if an error occurs


         ! Local Variables:

      REAL(ReKi)                                :: F0       (6)                  ! Total mooring line load acting on the support platform at the reference position (N, N-m)
      REAL(ReKi)                                :: Stff     (6,6)                ! Linear restoring matrix from all mooring lines (kg/s^2, kg-m/s^2, kg-m^2/s^2)
      REAL(ReKi)                                :: X0       (6)                  ! The reference position (m, rad)

      INTEGER                                   :: I                             ! Generic index.
      INTEGER                                   :: J                             ! Generic index.


      ErrStat = 0
      

         !-------------------------------------------------------------------------------------------         
         ! Some trivial values:
         !-------------------------------------------------------------------------------------------         
      F0  (1  ) = 0.0
      F0  (2  ) = 0.0
      F0  (3  ) = 0.0
      F0  (4  ) = 0.0
      F0  (5  ) = 0.0
      F0  (6  ) = 0.0
 
      X0  (1  ) = 0.0
      X0  (2  ) = 0.0
      X0  (3  ) = 0.0
      X0  (4  ) = 0.0
      X0  (5  ) = 0.0
      X0  (6  ) = 0.0

      Stff(1,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      Stff(2,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      Stff(3,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      Stff(4,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      Stff(5,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      Stff(6,:) = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
      

      
      
!JASON:         ! Use these values for the ITI Barge:
!JASON:      F0  (1  ) =         0.0
!JASON:      F0  (2  ) =         0.0
!JASON:      F0  (3  ) =  -2058000.0
!JASON:      F0  (4  ) =         0.0
!JASON:      F0  (5  ) =         0.0
!JASON:      F0  (6  ) =         0.0
!JASON:      Stff(1,:) = (/  15920.0,       0.0,     0.0,        0.0,   208200.0,        0.0 /)
!JASON:      Stff(2,:) = (/      0.0,   15920.0,     0.0,  -208200.0,        0.0,        0.0 /)
!JASON:      Stff(3,:) = (/      0.0,       0.0, 24930.0,        0.0,        0.0,        0.0 /)
!JASON:      Stff(4,:) = (/      0.0, -208300.0,     0.0, 26230000.0,        0.0,        0.0 /)
!JASON:      Stff(5,:) = (/ 208300.0,       0.0,     0.0,        0.0, 26230000.0,        0.0 /)
!JASON:      Stff(6,:) = (/      0.0,       0.0,     0.0,        0.0,        0.0, 24730000.0 /)
!JASON:      
!JASON:         ! Use these values for the OC3-Hywind Spar:
!JASON:      F0  (1  ) =        0.0
!JASON:      F0  (2  ) =        0.0
!JASON:      F0  (3  ) = -1607000.0
!JASON:      F0  (4  ) =        0.0
!JASON:      F0  (5  ) =        0.0
!JASON:      F0  (6  ) =        0.0
!JASON:      Stff(1,:) = (/    41180.0,       0.0,     0.0,         0.0,  -2821000.0,        0.0 /)
!JASON:      Stff(2,:) = (/        0.0,   41180.0,     0.0,   2821000.0,         0.0,        0.0 /)
!JASON:      Stff(3,:) = (/        0.0,       0.0, 11940.0,         0.0,         0.0,        0.0 /)
!JASON:      Stff(4,:) = (/        0.0, 2816000.0,     0.0, 311100000.0,         0.0,        0.0 /)
!JASON:      Stff(5,:) = (/ -2816000.0,       0.0,     0.0,         0.0, 311100000.0,        0.0 /)
!JASON:      Stff(6,:) = (/        0.0,       0.0,     0.0,         0.0,         0.0, 11560000.0 /)
!JASON:      
!JASON:         ! Use these values for the MIT/NREL TLP:
!JASON:      F0  (1  ) =         0.0
!JASON:      F0  (2  ) =         0.0
!JASON:      F0  (3  ) = -30690000.0
!JASON:      F0  (4  ) =         0.0
!JASON:      F0  (5  ) =         0.0
!JASON:      F0  (6  ) =         0.0
!JASON:      Stff(1,:) = (/   197700.0,       0.0,        0.0,           0.0,    -9468000.0,         0.0 /)
!JASON:      Stff(2,:) = (/        0.0,  197700.0,        0.0,     9468000.0,           0.0,         0.0 /)
!JASON:      Stff(3,:) = (/        0.0,       0.0, 79090000.0,           0.0,           0.0,         0.0 /)
!JASON:      Stff(4,:) = (/        0.0, 9469000.0,        0.0, 30750000000.0,           0.0,         0.0 /)
!JASON:      Stff(5,:) = (/ -9469000.0,       0.0,        0.0,           0.0, 30750000000.0,         0.0 /)
!JASON:      Stff(6,:) = (/        0.0,       0.0,        0.0,           0.0,           0.0, 144100000.0 /)
      
      
         !-------------------------------------------------------------------------------------------         
         ! Calculate the force
         !-------------------------------------------------------------------------------------------         
      

      DO I = 1,6     ! Loop through all mooring line forces and moments
         F(I) = F0(I)
         DO J = 1,6  ! Loop through all platform DOFs
            F(I) = F (I) - Stff(I,J)*( X(J) - X0(J) )
         ENDDO       ! J - All platform DOFs
      ENDDO          ! I - All mooring line forces and moments


      RETURN
   END SUBROUTINE UserLine
!----------------------------------------------------------------------------------------------------
   SUBROUTINE UserLine_GetOutputs( LineNum, LineNodes, FairHTen, FairVTen, &
                                  AnchHTen, AnchVTen , Nodesxi , Nodesyi , &
                                  Nodeszi , NodesTen,  UserLine_Data, ErrStat )
                            
   ! This is a dummy routine for calculating output of a user-specified mooring system. Modify
   ! this code to create your own model of an array of mooring lines. The local variables and
   ! associated calculations below provide a template for making this user-specified mooring
   ! system model. These are provided as an example only and can be modified as desired by the
   ! user without detriment to the interface (i.e., they are not necessary for the interface).
   ! This routine must compute:
   !
   !   FairHTen    - Effective horizontal tension at the fairlead of this mooring line
   !   FairVTen    - Effective vertical   tension at the fairlead of this mooring line
   !   AnchHTen    - Effective horizontal tension at the anchor   of this mooring line
   !   AnchVTen    - Effective vertical   tension at the anchor   of this mooring line
   !
   !   NodesTen(:) - Effective line tensions              at each node of this line
   !   Nodesxi (:) - xi-coordinates in the inertial frame of each node of this line
   !   Nodesyi (:) - yi-coordinates in the inertial frame of each node of this line
   !   Nodeszi (:) - zi-coordinates in the inertial frame of each node of this line
   !
   ! These secondary outputs are only used to determine what to output for the associated 
   ! parameters placed in the OutList from the primary input file.  The number of mooring lines
   ! where the fairlead and anchor tensions can be output and the number of nodes for each line
   ! where the mooring line position and tension can be output, NumLines and LineNodes, 
   ! respectively, must be set at initialization of this module, using the 
   ! UserLine_Init() subroutine.                             
   !-------------------------------------------------------------------------------------------------
   
         ! Passed variables
         
      INTEGER,    INTENT(IN )      :: LineNodes                                   ! Number of nodes on this line where the mooring line position and tension can be output, (-).
      INTEGER,    INTENT(IN )      :: LineNum                                     ! Number of current mooring line, where the fairlead and anchor tensions can be output, (-).

   
      REAL(ReKi), INTENT(OUT)      :: AnchHTen                                    ! Effective horizontal tension at the anchor   of this mooring line, N.
      REAL(ReKi), INTENT(OUT)      :: AnchVTen                                    ! Effective vertical   tension at the anchor   of this mooring line, N.
      REAL(ReKi), INTENT(OUT)      :: FairHTen                                    ! Effective horizontal tension at the fairlead of this mooring line, N.
      REAL(ReKi), INTENT(OUT)      :: FairVTen                                    ! Effective vertical   tension at the fairlead of this mooring line, N.
      REAL(ReKi), INTENT(OUT)      :: NodesTen (LineNodes)                        ! Effective line tensions              at each node of this line, N.
      REAL(ReKi), INTENT(OUT)      :: Nodesxi  (LineNodes)                        ! xi-coordinates in the inertial frame of each node of this line, meters.
      REAL(ReKi), INTENT(OUT)      :: Nodesyi  (LineNodes)                        ! yi-coordinates in the inertial frame of each node of this line, meters.
      REAL(ReKi), INTENT(OUT)      :: Nodeszi  (LineNodes)                        ! zi-coordinates in the inertial frame of each node of this line, meters.
   
      TYPE(UserLine_DataType), INTENT( INOUT )  :: UserLine_Data                  ! Internal data storage
   
      INTEGER,    INTENT(OUT)      :: ErrStat                                     ! a non-zero value indicates an error has occurred
   
   
         ! Local variable
      INTEGER                      :: J                                           ! Dummy loop counter
         

         !-------------------------------------------------------------------------------------------         
         ! initialize error status
         !-------------------------------------------------------------------------------------------                  
      ErrStat = 0


         !-------------------------------------------------------------------------------------------         
         ! return all the outputs associated with this line
         !-------------------------------------------------------------------------------------------         

      FairHTen = 0.0
      FairVTen = 0.0
      AnchHTen = 0.0
      AnchVTen = 0.0

      DO J = 1,LineNodes                                                          ! Loop through all nodes on this line where the line position and tension can be output

         Nodesxi (J) = 0.0
         Nodesyi (J) = 0.0
         Nodeszi (J) = 0.0
         NodesTen(J) = 0.0

      END DO                                                                      ! J - All nodes on this morring line where the line position and tension can be output

      RETURN
   
   END SUBROUTINE UserLine_GetOutputs
!----------------------------------------------------------------------------------------------------
   SUBROUTINE UserLine_Terminate (UserLine_Data, ErrStat)
   ! This termination function is used to deallocate the space for the UserLine_Data storage. It 
   ! should be called at the end of the simulation (or whenever the UserLine functionality is no
   ! longer needed for this instance of HydroDyn).
   ! If the user has opened files (e.g., to read input from a file), those files should be closed
   ! here, too.
   ! Unless there are files that could be open or the user has created new allocatable arrays, 
   ! this subroutine need not be changed.
   !-------------------------------------------------------------------------------------------------

         ! Passed variables
         
      TYPE(UserLine_DataType), INTENT( INOUT )  :: UserLine_Data                 ! Internal data storage
      INTEGER,                 INTENT(   OUT )  :: ErrStat                       ! A non-zero number if an error occurs


         ! Internal variables
         
      LOGICAL                                   :: Err                           ! determines if ANY error has occurred (while allowing us to continue deallocating arrays, etc.)         


         !-------------------------------------------------------------------------------------------         
         ! Initialize error status code
         !-------------------------------------------------------------------------------------------         
      ErrStat = 0
      Err     = .FALSE.

   
         !-------------------------------------------------------------------------------------------         
         ! close files, if applicable
         !-------------------------------------------------------------------------------------------         
         
         !-------------------------------------------------------------------------------------------         
         ! deallocate arrays in UserLine_Data
         !-------------------------------------------------------------------------------------------         
         

         !-------------------------------------------------------------------------------------------         
         ! determine if an error has occurred
         !-------------------------------------------------------------------------------------------         
         
      IF ( Err ) ErrStat = 1
      RETURN
   
   
   END SUBROUTINE UserLine_Terminate       
!----------------------------------------------------------------------------------------------------     
END MODULE UserLineModule
