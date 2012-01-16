!  TestInterface.f90 
!
!  FUNCTIONS:
!  TestInterface      - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: TestInterface
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

   program TestInterface

   use hydrodyn
   USE NWTC_Library
   USE                                 SharedDataTypes    !bjj: add this to NWTC_Library???? (right now it is just AeroDyn's source file)


   implicit none
   


      ! Variables

   integer                          :: ErrStat, i
!   character(1024)                  :: FileName, FileName2

   type(HD_DataType)                :: HydroDyn_data, hd_data2
   
   type(HD_InitDataType)            :: HydroDyn_InitData
   REAL(ReKi)                       :: CurrentTime


   TYPE(HydroConfig)                :: HD_ConfigMarkers
   TYPE(AllHydroMarkers)            :: HD_AllMarkers, HD_AllMarkers2
   TYPE(AllHydroLoads)              :: HD_AllLoads,   HD_AllLoads2


!   REAL(ReKi) :: dTheta, THETA1, THETA2, THETA3, euler(3,3), smlltrans(3,3), eang(3), stang(3)
!   integer    :: j, k, ErrStat2
!beep = .false.
!dtheta = 0.05;
!
!DO I=-10,10
!   THETA1 = I*dtheta
!   
!   DO J=-10,10
!!   DO J=0,1
!      THETA2 = J*dtheta
!   
!!      DO K=0,1
!      DO K=-10,10
!         THETA3 = K*dtheta
!
!            ! get the approx trans matrix with the 3 angles
!         CALL SmllRotTrans( 'Test1', Theta1, Theta2, Theta3, smllTrans )
!         
!            ! this should be the inverse
!         stAng = GetSmllRotAngs ( smllTrans, ErrStat )
!
!         
!            ! calculate euler 123 transformation:
!         euler(1,1) =     COS(Theta2)*COS(Theta3)
!         euler(1,2) =     COS(Theta2)*SIN(Theta3)            
!         euler(1,3) = -1.*SIN(Theta2)
!         
!         euler(2,1) = -1.*COS(Theta1)*SIN(Theta3) + SIN(Theta1)*SIN(Theta2)*COS(Theta3)
!         euler(2,2) =     COS(Theta1)*COS(Theta3) + SIN(Theta1)*SIN(Theta2)*SIN(Theta3)           
!         euler(2,3) =     SIN(Theta1)*COS(Theta2)
!
!         euler(3,1) =     SIN(Theta1)*SIN(Theta3) + COS(Theta1)*SIN(Theta2)*COS(Theta3)
!         euler(3,2) = -1.*SIN(Theta1)*COS(Theta3) + COS(Theta1)*SIN(Theta2)*SIN(Theta3)            
!         euler(3,3) =     COS(Theta1)*COS(Theta2)
!                     
!         eAng = GetSmllRotAngs ( euler, ErrStat2 )
!         
!         
!         WRITE( 10, '(3(F6.3,1X), 6(F8.4,1X), I2, 1X, I2 )' ) theta1, theta2, theta3, stang, eang, ErrStat, ErrStat2
!         
!      END DO      
!   END DO
!END DO
!
!

print *, epsilon(0.0_ReKi), SQRT( epsilon(0.0_ReKi) )
print *, epsilon(0.0_DbKi), SQRT( epsilon(0.0_DbKi) )

CALL ProgExit(1)




      ! set up initialization data
      
   HydroDyn_InitData%Gravity  = 9.81 ! m/s^2

   HD_ConfigMarkers%Substructure%Position = (/0., 0., 0./)
   
echo = .true.

CurrentTime = 0


pause

      ! Body of TestInterface
         
   HydroDyn_InitData%FileName = "D:\DATA\Fortran\IVF Projects\FAST\Release\CertTest\5MWTests\70100a_cases\NRELOffshrBsline5MW_Floating_TLP\NRELOffshrBsline5MW_HydroDyn_TLP.ipt"    
   CALL GetRoot( HydroDyn_InitData%FileName, HydroDyn_InitData%OutRootName )
   HydroDyn_InitData%OutRootName = TRIM(HydroDyn_InitData%OutRootName)//'_TEST'
   call hd_init(HydroDyn_InitData, HD_ConfigMarkers, HD_AllMarkers, HydroDyn_data, ErrStat)
   
!   IF (ErrStat /= 0) THEN  
!      CALL WrScr( 'Error initializing HydroDyn.' )      
!      call progexit( 1 )
!   END IF

PAUSE

   HydroDyn_InitData%FileName = "D:\DATA\Fortran\IVF Projects\FAST\Release\CertTest\NRELOffshrBsline5MW_HydroDyn_Monopile.dat"    
   CALL GetRoot( HydroDyn_InitData%FileName, HydroDyn_InitData%OutRootName )
   CALL hd_init(HydroDyn_InitData, HD_ConfigMarkers, HD_AllMarkers2, hd_data2, ErrStat)
   


      ! CALCULATE SOME LOADS HERE
do i = 1,4
      
   CALL HD_CalculateLoads( CurrentTime,  HD_AllMarkers,  HydroDyn_data, HD_AllLoads,  ErrStat )
PAUSE
      
   CALL HD_CalculateLoads( CurrentTime,  HD_AllMarkers2, hd_data2,      HD_AllLoads2, ErrStat )
      
PAUSE
   CurrentTime = CurrentTime + 1.2
end do
      ! close hydrodyn

   CALL HD_Terminate ( HydroDyn_data, ErrStat )
   
   IF (ErrStat /= 0) THEN  
      CALL WrScr( 'Error ending HydroDyn.' )      
   END IF
  
pause   

   CALL HD_Terminate ( hd_data2, ErrStat )


   end program TestInterface

