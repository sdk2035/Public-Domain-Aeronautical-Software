
  REAL,DIMENSION(10):: alt
  REAL,DIMENSION(10):: machno
  REAL,DIMENSION(10,10):: fnidl,wfidl
  REAL,DIMENSION(7):: alt24
  REAL,DIMENSION(13):: temp
  REAL,DIMENSION(13,7):: eprMax
  REAL,DIMENSION(10):: altFF
  REAL,DIMENSION(14):: epr
  REAL,DIMENSION(14,10):: fnMax
  REAL,DIMENSION(21):: decl
  REAL,DIMENSION(21,10):: decd
  REAL,DIMENSION(10):: mach24
  REAL,DIMENSION(9):: alt27
  REAL,DIMENSION(30,9):: thrust
  REAL,DIMENSION(8,9):: mach27
  REAL,DIMENSION(30,8,9):: tsfc
  REAL,DIMENSION(10):: machN1, machN2

  COMMON/BOEING/alt,MACHNO,FNIDL,WFIDL,ALT24,TEMP,EPRMAX,ALTFF,EPR,FNMAX, &
    DECL,DECD,MACH24,ALT27,THRUST,MACH27,TSFC,MACHN1,MACHN2

!      REAL MACHNO
!      REAL MACH24,MACH27,MACHN1,MACHN2
!      COMMON/BOEING/ALT(10),MACHNO(10),FNIDL(10,10),WFIDL(10,10)        &
!     &,ALT24(7),TEMP(13),EPRMAX(13,7),ALTFF(10)                         &
!     &,EPR(14),FNMAX(14,10)                                             &
!     &,DECL(21),DECD(21,10),MACH24(10)                                  &
!     &,ALT27(9),THRUST(30,9),MACH27(8,9),TSFC(30,8,9)                   &
!     &,MACHN1(10),MACHN2(10)




