!   THE CODE IS TO USE LATTICE BOLTZMANN METHOD TO SIMULATE
!   POISEUILLE FLOW
!   DEFINE VARIABLES


!   [1] COORDINATION

      MODULE CONST1   
      integer,parameter:: XN=300
      integer,parameter:: YN=101     
      REAL*8  X(XN),Y(YN),dim_x(XN),dim_y(YN)  
      INTEGER  SOLIDX(XN*YN),SOLIDY(XN*YN)
     &        ,SFACEX(XN*YN),SFACEY(XN*YN)     
      REAL*8  pEq(XN,YN,9),pIn(XN,YN,9),pOut(XN,YN,9)
     &       ,CU(XN,YN),SpIn(XN,YN,1)  
      REAL*8  U(XN,YN),V(XN,YN),P(XN,YN)
     &       ,dim_U(XN,YN),dim_V(XN,YN)
     &       ,ZETA(XN,YN),DX(XN),DY(YN)
     &       ,real_faceX(XN*YN),real_faceY(XN*YN)
     &       ,CPFACEX(XN*YN),CPFACEY(XN*YN)              
      REAL*8  CONVU(XN,YN),CONVV(XN,YN)
     &       ,OLDU(XN,YN),OLDV(XN,YN) 
      END MODULE CONST1
!   [2] VARIABLES
      MODULE CONST2
      INTEGER    L,TSTEP,z
     &          ,BLKLHS1,BLKRHS1,BLKTOP1,BLKBOT1
     &          ,CONT,EXTRCTIME,CONT1,EXTRCTIME1
     &          ,CONT2,EXTRCTIME2
     &          ,BLKCNTX1,BLKCNTY1
      real*8 BLKRAD1
      REAL*8  SOLIDT1   
      INTEGER  SFACET,fpx,fpy,fx,fy,real_face
      REAL*8  m,n,Delta,Chi,ubf,vbf,cub,cuf,sign2,xw,yw
      INTEGER  SOLIDT,ZT(3),ZB(3),ZR(3),ZL(3)
     &        ,CX(9),CY(9),OPP(9)
      REAL*8  CONVGENCEU,CONVGENCEV,TOL,W(9)
      DATA CX /0,1,0,-1,0,1,-1,-1,1/ !velocity vectors ei for D2Q9 (eq 2.25)
      DATA CY /0,0,1,0,-1,1,1,-1,-1/
      DATA OPP /1,4,5,2,3,8,9,6,7/
      DATA ZT /5,8,9/
      DATA ZB /3,6,7/
      DATA ZR /7,4,8/
      DATA ZL /2,6,9/
      END MODULE CONST2
!   [3] PHYSICAL CONSTANTS
       MODULE CONST3
      REAL*8  UMAX,DX_P,DT_P,RE_LB,DX_LB,TauF
     &       ,DT_LB,p0,C,Co,Cs,densityG,mu
      REAL*8  L_P,BLOK,NU_P,U_P,RE_P,RE_BLC,Ma_P,Ma_LB,BHR
     &       ,T_P,CONT_DTP
      END MODULE CONST3
!   [4] PARTICLES
      MODULE CONST4      
      integer,parameter:: sampnumb=1000000 !total particle number
      REAL*8  dp,densityP,Mg,Tg,R,lambdaG,Cunnin,f,Vp,A
     &       ,RelaxT,Stokes,BoltC,Fb,Gravity,pi,BMD,Pe
     &       ,Repx,Repy,Resx,Resy,Betax,Betay,slfcx,slfcy
      REAL*8  LeftBD,RightBD,TopBD,BotBD,Uxg,Uyg
     &       ,RhsUx,RhsUy,LhsUx,LhsUy,TopUx,TopUy,BotUx,BotUy
     &       ,Lx,Ly,Gx,Gy,Bx,By
     &       ,XPP,YPP,UxgT,UxgB,UygT,UygB,xr,yr
     &       ,sign1,U1,U2
      REAL*8  uxp(sampnumb),uyp(sampnumb),xp(sampnumb),yp(sampnumb)
     &       ,olduxp(sampnumb),olduyp(sampnumb)
     &       ,oldxp(sampnumb),oldyp(sampnumb)
     &       ,rdegree(sampnumb),depxp(sampnumb),depyp(sampnumb)
      INTEGER  np,npset,D,Lxloc,Rxloc,Byloc,Tyloc
      REAL*8  capnumb(60),capeff,ID
     &       ,Periodcunt  
      INTEGER  degree1(sampnumb),idegree(sampnumb),count(sampnumb)  
      INTEGER  denumber1(360),
      END MODULE CONST4

!===================================================================================================
      PROGRAM MAIN
!===================================================================================================
      USE CONST1
      USE CONST2
      USE CONST3
      USE CONST4
      OPEN(3,FILE='PARAMETERS.TXT')
      OPEN(4,FILE='INFO.TXT')
      OPEN(6,FILE='U.TXT')
      OPEN(7,FILE='tstep VELOCITY.TXT')
      OPEN(8,FILE='tstep Particle.TXT')
      OPEN(9,FILE='ZETA.TXT')
      OPEN(10,FILE='UVUV.TXT')
      OPEN(23,FILE='Vel & Particles at theta.TXT')
      OPEN(24,FILE='Vel & Particles at 1d5-theta.TXT')
      OPEN(28,FILE='Velocity flied.TXT')
      OPEN(39,FILE='DIVERGENT VELOCITY.TXT')
      OPEN(50,FILE='endParticle.TXT')
      OPEN(51,FILE='endVelocity.TXT')
      OPEN(52,FILE='endVortex.TXT')
      OPEN(53,FILE='depositionParticle.TXT')
      OPEN(56,FILE='Fibber1Degree.TXT')
      OPEN(311,FILE='step & capp eff.TXT')
      OPEN(312,FILE='np & capnumb.TXT')
      open(313,file='point.txt')
      open(319,file='tstep pressure.txt')
      open(320,file='cp.txt')
      open(321,file='curved velocity.txt')
!     FLUID MODEL 
!     Macroscopic properties!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      L_P = 88.0526E-6  ! Channel height (m)
      BHR = 0.5792 ! BLOCK/HEIGHT RATIO~
      BLOK = L_P*BHR !block height (m)
      U_P = 0.0365 !Physical Velocity (m/s)
      NU_P = 1.6E-5 !1.983E-5 !dynamic Viscosiy (m^2/s)
      mu = 1.983E-5 !  viscosity (kg/m-s)
      densityG = mu/NU_P !kg/m^3
      DX_P = L_P/real(YN-2)   ! (m) uniform grid size in dimensional form
      RE_BLC = U_P*BLOK/NU_P ! Re based on block
      RE_P = U_P*L_P/NU_P ! Re based on Channel

      pi = 3.141592653589793     
!   Microscopic properties
      C = 1.0  ! C = 1 (always) is the lattice speed
      cs = c/3**0.5 ! speed of sound (2.21)
      DX_LB = DX_P/L_P ! uniform grid size in LATTICE BOLZMANN UNITS
      DT_LB = DX_LB/C ! uniform time step size in LATTICE BOLZMANN UNITS
      do I=1,100000  ! find the characteristic sound speed that is suitable for converting dimensional to LBM units.
            Co = U_P*I !(m/s) 
            Ma_P = U_P/(cs*Co) !physical Mach no. 
            TauF = 0.5 + 3.0*( NU_P/(Co*L_P) )/(C**2*DT_LB) 
            if ((Ma_P<0.15) .and. (TauF<1)) then
                  WRITE(3,*) "auto I=",I
                  exit
            end if
      end do      
      UMAX = U_P/Co      !inlet velocity in LBM unit
      DT_P = L_P/Co*DT_LB !(s)  Physical time step   
      RE_LB = UMAX*DX_LB*real(YN-2)/( NU_P/(Co*L_P) ) !Reynolds number in LB unit
      Ma_LB = UMAX/CS !LB mach no.
      p0 = 1.d0/3.d0
      W(1) = 4.d0/9.d0 ! (2.26) in thesis
      W(2:5) = 1.d0/9.d0 ! (2.26) in thesis
      W(6:9) = 1.d0/36.d0 ! (2.26) in thesis
!   Block Parameters
      BLKRAD1 = real(YN-2)*BHR*0.5 !the radius of the cylinder 
      BLKCNTY1 =  (YN-2)*0.5 + 2 !the j index of the cylinder center
      BLKCNTX1 =  (XN-1)*0.5  !the i index of the cylinder center
      BLKLHS1 = BLKCNTX1 - BLKRAD1
      BLKRHS1 = BLKCNTX1 + BLKRAD1
      BLKBOT1 = BLKCNTY1 - BLKRAD1
      BLKTOP1 = BLKCNTY1 + BLKRAD1
      SOLIDT = 0
! CYLINDRICAL   
      DO I = 2,XN-1
            DO J  = 2,YN-1         
!                 1st CYLINDER  
                  SOLIDT1 = (I-BLKCNTX1)**2 + (J-BLKCNTY1)**2
                  IF ( SOLIDT1 <= BLKRAD1**2  )  THEN
                        SOLIDT = SOLIDT + 1 !number of nodes in the solid 
                        SOLIDX(SOLIDT) = I !the I index of nodes in the solid 
                        SOLIDY(SOLIDT) = J !the J index of nodes in the solid 
                  END IF
            END DO
      END DO
      write(313,*) 'VARIABLES = "X","Y","p","c"'
      write(313,*) 'ZONE T="present"'
      write(313,*) 'I=',SOLIDT,', J=',1,',F=POINT'
      do i=1,SOLIDT
            write(313,*) (SOLIDX(i)-1.0)/(YN-2),(SOLIDY(i)-1.5)/(YN-2),0,0
      end do
      SFACET=0
!     1st BALL SURFACE
      do I = 1,1440
            fx=int(BLKRAD1*cos(I*pi/720))+BLKCNTX1
            fy=int(BLKRAD1*sin(I*pi/720))+BLKCNTY1
            SFACET = SFACET + 1    !number of points on the cylinder surface       
            SFACEX(SFACET) = fx     !the I index of points on the cylinder surface 
            SFACEY(SFACET) = fy     !the J index of points on the cylinder surface
      end do    
      write(313,*) 'VARIABLES = "X","Y","p","c"'
      write(313,*) 'ZONE T="present"'
      write(313,*) 'I=',SFACET,', J=',1,',F=POINT'
      do i=1,SFACET
            write(313,*) (SFACEX(i)-1.0)/(YN-2),(SFACEY(i)-1.5)/(YN-2),0,0
      end do
!     Particle parameters!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      dp = 10E-9 ! particle diameter (m)
      ID = dp/BLOK     !20E-6,10E-6,5E-6,1E-6,50E-9
      np = 0          !501 Particle number
      npset = 50      !501 Particle number in one set
      densityP = 10490.0 ! kg/m^3
      Mg = 28.966     ! AIR MOLECULAR WEIGHT (kg/kmol)
      Tg = 300.0        ! gas temperature (K)
      R = 8.31451     ! kJ/kmol-K
      lambdaG = NU_P*(pi*Mg/2.0/R/Tg/1000.0)**0.5    ! (m) molecule mean free path (eq 2.55)
      Cunnin = 1.0 + 2.0*lambdaG/dp* !cunningham correction factor (eq 2.54)
     &      (1.252 + 0.399*exp(-1.1*dp/2.0/lambdaG))
      f = 3.0*pi*mu*dp/Cunnin    !(kg/s) drag force terms (eq 2.53)
      Vp = 4.0/3.0*pi*(dp/2.0)**3    ! particle volume (m^3)
      A = densityP*Vp      !mass of particle (kg)
      RelaxT = A/f     ! Relaxation time(s)   
!     Stokes = RelaxT*U_P/BLOK
      Stokes=Cunnin*(dp**2)*(densityP/densityG)* !stokes number
     &  (RE_P*BHR)/(18.0*(BLOK**2))
!     BROWNIAN 
      BoltC = 1.38E-23 !J/K Boltzmann constant (eq 2.57)
      Fb = (6.0*pi*mu*dp*BoltC*Tg/Cunnin/DT_P)**0.5 !Brounign motion terms (eq 2.57)
      BMD = (BoltC*Tg*Cunnin)/(3.0*pi*mu*dp) !(eq 2.80)
      Pe = U_P*BLOK/BMD !Peclet number (eq 2.80)
!     Gravigy 
      Gravity = (densityP-densityG)*Vp*(-9.8) !Gravity terms (eq 2.56)
!     Gravity = 0
      CALL INPUT
      call random_seed()
!     dimensional location
      dim_x(1:XN) = X(1:XN)*L_P !physical x 
      dim_y(1:YN) = Y(1:YN)*L_P !physical y
!     Collision BC
!     CHANNEL WALLS
      LeftBD = dim_x(1) + dp/2.0 !the particle can go in x direction of the compusational domain
      RightBD = dim_x(XN) - dp/2.0
      TopBD = dim_y(YN)
      BotBD = dim_y(1)
!     OBSTACLES

  


!     Initial Setting
!     TO EXTRACT FLUID & PARTICLE DATA
      CONT = 0
      EXTRCTIME = 50000   
      T_P = 0
      CONT_DTP = 0.0    
!     to extract steady timestep
      CONT1 = 0  
      EXTRCTIME1 = 5000   
      CONVGENCEU = 1
      CONVGENCEV = 1
!     TO EXTRACT PARTICLE TRACE
      CONT2 = 0
      EXTRCTIME2 = 1 !100
!     TO EXTRACT PERIODICAL DATA
      Periodcunt = 1.0
!     ITERATION, TOL, CONSTS
      TSTEP = 100000000 
      TOL = 1E-6 !residual criteria
      WRITE(3,*) "uIn=",UMAX
      WRITE(3,*) "TSTEP=",TSTEP
      WRITE(3,*) "TAUF=",TauF
      WRITE(3,*) "Dimensional Channel height (m) =",L_P
      WRITE(3,*) "Dimensional Inlet Flow Speed (m/s) =",U_P
      WRITE(3,*) "Physical Time Step (s) =",DT_P
      WRITE(3,*) "dimensional sound speed Co=",Co
      WRITE(3,*) "dimensionless sound speed C=",C
      WRITE(3,*) "RE_CHANNEL=",RE_P
      WRITE(3,*) "RE_CHANNEL_LB=",RE_LB
      WRITE(3,*) "RE_BLOCK=",RE_BLC
      WRITE(3,*) "----- PE_paticle ----- =",Pe
      WRITE(3,*) "Normalized BLKLHS=",REAL(BLKLHS1)/REAL((YN-2))
      WRITE(3,*) "Normalized BLKRHS=",REAL(BLKRHS1)/REAL((YN-2))
      WRITE(3,*) "Normalized BLKBOT=",REAL(BLKBOT1)/REAL((YN-2))
      WRITE(3,*) "Normalized BLKTOP=",REAL(BLKTOP1)/REAL((YN-2))
      WRITE(3,*) "Physical Mach No.=",Ma_P
      WRITE(3,*) "LB Mach No.=",Ma_LB
      WRITE(3,*) "B/H=", (REAL(BLKTOP1)/REAL((YN-2))-
     &        REAL(BLKBOT1)/REAL((YN-2)) )/1.0
      WRITE(3,*) "GRID",XN,"X",YN
      WRITE(3,*) "Particle diameter =",dp
      WRITE(3,*) "dim_x(1) + dp/2 =",LeftBD
      WRITE(3,*) "particle release position =",dim_x(3)
      WRITE(3,*) "dim_x(3)-dim_x(2) =",dim_x(3)-dim_x(2)
      WRITE(3,*) "Should be larger than 1 =",dim_x(3)/LeftBD
      WRITE(3,*) "interception number (R) =",ID
      WRITE(3,*) "Particle Relaxation Time =",RelaxT
      WRITE(3,*) "Stokes Number=",Stokes
      WRITE(3,*) "DT_LB=",DT_LB
!    STOP
    DO L = 1,TSTEP
            IF ( (CONVGENCEU>=TOL) .or. (CONVGENCEV>=TOL) )THEN !!!steady state
!           if((L<TSTEP))then
				CALL FLUID    
				CONT1 = CONT1 + 1      
			end if
!      particle!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			IF ( (CONVGENCEU<TOL) .AND. (CONVGENCEV<TOL) )THEN !!!steady state
!           IF ( (L>=1) ) THEN
!           IF ( (L>TSTEP) ) THEN
				CONT_DTP = CONT_DTP + DT_P !count time when released particle
!           PARTICLE SETS RELEASED
                IF ( (CONT_DTP>=5*DT_P) .and. (np<sampnumb) )THEN 
                    np = np + npset
!                   Initial Particle velocity for a new particle set
                    uxp( (np-npset+1) : np ) = U_P  !velocity in x direction of particle
                    uyp( (np-npset+1) : np ) = 0.000 !velocity in y direction of particle
!                   Initial Particle position for a new particle set
                    xp( (np-npset+1) : np) = dim_x(3) !physcial x location of particle
                    call random_number( yp( (np-npset+1):(np) ) ) !physical y location of particle
                    yp( (np-npset+1):(np) ) = yp( (np-npset+1):(np) )*L_P
                    CONT_DTP = 0.0
                END IF         
                CALL PARTICLES
            END IF
!           particle!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            CONVU(2:XN-1,2:YN-1) = 
            &     (U(2:XN-1,2:YN-1)- OLDU(2:XN-1,2:YN-1))/
            &                        (MAXVAL(U)-MINVAL(U))
            CONVGENCEU = MAXVAL(DABS(CONVU))
            CONVV(2:XN-1,2:YN-1) = 
            &     (V(2:XN-1,2:YN-1)- OLDV(2:XN-1,2:YN-1))/
            &                        (MAXVAL(V)-MINVAL(V))
            CONVGENCEV = MAXVAL(DABS(CONVV))
            IF ((ABS(MAXVAL(U))>10) .OR. (ABS(MAXVAL(V))>10))THEN
                WRITE(4,*) 'DIVERGENT, RESET THE PARAMETERS!!!'
                WRITE(4,*) 'MAXU=',MAXVAL(U)
                WRITE(4,*) 'MAXV=',MAXVAL(V)
                WRITE(4,*) 'L=',L
                DO I = 1,XN
                    DO J = 1,YN
                         WRITE(39,7) X(I),Y(J),U(I,J),V(I,J)
                    END DO
                END DO
                STOP
            END IF
            IF ((L==TSTEP))THEN
                WRITE(4,*) 'CONVERGENCE FOR U,V,T HASNT BEEN ENOUGH !!!'
                WRITE(4,*) 'L=',L
                WRITE(4,*) 'CONVGENCEU=',CONVGENCEU
                WRITE(4,*) 'CONVGENCEV=',CONVGENCEV
                CALL VORTNSTM
            END IF
			T_P = T_P + DT_P
			CONT = CONT + 1
			if (CONT==EXTRCTIME) then      
                close(6)
                close(7)
                close(8)
                close(319)
                OPEN(6,FILE='U.TXT')
                OPEN(7,FILE='tstep VELOCITY.TXT')
                OPEN(8,FILE='tstep Particle.TXT')
                open(319,file='tstep pressure.txt')
                write(7,*) 'VARIABLES = "X","Y","U","V"'
                 write(7,*) 'ZONE T="present"'
                  write(7,*) 'I=',YN,', J=',XN,',F=POINT'
                  DO I = 1,XN
                        DO J = 1,YN
                              WRITE(7,*) X(I),Y(J),U(I,J)*Co/U_P,V(I,J)*Co/U_P
                        END DO
                  END DO
                  write(8,*) 'VARIABLES = "X","Y","P","C"'
                  write(8,*) 'ZONE T="present"'
                  write(8,*) 'I=',sampnumb,', J=',1,',F=POINT'
                  DO I = 1,np
                        WRITE(8,*) xp(I)/L_P,yp(I)/L_P,0,0
                  END DO
                  write(319,*) 'VARIABLES = "X","Y","P","C"'
                  write(319,*) 'ZONE T="present"'
                  write(319,*) 'I=',YN,', J=',XN,',F=POINT'
                  DO I = 1,XN
                        DO J = 1,YN
                              WRITE(319,*) X(I),Y(J),P(I,J),0
                        END DO
                  END DO
                  write(6,*) 'VARIABLES = "Y","U"'
                  write(6,*) 'ZONE T="present"'
                  write(6,*) 'I=',YN,', J=',1,',F=POINT'
                  DO J = 1,YN
                        WRITE(6,*) Y(J),U(XN/2,J)*Co/U_P
                  END DO
                  write(6,*) 'VARIABLES = "Y","U"'
                  write(6,*) 'ZONE T="present"'
                  write(6,*) 'I=',YN,', J=',1,',F=POINT'
                  DO J = 1,YN
                        WRITE(6,*) Y(J),U(XN/3,J)*Co/U_P
                END DO
                write(6,*) 'VARIABLES = "Y","U"'
                write(6,*) 'ZONE T="present"'
                write(6,*) 'I=',YN,', J=',1,',F=POINT'
                DO J = 1,YN
                    WRITE(6,*) Y(J),U(XN*2/3,J)*Co/U_P
                END DO
                if(capnumb(2)>0)then
                    capeff = (  capnumb(3) )/(capnumb(3)+capnumb(2))
                    write(311,*)CONT1,L,capeff,CONVGENCEU,CONVGENCEV
                    write(312,*)np,capnumb(1),capnumb(2),capnumb(3)
                end if
                CONT = 0
		end if
        OLDU(1:XN,1:YN) = U(1:XN,1:YN)
        OLDV(1:XN,1:YN) = V(1:XN,1:YN)
    END DO 
    call Cp_dimensionless_number
1     FORMAT(3X,F10.5)
7     FORMAT(3X,F10.5,3X,F10.5,3X,F10.5,3X,F10.5,3X,F10.3,3X,F10.3)
8     FORMAT(3X,F10.5,3X,F10.5,3X,F10.3,5X,F10.3,5X,F10.3,5X,F10.3)
      END PROGRAM MAIN
!===================================================================================================
      SUBROUTINE INPUT 
!===================================================================================================
!     INITIAL VELOCITY, DENSITY AND PRESSURE IN THE 2D CHANNEL
      USE CONST1
      USE CONST2
      USE CONST3
      USE CONST4  
      U(1:XN,1:YN) = 0.0000
      V(1:XN,1:YN) = 0.0000
      P(1:XN,1:YN) = p0
      OLDU(1:XN,1:YN) = 0.0000
      OLDV(1:XN,1:YN) = 0.0000
      DO I = 1,9    !2.21
            CU(1:XN,1:YN) = 3.0*( CX(I)*U(1:XN,1:YN) + CY(I)*V(1:XN,1:YN) ) ! 3*ei*u (eq 2.21)
            pIn(1:XN,1:YN,I) = W(I)*P(1:XN,1:YN)            !f^(eq)_(i) (eq 2.21)
            &           *( 1 + CU(1:XN,1:YN)/C + 0.5*CU(1:XN,1:YN)**2/C**2
            &           - 1.5*(U(1:XN,1:YN)**2 + V(1:XN,1:YN)**2)/C**2 )
      END DO
      Y(1) = 0
      Y(2) = 0.5
      Y(YN) = YN-2
      DO J = 2,YN-2 
            Y(J+1) = Y(J) +1
      ENDDO 
      X(1) = 0
      DO I = 1,XN-1 
            X(I+1) = X(I) +1
      ENDDO 
      Y(1:YN) = Y(1:YN)/real(YN-2)
      X(1:XN) = X(1:XN)/real(YN-2)
      END SUBROUTINE INPUT
!===================================================================================================
      SUBROUTINE MacroBC 
!===================================================================================================
      USE CONST1
      USE CONST2
      USE CONST3
      USE CONST4
      DO I = 1,SOLIDT
            U( SOLIDX(I),SOLIDY(I) ) = 0.000
            V( SOLIDX(I),SOLIDY(I) ) = 0.000
      end do
!     MACROSCOPIC BOUNDARY CONDITIONS 
!     Left Inlet
      U(1,1:YN) = UMAX
      V(1,1:YN) = 0.000
!     Right Inlet
      U(XN,1:YN) = 4.0/3.0*U(XN-1,1:YN)-1.0/3.0*U(XN-2,1:YN) !outlet boundary eq(2.38)
      V(XN,1:YN) = 4.0/3.0*V(XN-1,1:YN)-1.0/3.0*V(XN-2,1:YN) !outlet boundary eq(2.39)
      END SUBROUTINE MacroBC
!===================================================================================================
      SUBROUTINE MicroBC 
!===================================================================================================
      USE CONST1
      USE CONST2
      USE CONST3
      USE CONST4
!           INLET
      P(1,1:YN)=( pin(1,1:YN,1)+pin(1,1:YN,3)+pin(1,1:YN,5)       !chapter 3 of On pressure and velocity flow boundary conditions and bounceback for the lattice Boltzmann BGK model, Zou and He, 1997
      &          +2*(pin(1,1:YN,7)+pin(1,1:YN,4)+pin(1,1:YN,8)) )
      &             /((1-U(1,1:YN))/C)
      pIn(1,1:YN,2) = pIn(1,1:YN,4)+2.0/3.0*P(1,1:YN)/C*U(1,1:YN)    
      pIn(1,1:YN,6) = pIn(1,1:YN,8)+ 
      &                  ( pIn(1,1:YN,5)-pIn(1,1:YN,3))/2.0
      &                    + P(1,1:YN)*U(1,1:YN)/6.0/C
      &                    + P(1,1:YN)*V(1,1:YN)/2.0/C
      pIn(1,1:YN,9) = pIn(1,1:YN,7)+ 
      &                  ( pIn(1,1:YN,3)-pIn(1,1:YN,5))/2.0
      &                    + P(1,1:YN)*U(1,1:YN)/6.0/C
      &                    - P(1,1:YN)*V(1,1:YN)/2.0/C
!           OUTLET
      P(XN,1:YN)=( pin(XN,1:YN,1)+pin(XN,1:YN,3)+pin(XN,1:YN,5)
      &          +2*(pin(XN,1:YN,6)+pin(XN,1:YN,2)+pin(XN,1:YN,9)) )
      &             /((1+U(XN,1:YN))/C)
      pIn(XN,1:YN,4) = pIn(XN,1:YN,2)-2.0/3.0*P(XN,1:YN)/C*U(XN,1:YN)     
      pIn(XN,1:YN,8) = pIn(XN,1:YN,6)- 
      &                  ( pIn(XN,1:YN,5)-pIn(XN,1:YN,3))/2.0
      &                    - P(XN,1:YN)*U(XN,1:YN)/6.0/C
      &                    - P(XN,1:YN)*V(XN,1:YN)/2.0/C
      pIn(XN,1:YN,7) = pIn(XN,1:YN,9)- 
      &                  ( pIn(XN,1:YN,3)-pIn(XN,1:YN,5))/2.0
      &                    - P(XN,1:YN)*U(XN,1:YN)/6.0/C
      &                    + P(XN,1:YN)*V(XN,1:YN)/2.0/C
      END SUBROUTINE MicroBC

!===================================================================================================
      SUBROUTINE ConerBC 
!===================================================================================================
      USE CONST1
      USE CONST2
      USE CONST3
      USE CONST4
! !     Left Top Coner
!         pIn(1,YN,2) = pIn(1,YN,4)
!         pIn(1,YN,5) = pIn(1,YN,3)
!         pIn(1,YN,9) = pIn(1,YN,7)
!         pIn(1,YN,6) = 0.5*(P(1,YN)-pIn(1,YN,1))-
!      &                     (pIn(1,YN,4)+pIn(1,YN,7)+pIn(1,YN,3))
!             pIn(1,YN,8) = pIn(1,YN,6)
! !     Left Bot Coner
!         pIn(1,1,2) = pIn(1,1,4)
!         pIn(1,1,3) = pIn(1,1,5)
!         pIn(1,1,6) = pIn(1,1,8)
!         pIn(1,1,7) = 0.5*(P(1,1)-pIn(1,1,1))-
!      &                     (pIn(1,1,4)+pIn(1,1,8)+pIn(1,1,5))
!             pIn(1,1,9) = pIn(1,1,7)
! !     Right Top Coner
!         pIn(XN,YN,4) = pIn(XN,YN,2)
!         pIn(XN,YN,5) = pIn(XN,YN,3)
!         pIn(XN,YN,8) = pIn(XN,YN,6)
!         pIn(XN,YN,7) = 0.5*(P(XN,YN)-pIn(XN,YN,1))-
!      &                     (pIn(XN,YN,3)+pIn(XN,YN,6)+pIn(XN,YN,2))
!             pIn(XN,YN,9) = pIn(XN,YN,7)
! !     Right Bot Coner
!         pIn(XN,1,4) = pIn(XN,1,2)
!         pIn(XN,1,3) = pIn(XN,1,5)
!         pIn(XN,1,7) = pIn(XN,1,9)
!         pIn(XN,1,6) = 0.5*(P(XN,1)-pIn(XN,1,1))-
!      &                     (pIn(XN,1,2)+pIn(XN,1,9)+pIn(XN,1,5))
!           pIn(XN,1,8) = pIn(XN,1,6)
      END SUBROUTINE ConerBC

!===================================================================================================
      SUBROUTINE FLUID 
!===================================================================================================
      USE CONST1
      USE CONST2
      USE CONST3
      USE CONST4
      CALL MacroBC
      CALL MicroBC
      CALL ConerBC
!     FLUID COLLISION STEP
      DO J = 1,9
            CU(1:XN,1:YN) = 3.0*(CX(J)*U(1:XN,1:YN)+CY(J)*V(1:XN,1:YN))
            pEq(1:XN,1:YN,J) = W(J)*P(1:XN,1:YN)                              !eq(2.21)
            &     *( 1 + CU(1:XN,1:YN)/C + 0.5*CU(1:XN,1:YN)**2/C**2
            &       - 1.5*(U(1:XN,1:YN)**2 + V(1:XN,1:YN)**2)/C**2 ) 
            pOut(1:XN,1:YN,J) = pIn(1:XN,1:YN,J)-                             !eq(2.13)
            &               (1.0/TauF)*(pIn(1:XN,1:YN,J)-pEq(1:XN,1:YN,J))
      END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    REFRENCE : Lattice Boltzmann method for 3-D  Flows with Curved BOundary
      sign2=1
      real_face=0      
      DO I = 1,SFACET
            m = SFACEX(I)-BLKCNTX1
            n = SFACEY(I)-BLKCNTY1
!           The nodes are nearby the wall surface.(the solid part)
            fpx=sign(sign2,m)
            fpy=sign(sign2,n)
            if((m<BLKRAD1*0.15) .and. (m>-BLKRAD1*0.15)) then
                  fpx=0
            else if((n<BLKRAD1*0.15) .and. (n>-BLKRAD1*0.15))then
                  fpy=0
            end if
            do j=1,100
                  xw=SFACEX(I)+0.01*j*fpx
                  yw=SFACEY(I)+0.01*j*fpy
                  if (((xw-BLKCNTX1)**2+(yw-BLKCNTY1)**2)>=BLKRAD1**2) then
                        real_face=real_face+1
                        real_faceX(real_face)=XW
                        real_faceY(real_face)=YW
                        exit
                  end if
            end do    
            Delta = ( (SFACEX(I)+fpx-xw)**2+(SFACEY(I)+fpy-yw)**2 )**0.5 !eq(2.46)
            &       /( fpx**2+fpy**2 )**0.5
            IF (Delta >= 1.0/2.0) THEN    !eq(2.49)
                  Chi = (2.0*Delta-1.0)/(TauF+1.0/2.0)			!different with thesis?
                  Ubf = (1.0-3.0/(2.0*delta))*U(SFACEX(I)+fpx,SFACEY(I)+fpy)
                  Vbf = (1.0-3.0/(2.0*delta))*V(SFACEX(I)+fpx,SFACEY(I)+fpy)
            ElSE                                            !eq(2.50)
                  Chi = (2.0*Delta-1.0)/(TauF-2.0)    
                  Ubf = U(SFACEX(I)+fpx*2,SFACEY(I)+fpy*2)
                  Vbf = V(SFACEX(I)+fpx*2,SFACEY(I)+fpy*2)
            END IF
            DO J = 1,9
                  CUb = 3.0*(CX(J)*Ubf + CY(J)*Vbf)
                  CUf = 3.0*(CX(J)*U(SFACEX(I)+fpx,SFACEY(I)+fpy)
                  &         + CY(J)*V(SFACEX(I)+fpx,SFACEY(I)+fpy))
                  pEq(SFACEX(I),SFACEY(I),J) = W(J)*P(SFACEX(I)+fpx,SFACEY(I)+fpy)  !eq(2.48)
                  &      *( 1.0 + CUb/C + 0.5 *CUf*CUf/C**2
                  &     - 1.5*(U(SFACEX(I)+fpx,SFACEY(I)+fpy)**2 
                  &          + V(SFACEX(I)+fpx,SFACEY(I)+fpy)**2)/C**2  )
                  pOut(SFACEX(I),SFACEY(I),opp(J)) =								!eq(2.47)
                  &              (1.0-Chi)* pOut(SFACEX(I)+fpx,SFACEY(I)+fpy,J)
                  &                 + Chi * pEq(SFACEX(I),SFACEY(I),J)
            END DO
      End Do  
      if(L==1)then
            write(313,*) 'VARIABLES = "X","Y","p","c"'
            write(313,*) 'ZONE T="present"'
            write(313,*) 'I=',real_face,', J=',1,',F=POINT'
            do I=1,real_face
                  write(313,*) (real_faceX(i)-1.0)/(YN-2)
                  &        ,(real_faceY(I)-1.5)/(YN-2),0,0
            end do
      end if
!     BOUNCE BACK 
      pIn(1:XN,YN,5) = pOut(1:XN,1,5)
      pIn(1:XN,YN,8) = pOut(1:XN,1,8)
      pIn(1:XN,YN,9) = pOut(1:XN,1,9)
      pIn(1:XN,1,7) = pOut(1:XN,YN,7)
      pIn(1:XN,1,3) = pOut(1:XN,YN,3)
      pIn(1:XN,1,6) = pOut(1:XN,YN,6)

!     FLUID STREAMING STEP  
       DO J = 1,XN                          !  C0
            DO K = 1,YN
                  pIn(J,K,1) = pOut(J,K,1)
            ENDDO
      ENDDO  
      DO J = 2,XN                         ! --C2
            DO K = 1,YN                   !
                  pIn(J,K,2) = pOut(J-1,K,2)
            ENDDO
      ENDDO  
      DO J = 1,XN                               !   C3
            DO K = 2,YN                         !   |
                  pIn(J,K,3) = pOut(J,K-1,3)    !   |
            ENDDO
      ENDDO
      DO J = 1,XN-1                             !   C4--
            DO K = 1,YN                         !
                    pIn(J,K,4) = pOut(J+1,K,4)
            ENDDO
      ENDDO
      DO J = 1,XN                             !   |
           DO K = 1,YN-1                      !   |
                  pIn(J,K,5) = pOut(J,K+1,5)  !   C5
           ENDDO
      ENDDO
      DO J = 2,XN
          DO K = 2,YN                                 !   C6
                  pIn(J,K,6) = pOut(J-1,K-1,6)        !  /
          ENDDO
      ENDDO  
      DO J = 1,XN-1
          DO K = 2,YN                     !   C7
                  pIn(J,K,7) = pOut(J+1,K-1,7) !     \ 
          ENDDO
      ENDDO
      DO J = 1,XN-1
           DO K = 1,YN-1                   !     /
                  pIn(J,K,8) = pOut(J+1,K+1,8) !   C8
           ENDDO
      ENDDO 
      DO J = 2,XN
           DO K = 1,YN-1                   !  \ 
                  pIn(J,K,9) = pOut(J-1,K+1,9) !   C9
           ENDDO
      ENDDO
!         ASSEMBLE MACROSCOPIC VARIABLES
!         density
      SpIn(1:XN,1:YN,1) = pIn(1:XN,1:YN,1)
      DO J = 2,9
           SpIn(1:XN,1:YN,1) = SpIn(1:XN,1:YN,1)
           &        + pIn(1:XN,1:YN,J)
      END DO
      P(1:XN,1:YN) = SpIn(1:XN,1:YN,1) 
!     U COMPONENT
      SpIn(1:XN,1:YN,1) = CX(1)*pIn(1:XN,1:YN,1)
      DO J = 2,9   
            SpIn(1:XN,1:YN,1) = SpIn(1:XN,1:YN,1)
            &        +CX(J)*pIn(1:XN,1:YN,J)
      END DO
      U(1:XN,1:YN) = SpIn(1:XN,1:YN,1)           
      U(1:XN,1:YN) = C*U(1:XN,1:YN)/P(1:XN,1:YN) 
!     V COMPONENT
      SpIn(1:XN,1:YN,1) = CY(1)*pIn(1:XN,1:YN,1)
      DO J = 2,9   
            SpIn(1:XN,1:YN,1) = SpIn(1:XN,1:YN,1)
            &        +CY(J)*pIn(1:XN,1:YN,J)
      END DO
      V(1:XN,1:YN) = SpIn(1:XN,1:YN,1)
      V(1:XN,1:YN) = C*V(1:XN,1:YN)/P(1:XN,1:YN) 
!       CALL MacroBC
!       CALL MicroBC
!       CALL ConerBC
      END SUBROUTINE FLUID

!===================================================================================================
      SUBROUTINE PARTICLES 
!===================================================================================================
      USE CONST1
      USE CONST2
      USE CONST3
      USE CONST4
!     PARTICLE EQUATION
!     Dimensional velocity
      capnumb(1:60) = 0.0
      count(1:sampnumb) = 0
      denumber1(1:360)=0
      D = 0
      dim_U(1:XN,1:YN) = U(1:XN,1:YN)*Co
      dim_V(1:XN,1:YN) = V(1:XN,1:YN)*Co
      DO I = 1,np  
            SOLIDT1 = ( xp(I)-dim_x(BLKCNTX1) )**2 + 
            &          ( yp(I)-dim_y(BLKCNTY1) )**2
            IF (   yp(I)< BotBD  )THEN    ! line695-700 periodic BC
                  
                  yp(I) = yp(I) + (TopBD - BotBD) !yp(I) 
            ELSEIF (  yp(I)> TopBD  )THEN
                  
                  yp(I) = yp(I) - (TopBD - BotBD) !yp(I)      
            ELSEIF  ( xp(I)<= LeftBD )  THEN
                  xp(I) = xp(I)+ dim_x(4)
                   
                  IF (I<=sampnumb)THEN   !meaningless ?
                        capnumb(1) = capnumb(1) + 1 
                  ENDIF   
            ELSEIF (  xp(I)>= RightBD   )THEN
                  
                  
                  IF (I<=sampnumb)THEN
                        capnumb(2) = capnumb(2) + 1 
                  ENDIF
!           CYLINDER 1
            ELSEIF ( SOLIDT1 <= ( (BLOK + dp)/2.0 )**2 )THEN
				IF (I<=sampnumb)THEN
					capnumb(3) = capnumb(3) + 1  !TOTAL TRAPPED PARTICLES 
				END IF   
            
            
				If (count(I) == 0) THEN  !meaningless ?
					rdegree(I)=atan2(yp(I)-dim_y(BLKCNTY1),xp(I)-dim_x(BLKCNTX1))
					idegree(I) = int(rdegree(I)/pi*180)
					if (yp(I)<dim_y(BLKCNTY1))then
                        idegree(I) = int(rdegree(I)/pi*180+360)
					end if   
					degree1(I) = idegree(I)
					denumber1(degree1(I))=denumber1(degree1(I))+1
					D = D + 1
					depxp(D) = xp(I)
					depyp(D) = yp(I)
					count(I) = 1
				END IF
            ELSE
!               Fluid velocity at RHS of particle
                XPP = xp(I) + dp/2.0
                YPP = yp(I)
                CALL PVELOCITY
                RhsUx = Uxg
                RhsUy = Uyg
!               Fluid velocity at LHS of particle
                XPP = xp(I) - dp/2.0
                YPP = yp(I)
				CALL PVELOCITY
                LhsUx = Uxg
                LhsUy = Uyg
!                Fluid velocity at Top of particle
                XPP = xp(I)
                YPP = yp(I)+ dp/2.0
                CALL PVELOCITY
                TopUx = Uxg
                TopUy = Uyg
!               Fluid velocity at Bottom of particle  
                XPP = xp(I)
                YPP = yp(I)- dp/2.0
                CALL PVELOCITY
                BotUx = Uxg
                BotUy = Uyg
!               Fluid velocity at the position of the particle
                XPP = xp(I)
                YPP = yp(I)
                CALL PVELOCITY
!==============================================================================================
!==============================================================================================
                Repx = dp/Nu_P*abs( Uxg-uxp(I))				!(eq 2.64)
                Repy = dp/Nu_P*abs( Uyg-uyp(I))
                Resx = dp**2/Nu_P*( 0.5*abs(RhsUy - LhsUy)/dp )	!(eq 2.65)
                Resy = dp**2/NU_P*( 0.5*abs(TopUx - BotUx)/dp )
                Betax = 0.5*Resx/Repx						!(eq 2.67)
                Betay = 0.5*Resy/Repy
                if (Repx<=40) THEN							!(eq 2.62)
                    slfcx = (1.0-0.3314*Betax**0.5)*exp(-0.1*Repx)+0.3314*Betax**0.5
                else
                    slfcx = 0.0524*(0.5*Resx)**0.5
                End If
                IF (Repy<=40) THEN							!(eq 2.63)
                    slfcy = (1.0-0.3314*Betay**0.5)*exp(-0.1*Repy)+0.3314*Betay**0.5
                else
                    slfcy = 0.0524*(0.5*Resy)**0.5
                End If
!==============================================================================================
!==============================================================================================
                sign1 = 1
!               Fluid lift force without velocity term (kg/s)
                Lx = 1.615*densityG*NU_P**0.5*dp**2				!(eq 2.72)
                & *( abs(RhsUy - LhsUy)/dp )**0.5
                &   *sign( sign1, (RhsUy - LhsUy)/dp )*slfcx
                Ly = 1.615*densityG*NU_P**0.5*dp**2				!(eq 2.73)
                & *( abs(TopUx - BotUx)/dp )**0.5
                &   *sign( sign1, (TopUx - BotUx)/dp )*slfcy
!               Bforce    (N)
                call random_number(U1 )
                call random_number(U2 )
                Gx = ( -2.0*log(U1) )**0.5*cos(2.0*pi*U2)		!(eq 2.58)
                Gy = ( -2.0*log(U1) )**0.5*sin(2.0*pi*U2)
                Bx = Fb*Gx
                By = Gravity + Fb*Gy      
!               Particle velocity (m/s)
                olduxp(I) = uxp(I)
                uxp(I) = (olduxp(I) - ( Uxg + Bx/(f+Lx) ) )*exp(-(f+Lx)/A*DT_P) 		!(eq 2.76)
                &       + Uxg + Bx/(f+Lx)
                olduyp(I) = uyp(I)
                uyp(I) = (olduyp(I) - ( Uyg + By/(f+Ly) ) )*exp(-(f+Ly)/A*DT_P)
                &       + Uyg + By/(f+Ly)
!               Particle position (m)
                oldxp(I) = xp(I)
                xp(I) = oldxp(I) + (Uxg + Bx/(f+Lx) )*DT_P + A/(f+Lx)*					!(eq 2.77)
                &( olduxp(I) - Uxg - Bx/(f+Lx) )*(1.0)*( 1.0-exp(-(f+Lx)*DT_P/A) )
                oldyp(I) = yp(I)
                yp(I) = oldyp(I) + (Uyg + By/(f+Ly) )*DT_P + A/(f+Ly)*
                &( olduyp(I) - Uyg - By/(f+Ly) )*(1.0)*( 1.0-exp(-(f+Ly)*DT_P/A) )
            ENDIF
		ENDDO
		IF(capnumb(3)+capnumb(2)>=sampnumb)THEN
            WRITE(4,*) 'capnumb AT left WALL=',capnumb(1)
            WRITE(4,*) 'capnumb AT right WALL=',capnumb(2)
            WRITE(4,*) 'capnumb AT CYLINDER 1',capnumb(3)
            WRITE(4,*) 'capnumb AT CYLINDER 2',capnumb(8)
            WRITE(4,*) 'capnumb AT CYLINDER 3',capnumb(13)
            WRITE(4,*) 'sampnumb=',sampnumb             
            capeff = ( capnumb(1) + capnumb(2) + capnumb(3)
            &             + capnumb(8) + capnumb(13) + capnumb(18)
            &             + capnumb(23) + capnumb(28) + capnumb(33) 
            &             + capnumb(38) + capnumb(43) + capnumb(48) )/sampnumb
            WRITE(4,*) 'CAPTURE EFF=',capeff               
            capeff = ( capnumb(1)  )/sampnumb
            WRITE(4,*) 'inlet CAPTURE EFF=',capeff              
            capeff = ( capnumb(2)  )/sampnumb
            WRITE(4,*) 'outlet CAPTURE EFF=',capeff             
            capeff = (  capnumb(3) + capnumb(8) + capnumb(13)
            &               + capnumb(18) + capnumb(23) + capnumb(28)
            &               + capnumb(33) + capnumb(38) + capnumb(43) 
            &               + capnumb(48))/sampnumb
            WRITE(4,*) 'OBST CAPTURE EFF=',capeff 
            WRITE(4,*) 'Total time step is',L
            WRITE(4,*) 'Total physical time is',T_P
            DO I = 1,XN
                  DO J = 1,YN
                        WRITE(51,4) X(I),Y(J),U(I,J),V(I,J)
                  END DO
            END DO
            DO I = 1,np
                  WRITE(50,*) xp(I)/L_P,yp(I)/L_P,0,0
            END DO
            CALL VORTNSTM
            DO I = 1,XN
                  DO J = 1,YN
                        WRITE(52,*) X(I),Y(J),ZETA(I,J)
                  END DO
            END DO             
            DO I = 1,D
                  write(53,*)depxp(I)/L_P,depyp(I)/L_P,0
            END DO           
            DO I = 1,360
                  write(56,*)I,denumber1(I)
            END DO
            call Cp_dimensionless_number
            STOP
      END IF
4     FORMAT(3X,F10.5,3X,F10.5,3X,F10.5,3X,F10.5,3X,F10.3,3X,F10.3)
      END SUBROUTINE PARTICLES
!===================================================================================================
      SUBROUTINE PVELOCITY
!===================================================================================================
      USE CONST1
      USE CONST2
      USE CONST3
      USE CONST4     
      Lxloc = int(XPP/L_P*real(YN-2))+1
      Rxloc = Lxloc + 1     
      Byloc = int(YPP/L_P*real(YN-2)+0.5)+1
      Tyloc = Byloc + 1
      xr = ( XPP - dim_x(Lxloc) )/( dim_x(Rxloc) - dim_x(Lxloc)  )
      UxgT = xr*( dim_U(Rxloc, Tyloc ) - dim_U(Lxloc, Tyloc ) )
     &     + dim_U(Lxloc, Tyloc )
      UxgB = xr*( dim_U(Rxloc, Byloc ) - dim_U(Lxloc, Byloc ) )
     &     + dim_U(Lxloc, Byloc ) 
      UygT = xr*( dim_V(Rxloc, Tyloc ) - dim_V(Lxloc, Tyloc ) )
     &     + dim_V(Lxloc, Tyloc )    
      UygB = xr*( dim_V(Rxloc, Byloc ) - dim_V(Lxloc, Byloc ) )
     &     + dim_V(Lxloc, Byloc )
      yr = ( YPP - dim_y(Byloc) )/( dim_y(Tyloc) - dim_y(Byloc)  )
      Uxg = yr*(UxgT - UxgB) + UxgB
      Uyg = yr*(UygT - UygB) + UygB
      END SUBROUTINE PVELOCITY

!===================================================================================================
      SUBROUTINE VORTNSTM 
!===================================================================================================
      USE CONST1
      USE CONST2
      USE CONST3
      USE CONST4
      DO I = 1,XN-1
            DX(I) = X(I+1) - X(I)
      END DO
      DO J = 1,YN-1
            DY(J) = Y(J+1) - Y(J)
      END DO
      DO I = 2,XN-1
            DO J = 2,YN-1
                  ZETA(I,J) = (V(I+1,J)-V(I-1,J))/(DX_LB+DX_LB)
                  &             -(U(I,J+1)-U(I,J-1))/(DX_LB+DX_LB)
            END DO
      END DO         
!     BC LHS & RHS
      DO J = 2,YN-1   
            ZETA(1,J) = (V(2,J)-V(1,J))/DX_LB
            &             -(U(1,J+1)-U(1,J-1))/(DX_LB+DX_LB)
            ZETA(XN,J) = (V(XN,J)-V(XN-1,J))/DX_LB
            &             -(U(XN,J+1)-U(XN,J-1))/(DX_LB+DX_LB)
      END DO
!     BC TOP & BOT
      DO I = 2,XN-1   
            ZETA(I,1) = (V(I+1,1)-V(I-1,1))/(DX_LB+DX_LB)
            &             -(U(I,2)-U(I,1))/DX_LB
            ZETA(I,YN) = (V(I+1,YN)-V(I-1,YN))/(DX_LB+DX_LB)
            &             -(U(I,YN)-U(I,YN-1))/DX_LB
      END DO
      DO I = 1,XN
            DO J = 1,YN
                  WRITE(28,*) X(I),Y(J),U(I,J),V(I,J)
                  WRITE(9,*) X(I),Y(J),ZETA(I,J)
            END DO
      END DO
      END SUBROUTINE VORTNSTM
      
!===================================================================================================				   				   
      SUBROUTINE Cp_dimensionless_number 
!===================================================================================================				   				   
      USE CONST1
      USE CONST2
      USE CONST3
      USE CONST4
      integer SFACET_cp1
      real*8 fx1,fy1,delta_DEGREE,PX_p(1440),PY_p(1440),cp
      do j=1,10    
            SFACET_cp1=0
            do I = 1,1440
                  fx1=(BLKRAD1*(1+0.01*j))*cos(I*pi/720)+BLKCNTX1
                  fy1=(BLKRAD1*(1+0.01*j))*sin(I*pi/720)+BLKCNTY1
                  SFACET_cp1 = SFACET_cp1 + 1
                  CPFACEX(SFACET_cp1) = fx1
                  CPFACEY(SFACET_cp1) = fy1
            end do     
            write(313,*) 'VARIABLES = "X","Y","p","c"'
            write(313,*) 'ZONE T="present"'
            write(313,*) 'I=',SFACET_cp1,', J=',1,',F=POINT'
            do I=1,SFACET_cp1
                  write(313,*) (CPFACEX(i)-1.0)/(YN-2)
                  &        ,(CPFACEY(I)-1.5)/(YN-2),0,0
            end do
            write(320,*)'VARIABLES = "degree","Cp"'
            write(320,*)'ZONE T="Cp"'
            write(320,*)'I=1, J=',SFACET_cp1,',F=POINT'     
            write(321,*)'VARIABLES = "degree","velocity"'
            write(321,*)'ZONE T="surface velocity"'
            write(321,*)'I=1, J=',SFACET_cp1,',F=POINT'     
            do I=1,SFACET_cp1
                  DO Z = 2,XN     
                        IF ( (CPFACEX(i)-1.0)/(YN-2) <= x(Z) ) THEN
                              Lxloc = Z-1
                              Rxloc = Z     
                        EXIT
                        END IF     
                  END DO
                  DO Z = 2,YN
                        IF ( (CPFACEY(I)-1.5)/(YN-2) <= y(Z) ) THEN
                              Byloc = Z-1
                              Tyloc = Z
                        EXIT
                        END IF     
                  END DO
                  xr = ( (CPFACEX(i)-1.0)/(YN-2) - x(Lxloc) )
                  &      /( x(Rxloc) - x(Lxloc)  )
                  UxgT = xr*( U(Rxloc, Tyloc ) - U(Lxloc, Tyloc ) )
                  &     + U(Lxloc, Tyloc )
                  UxgB = xr*( U(Rxloc, Byloc ) - U(Lxloc, Byloc ) )
                  &     + U(Lxloc, Byloc ) 
                  UygT = xr*( V(Rxloc, Tyloc ) - V(Lxloc, Tyloc ) )
                  &     + V(Lxloc, Tyloc )   
                  UygB = xr*( V(Rxloc, Byloc ) - V(Lxloc, Byloc ) )
                  &     + V(Lxloc, Byloc )
                  yr = ( (CPFACEY(I)-1.5)/(YN-2) - y(Byloc) )
                  &      /( y(Tyloc) - y(Byloc)  ) 
                  Uxg = yr*(UxgT - UxgB) + UxgB  
                  Uyg = yr*(UygT - UygB) + UygB     
                  Cp = 1-(Uxg**2+Uyg**2)/UMAX**2 
      !          Cp = (P(SFACEX(i),SFACEY(i))-P(2,(YN-1)/2))/3.0
      !           &         /(0.5*P(2,(YN-1)/2)*UMAX**2)      
                  PX_p(I)=CPFACEX(i)-BLKCNTX1
                  PY_p(I)=CPFACEY(i)-BLKCNTY1
                  delta_DEGREE= ATAN2(PY_p(i),PX_p(i))/PI*180
                  if (PY_p(I)<=0) then
                        delta_DEGREE= ATAN2(PY_p(i),PX_p(i))/PI*180+360
                  end if  
                  write(321,*) delta_DEGREE,(Uxg**2+Uyg**2)**0.5
                  write(320,*) delta_DEGREE,Cp
            end do
      end do
      END SUBROUTINE Cp_dimensionless_number        !pressure drop (P(3,(YN-1)/2)-P(XN-2,(YN-1)/2))*densityG*Co**2
