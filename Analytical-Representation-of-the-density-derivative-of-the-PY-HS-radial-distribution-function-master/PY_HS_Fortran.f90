!*==AA0001.spg  processed by SPAG 6.72Dc at 03:26 on  5 Apr 2020
      IMPLICIT NONE
!*--AA00013
!*** Start of declarations inserted by SPAG
      REAL*8 d , dr , DRRG , dy , flag , GR , one , PI , r , rd , rho , &
           & rhod , y
      INTEGER j , nr
!*** End of declarations inserted by SPAG
      COMMON /PIE   / PI
      
!     If you find this code useful, please cite as: 
      
!     Braden D. Kelly and William R. Smith and Douglas Henderson, 
!     Analytical representation of the density derivative of the 
!     Percus–Yevick hard-sphere radial distribution function,
!     Molecular Physics, 114,16-17,2446-2450,2016
!     https://doi.org/10.1080/00268976.2016.1164908

!     THIS PROGRAM CALCULATES BOTH Y(R), WHICH IS EQUAL TO G(R) FOR R>1.
!     AS WELL AS THE DENSITY DERIVATIVE OF Y(R). THE DISTANCES ARE
!     EXPECTED TO BE IN UNITS OF R/D AND THE DENSITY IS EXPECTED TO BE
!     IN UNITS OF RHO*(D**3) WHERE D IS THE HARD SPHERE DIAMETER. YOU
!     PROCEED IN THE FOLLOWING FASHION.
      one = 1.D0
      PI = 4.D0*DATAN(one)
      d = 1.D0
!     YOU CAN USE ANY UNITS YOU WANT BUT IF D=1 DISTANCES ARE IN UNITS OF R/D
      dr = 0.1D0
      nr = 1
      rho = 1.0D0
      rhod = rho*(d**3)
!     GRPRLM GIVES THE PY APPROXIMATION;
      CALL GRPRLM(rhod)
!
      PRINT 99001 , rho
99001 FORMAT ("Density (RHO)",F16.1)
!
!     YOU HAVE NOW SET UP THE CONSTANTS WHICH DEPEND UPON RHO BUT ARE
!     INDEPENDANT OF R.
!
!     WHEN YOU COME TO THE PART OF THE PROGRAM WHICH NEEDS Y(R) OR G(R)
!     DO SOMETHING LIKE THE FOLLOWING
      flag = 1.D0
!     IF FLAG >=1 THEN G(R) AS WELL AS THE DERIVATIVE OF G(R) WILL BE
!     CALCULATED. OTHERWISE ONLY G(R) WILL BE CALCULATED.
!
      IF ( flag.GE.1.D0 ) THEN
         DO j = 1 , 101
            r = (j-1)*dr
            rd = r/d
            y = GR(rd)
            CALL GRPRLD(rhod,rd)
            dy = DRRG(rd)
!      IF(RD.LT.1.D0)Y=0.D0
!      IF(RD.LT.1.D0)DY=0.D0
            PRINT 99002 , rd , y , dy
99002       FORMAT (F16.2,2F16.6)
         ENDDO
      ELSE
         DO j = 1 , 101
            r = (j-1)*dr
            rd = r/d
            y = GR(rd)
!      IF(RD.LT.1.D0)Y=0.D0
            PRINT 99003 , rd , y
99003       FORMAT (F16.2,F16.6)
         ENDDO
      ENDIF
      END
!*==GRPRLM.spg  processed by SPAG 6.72Dc at 03:26 on  5 Apr 2020
!
!     !%--------------------------------------------------------------%!
!
      SUBROUTINE GRPRLM(Rho)
      IMPLICIT NONE
!*--GRPRLM67
!*** Start of declarations inserted by SPAG
      REAL*8 bot1 , bot2 , C1 , C2 , C3 , ETA , eta12 , eta3 , ETAm ,   &
           & EX11 , ff , IG101 , IG201 , IG211 , IG221 , IG301 , IG311 ,&
           & IG321 , IG331 , IG401
      REAL*8 IG411 , IG421 , IG431 , IG441 , IG501 , IG511 , IG521 ,    &
           & IG531 , IG541 , IG551 , ILT , ISDp , ISP , IT1 , IXL1 ,    &
           & IXL2 , IXL3 , IXL4 , IXL5 , ixx1
      REAL*8 lp , LT1 , LT2 , par , PI , Rho , SDP , sdp1 , sdp2 , sp1 ,&
           & sp2 , sp3 , stp , top , top1 , top2 , topp , xeta , XF1 ,  &
           & XF2
      REAL*8 XF3 , XF4 , XF5 , yminus , yplus
!*** End of declarations inserted by SPAG
      COMMON /PIE   / PI
      COMMON /BH    / C1 , C2 , C3
      COMMON /RGR   / IG102 , IG202 , IG212 , IG222 , IG302 , IG312 ,   &
                    & IG322 , IG332 , IG402 , IG412 , IG422 , IG432 ,   &
                    & IG442 , IG502 , IG512 , IG522 , IG532 , IG542 ,   &
                    & IG552 , IT2 , EX12 , IG101 , IG201 , IG211 ,      &
                    & IG221 , IG301 , IG311 , IG321 , IG331 , IG401 ,   &
                    & IG411 , IG421 , IG431 , IG441 , IG501 , IG511 ,   &
                    & IG521 , IG531 , IG541 , IG551 , IT1 , EX11 , ETA
!     Common variables to pass to derivative subroutine GRPRLD
      COMMON /DERIV / ISP , XF1 , XF2 , XF3 , XF4 , XF5 , ETAm , LT1 ,  &
                    & LT2 , ILT , ISPc , ILTc , SDP , ISDp , IXL1 ,     &
                    & IXL2 , IXL3 , IXL4 , IXL5
      COMPLEX*16 EX12 , IG102 , IG202 , IG212 , IG222 , IG302 , IG312 , &
               & IG322 , IG332 , IG402 , IG412 , IG422 , IG432 , IG442 ,&
               & IG502 , IG512 , IG522 , IG532 , IG542 , IG552 , IT2
      COMPLEX*16 DCMPLX , CDEXP , ixl1c , ixl2c , ixl3c , ixl4c ,       &
               & ixl5c , ILTc , jx , ISPc , isdpc , ixx1c
!
      ETA = PI*Rho/6.D0
      ETAm = 1.D0 - ETA
      bot1 = ETAm**4
      top = 1.D0 + 2.D0*ETA
      top1 = top*top
      C1 = top1/bot1
      bot2 = 4.D0*bot1
      topp = 2.D0 + ETA
      top2 = topp*topp
      C2 = -top2/bot2
      C2 = 6.D0*ETA*C2
      C3 = ETA*C1/2.D0
      eta12 = 12.D0*ETA
      LT1 = 1.D0 + ETA/2.D0
      LT2 = 1.D0 + 2.D0*ETA
      lp = LT1
      xeta = 1.D0/(1.D0-ETA)
      eta3 = (1.D0-ETA)**2
      sp1 = 1.D0
      sp2 = 4.D0*ETA*xeta
      sp3 = 6.D0*(ETA*xeta)**2
      sdp1 = 2.D0
      sdp2 = sp2
      stp = 2.D0
      ff = (-ETA+3.D0)*ETA + 3.D0
      par = DSQRT(1.D0+2.D0*(ETA**2/ff)**2)
      yplus = (1.D0+par)**(1.D0/3.D0)
      yminus = -(-1.D0+par)**(1.D0/3.D0)
      jx = CDEXP(DCMPLX(0.D0,2.D0*PI/3.D0))
      par = (2.D0*ETA*ff)**(1.D0/3.D0)
      IT1 = xeta*(-2.D0*ETA+par*(yplus+yminus))
      IT2 = xeta*(-2.D0*ETA+par*(yplus*jx+yminus/jx))
      XF1 = 3.D0*eta3
      par = eta3**2/ETA
      XF2 = par*3.D0/4.D0
      par = par*eta3/ETA
      XF3 = par*3.D0/8.D0
      par = par*eta3/ETA
      XF4 = par*9.D0/32.D0
!     !%----------------- Variables for First Root -------------------%!
      ILT = LT1*IT1 + LT2
      ISDp = sdp1*IT1 + sdp2
      ISP = (sp1*IT1+sp2)*IT1 + sp3
      EX11 = DEXP(-IT1)
      IXL1 = (ILT/ISP)**2
      IXL2 = 15.D0*ISDp**2 - 4.D0*ISP*stp
      IXL3 = ILT + 4.D0*IT1*lp
      IXL4 = ILT*ISDp/ISP
      IXL5 = 2.D0*ILT + 3.D0*IT1*lp
      IG101 = IT1*ILT/XF1/ISP
      IG201 = -ILT/XF2/ISP**2
      IG211 = ILT*(1.D0-IT1*(2.D0+ISDp/ISP)) + 2.D0*lp*IT1
      IG221 = ILT*IT1
      IG301 = ILT/XF3/ISP**3
      IG311 = ILT**2*IT1/ISP**2*(3.D0*ISDp**2-ISP*stp)                  &
            & - 3.D0*ILT*ISDp/ISP*(ILT+3.D0*IT1*lp)                     &
            & + 6.D0*lp*(ILT+lp*IT1)
      IG321 = ILT*(6.D0*lp*IT1+ILT*(2.D0-3.D0*ISDp*IT1/ISP))
      IG331 = IG221*ILT
      IG401 = -ILT/XF4/ISP**4
      IG411 = 5.D0*IT1*IXL1*ILT/ISP*ISDp*(2.D0*ISP*stp-3.D0*ISDp**2)    &
            & + IXL1*IXL2*IXL3 - 24.D0*lp*IXL4*IXL5 +                   &
            & 12.D0*lp**2*(3.D0*ILT+2.D0*IT1*lp)
      IG421 = (IT1*IXL1*IXL2-12.D0*(IXL4*IXL3-lp*IXL5))*ILT
      IG431 = (-6.D0*IT1*IXL4+3.D0*IXL3)*ILT**2
      IG441 = IG331*ILT
!
      ISP = ISP*XF1
      ISDp = ISDp*XF1
      stp = stp*XF1
      ixx1 = ILT + 5.D0*IT1*lp
      IG501 = 864.D0*ETA**4.D0*ILT/(ISP**5)
!
      IG511 = 120.D0*LT1**3*(2.D0*ILT+LT1*IT1)                          &
            & + (10.D0*ILT**2*(6.D0*ILT**2*ISDp*stp+IT1*ILT**2*stp**2+  &
            & 45.D0*ILT*LT1*ISDp**2+30.D0*IT1*ILT*LT1*ISDp*stp+         &
            & 90.D0*IT1*LT1**2.D0*ISDp**2))                             &
            & /ISP**2 - (100.D0*ILT*LT1*(ILT**2*stp+6.D0*ILT*LT1*ISDp+  &
            & 6.D0*LT1**2*ISDp*IT1+2.D0*ILT*LT1*stp*IT1))               &
            & /ISP - (105.D0*ILT**3*ISDp**2*(ILT*ISDp+ILT*stp*IT1+      &
            & 5.D0*LT1*ISDp*IT1))/ISP**3 + (105.D0*ILT**4*ISDp**4*IT1)  &
            & /ISP**4
!
      IG521 = 240.D0*ILT*LT1**2*(ILT+LT1*IT1)                           &
            & - (20.D0*ILT**2*(ILT**2.D0*stp+15.D0*ILT*LT1*ISDp+        &
            & 30.D0*LT1**2*ISDp*IT1+5.D0*ILT*LT1*stp*IT1))              &
            & /ISP + (30.D0*ILT**3*ISDp*(3.D0*ILT*ISDp+2.D0*ILT*stp*IT1+&
            & 15.D0*LT1*ISDp*IT1))/ISP**2 - (105.D0*ILT**4*ISDp**3*IT1) &
            & /ISP**3
!
      IG531 = 5.D0*ILT**2*(24.D0*IT1*LT1**2+12.D0*ILT*LT1)              &
            & + (45.D0*ILT**4*ISDp**2*IT1-                              &
            & 5.D0*ILT**2*ISP*(6.D0*ILT**2*ISDp+2.D0*ILT**2*stp*IT1+    &
            & 30.D0*ILT*LT1*ISDp*IT1))/ISP**2
!
      IG541 = 4.D0*ILT**3*ixx1 - 10.D0*IT1*ILT**3*IXL4
      IG551 = IT1*ILT**4
!
      ISP = ISP/XF1
      ISDp = ISDp/XF1
      stp = stp/XF1
!     !%--------------- Variables for Second Root --------------------%!
      ILTc = LT1*IT2 + LT2
      isdpc = sdp1*IT2 + sdp2
      ISPc = (sp1*IT2+sp2)*IT2 + sp3
      EX12 = CDEXP((0.D0,0.D0)-IT2)
      ixl1c = (ILTc/ISPc)**2
      ixl2c = 15.D0*isdpc**2 - 4.D0*ISPc*stp
      ixl3c = ILTc + 4.D0*IT2*lp
      ixl4c = ILTc*isdpc/ISPc
      ixl5c = 2.D0*ILTc + 3.D0*IT2*lp
      IG102 = IT2*ILTc/XF1/ISPc
      IG202 = -ILTc/XF2/ISPc**2
      IG212 = ILTc*(1.D0-IT2*(2.D0+isdpc/ISPc)) + 2.D0*lp*IT2
      IG222 = ILTc*IT2
      IG302 = ILTc/XF3/ISPc**3
      IG312 = ILTc**2*IT2/ISPc**2*(3.D0*isdpc**2-ISPc*stp)              &
            & - 3.D0*ILTc*isdpc/ISPc*(ILTc+3.D0*IT2*lp)                 &
            & + 6.D0*lp*(ILTc+lp*IT2)
      IG322 = ILTc*(6.D0*lp*IT2+ILTc*(2.D0-3.D0*isdpc*IT2/ISPc))
      IG332 = IG222*ILTc
      IG402 = -ILTc/XF4/ISPc**4
      IG412 = 5.D0*IT2*ixl1c*ILTc/ISPc*isdpc*                           &
            & (2.D0*ISPc*stp-3.D0*isdpc**2) + ixl1c*ixl2c*ixl3c -       &
            & 24.D0*lp*ixl4c*ixl5c + 12.D0*lp**2*(3.D0*ILTc+2.D0*IT2*lp)
      IG422 = (IT2*ixl1c*ixl2c-12.D0*(ixl4c*ixl3c-lp*ixl5c))*ILTc
      IG432 = (-6.D0*IT2*ixl4c+3.D0*ixl3c)*ILTc**2
      IG442 = IG332*ILTc
!
      ISPc = ISPc*XF1
      isdpc = isdpc*XF1
      stp = stp*XF1
      ixx1c = ILTc + 5.D0*IT2*lp
      IG502 = 864.D0*ETA**4*ILTc/ISPc**5
      IG512 = 120.D0*LT1**3*(2.D0*ILTc+LT1*IT2)                         &
            & + (10.D0*ILTc**2*(6.D0*ILTc**2*isdpc*stp+                 &
            & IT2*ILTc**2*stp**2.D0+45.D0*ILTc*LT1*isdpc**2+            &
            & 30.D0*IT2*ILTc*LT1*isdpc*stp+90.D0*IT2*LT1**2*isdpc**2))  &
            & /ISPc**2 -                                                &
            & (100.D0*ILTc*LT1*(ILTc**2*stp+6.D0*ILTc*LT1*isdpc+        &
            & 6.D0*LT1**2.D0*isdpc*IT2+2.D0*ILTc*LT1*stp*IT2))          &
            & /ISPc - (105.D0*ILTc**3*isdpc**2*(ILTc*isdpc+ILTc*stp*IT2+&
            & 5.D0*LT1*isdpc*IT2))/ISPc**3 +                            &
            & (105.D0*ILTc**4*isdpc**4*IT2)/ISPc**4
      IG522 = 240.D0*ILTc*LT1**2*(ILTc+LT1*IT2)                         &
            & - (20.D0*ILTc**2*(ILTc**2*stp+15.D0*ILTc*LT1*isdpc+       &
            & 30.D0*LT1**2*isdpc*IT2+5.D0*ILTc*LT1*stp*IT2))            &
            & /ISPc + (30.D0*ILTc**3*isdpc*(3.D0*ILTc*isdpc+            &
            & 2.D0*ILTc*stp*IT2+15.D0*LT1*isdpc*IT2))/ISPc**2 -         &
            & (105.D0*ILTc**4*isdpc**3*IT2)/ISPc**3
      IG532 = 5.D0*ILTc**2*(24.D0*IT2*LT1**2+12.D0*ILTc*LT1)            &
            & + (45.D0*ILTc**4*isdpc**2*IT2-                            &
            & 5.D0*ILTc**2*ISPc*(6.D0*ILTc**2*isdpc+                    &
            & 2.D0*ILTc**2*stp*IT2+30.D0*ILTc*LT1*isdpc*IT2))/ISPc**2
      IG542 = 4.D0*ILTc**3*ixx1c - 10.D0*IT2*ILTc**3*ixl4c
      IG552 = IT2*ILTc**4
!
      ISPc = ISPc/XF1
      isdpc = isdpc/XF1
      stp = stp/XF1
      END
!*==GRPRLD.spg  processed by SPAG 6.72Dc at 03:26 on  5 Apr 2020
!
!     !%--------------------------------------------------------------%!
!
!     Calculates the analytical derivatives of the G(R) components
!
      SUBROUTINE GRPRLD(Rho,R)
99001 FORMAT (2X,2F16.6)
      IMPLICIT NONE
!*--GRPRLD269
!*** Start of declarations inserted by SPAG
      REAL*8 a1 , a10 , a11 , a12 , a13 , a14 , a1conv , a2 , a3 , a4 , &
           & a5 , a6 , a7 , a8 , a9 , C1 , C2 , C3 , da0 , da0dt
      REAL*8 da1 , da10 , da11 , da12 , da13 , da14 , da1de , da2 ,     &
           & da3 , da4 , da5 , da6 , da7 , da8 , da9 , db0de , db1 ,    &
           & db1de , db1dt , db2
      REAL*8 db2de , dc0de , dc11 , dc12 , dc1de , dc1dt , dc2 , dc2de ,&
           & dc2dt , dc3 , dc3de , dc3dt , dd0de , dd1 , dd1de , dd1dt ,&
           & dd2 , dd2de , dd2dt , dd3
      REAL*8 dd3de , dd3dt , dd4 , dd4de , dd4dt , de0de , de1 , de1a , &
           & de1b , de1c , de1d , de1de , de1dt , de2 , de2a , de2b ,   &
           & de2c , de2de , de2dt , de3
      REAL*8 de3a , de3b , de3de , de3dt , de4 , de4de , de4dt , de5 ,  &
           & de5de , de5dt , DEC1 , DEC2 , DEC3 , DF1de , DF2de ,       &
           & DF3de , DF4de , DF5de , dil1de , dis2dt
      REAL*8 dl1de , dlde , ds1 , ds1de , ds1dt , ds2 , ds2de , ds2dt , &
           & ds3de , dsde , dsp , DTE , dter , ETA , ETAm , EX11 ,      &
           & IG101 , IG201 , IG211 , IG221
      REAL*8 IG301 , IG311 , IG321 , IG331 , IG401 , IG411 , IG421 ,    &
           & IG431 , IG441 , IG501 , IG511 , IG521 , IG531 , IG541 ,    &
           & IG551 , ILT , ISDp , ISP , IT1 , IXL1
      REAL*8 IXL2 , IXL3 , IXL4 , IXL5 , LT1 , LT2 , PI , R , Rho ,     &
           & SDP , sinv , XF1 , XF2 , XF3 , XF4 , XF5
!*** End of declarations inserted by SPAG
      COMMON /PIE   / PI
      COMMON /BH    / C1 , C2 , C3
      COMMON /DBH   / DEC1 , DEC2 , DEC3
      COMMON /RGR   / IG102 , IG202 , IG212 , IG222 , IG302 , IG312 ,   &
                    & IG322 , IG332 , IG402 , IG412 , IG422 , IG432 ,   &
                    & IG442 , IG502 , IG512 , IG522 , IG532 , IG542 ,   &
                    & IG552 , IT2 , EX12 , IG101 , IG201 , IG211 ,      &
                    & IG221 , IG301 , IG311 , IG321 , IG331 , IG401 ,   &
                    & IG411 , IG421 , IG431 , IG441 , IG501 , IG511 ,   &
                    & IG521 , IG531 , IG541 , IG551 , IT1 , EX11 , ETA
!     Passed from subroutine GRPRLM to subroutine GRPRLD
      COMMON /DERIV / ISP , XF1 , XF2 , XF3 , XF4 , XF5 , ETAm , LT1 ,  &
                    & LT2 , ILT , ISPc , ILTc , SDP , ISDp , IXL1 ,     &
                    & IXL2 , IXL3 , IXL4 , IXL5
!     Common block for derivatives
      COMMON /DRGR1 / DF1de , DIF1de , DTE , DITe , DF2de , DIF2de ,    &
                    & DF3de , DIF3de , DF4de , DIF4de , DF5de , DIF5de
      COMPLEX*16 EX12 , IG102 , IG202 , IG212 , IG222 , IG302 , IG312 , &
               & IG322 , IG332 , IG402 , IG412 , IG422 , IG432 , IG442 ,&
               & IG502 , IG512 , IG522 , IG532 , IG542 , IG552 , IT2
      COMPLEX*16 ILTc , ISPc , siinv , dis3de
!
!     Derivative Variables
!
!     Complex variables Shells One to Three
      COMPLEX*16 disp , disde , DITe , dis1 , dis1dt , dis1de , dilde , &
               & dia0dt , dia0 , DIF1de , dib0de , dib1 , dib1dt ,      &
               & dib1de , dib2 , dib2de , DIF2de , dis2de , dis2
      COMPLEX*16 dia1de , dic0de , dic1dt , dic1de , dic2 , dic2dt ,    &
               & dic2de , dic3 , dic3dt , dic3de , DIF3de
!         --Fourth Shell
      COMPLEX*16 DIF4de , did0de , did1 , did1dt , did1de , did2 ,      &
               & did2dt , did2de , did3 , did3dt , did3de , did4 ,      &
               & did4dt , did4de
      COMPLEX*16 ai1 , ai2 , ai3 , ai4 , ai5 , dai1 , dai2 , dai3 ,     &
               & dai4 , dai5
!        --Fifth Shell
      COMPLEX*16 die0de , die1 , die1dt , die1de , die2 , die2dt ,      &
               & die2de , die3 , die3dt , die3de , die4 , die4dt ,      &
               & die4de , die5 , die5dt , die5de
      COMPLEX*16 die1a , die1b , die1c , die1d , die2a , die2b , die2c ,&
               & die3a , die3b , DIF5de
!
!     !%------------------------DERIVATIVES---------------------------%!
!      Derivatives of Y(R) for R<=1.D0
      dter = (1.D0-ETA)**5
      DEC1 = (4.D0*(1.D0+2.D0*ETA)*(2.D0+ETA))/dter
      DEC2 = -(3.D0*(2.D0+ETA)*(ETA*ETA+9.D0*ETA+2.D0))/(2.D0*dter)
      DEC3 = (1.D0+2*ETA)*(2.D0*ETA*ETA+9.D0*ETA+1.D0)/(2.D0*dter)
!     !%--------------------------------------------------------------%!
!     Derivatives of the GR(R) shells with respect to ETA for R >=1.D0
!     Note: in subroutine GRPRLM ISP = 3(1-ETA)^2*S1, whereas in
!     subroutine GRPRLD,DSP = S1,The "D" represents S1 being used for
!     the derivs
      dsp = ISP*XF1
      dsde = -2.D0*ETAm*IT1**3 + 6.D0*IT1**2*(1.D0-2.D0*ETA)            &
           & + 36.D0*ETA*IT1 - 48.D0*ETA - 12.D0
      DTE = -(dsde/dsp)
!     Derivatives of S1:
      ds1 = 36.D0*ETA - 12.D0*ETA*IT1 + 3.D0*IT1**2*(2.D0*ETA-2.D0)     &
          & - 12.D0*IT1*(ETA-1.D0)
      ds1dt = 6.D0*IT1*(ETA-1.D0)**2 + 12.D0*ETA*(1.D0-ETA)
      ds1de = ds1 + ds1dt*DTE
!     Derivatives of S2:
      ds2 = 6.D0*(2.D0*ETA-2.D0)*IT1 - 12.D0*(2.D0*ETA-1.D0)
      ds2dt = 6.D0*(1.D0-ETA)**2
      ds2de = ds2 + 6.D0*(ETA-1.D0)**2*DTE
!     Derivative of S3
      ds3de = 12.D0*(ETA-1.D0)
!     Derivative of L with respect to density (ETA) and dL1/dETA
      dlde = IT1/2.D0 + DTE*LT1 + 2.D0
      dl1de = 1.D0/2.D0
!     Derivatives of the function f1 ("0" = zero, not "oh")
      da0dt = ILT/dsp
      da0 = IT1/dsp*dlde - IT1*ILT/(dsp**2)*ds1de
      DF1de = da0 + da0dt*DTE
!     Derivatives of the function f2
      db0de = 12.D0/(dsp**2)*(-ILT-ETA*dlde+2.D0*ETA*ILT*ds1de/dsp)
      db1 = dlde*(1.D0-IT1*ds1dt/dsp)                                   &
          & + IT1*ILT/dsp*(-ds2de+ds1de*ds1dt/dsp) + 2.D0*IT1*dl1de
      db1dt = -ds1dt*ILT/dsp + 2.D0*LT1
      db1de = db1 + db1dt*DTE
      db2 = IT1*dlde
      db2de = db2 + ILT*DTE
      DF2de = db0de*(IG211+IG221*(R)) + IG201*(db1de+db2de*(R-2.D0))
!
!     Derivatives of the function f3
!
      da1de = 2.D0*ILT*dlde/(dsp**2) - (2.D0*(ILT**2)*ds1de)/(dsp**3)
      dc0de = 72.D0*ETA/(dsp**3)                                        &
            & *(2.D0*ILT+ETA*dlde-3.D0*ETA*ILT*ds1de/dsp)
      a1conv = 9.D0*(1.D0-ETA)**4
      sinv = 1.D0/dsp
      dc11 = 6.D0*dl1de*(ILT+IT1*LT1) + 6.D0*LT1*(dlde+IT1*dl1de)       &
           & - sinv*(3.D0*ILT*ds2de*(ILT+3.D0*IT1*LT1))                 &
           & - sinv*(3.D0*ds1dt*dlde*(ILT+3.D0*IT1*LT1))                &
           & - sinv*(3.D0*ILT*ds1dt*(dlde+3.D0*IT1*dl1de))              &
           & - sinv**2*(IT1*ILT**2*(dsp*ds3de-6.D0*ds1dt*ds2de+         &
           & ds2dt*ds1de))                                              &
           & - sinv**2*(2.D0*IT1*ILT*dlde*(dsp*ds2dt-3.D0*ds1dt**2))    &
           & + sinv**2*(3.D0*ILT*ds1dt*ds1de*(ILT+3.D0*IT1*LT1))        &
           & + sinv**3*(2.D0*IT1*ILT**2*ds1de*(dsp*ds2dt-3.D0*ds1dt**2))
      dc1dt = 6.D0*LT1**2 - sinv**2.D0*(ILT**2*(dsp*ds2dt-3.D0*ds1dt**2)&
            & ) - 9.D0*ILT*LT1*ds1dt/dsp
      dc1de = dc11 + dc1dt*DTE
      dc2 = dlde*(6.D0*IT1*LT1-ILT*((3.D0*IT1*ds1dt)/dsp-2.D0))         &
          & - ILT*(ILT*((3.D0*IT1*ds2de)/dsp-(3.D0*IT1*ds1dt*ds1de)     &
          & /dsp**2)-6.D0*IT1*dl1de+dlde*((3.D0*IT1*ds1dt)/dsp-2.D0))
      dc2dt = ILT*(6.D0*LT1-3.D0*ILT*ds1dt/dsp)
      dc2de = dc2 + dc2dt*DTE
      dc3 = 2.D0*IT1*ILT*dlde
      dc3dt = ILT**2
      dc3de = dc3 + dc3dt*DTE
      DF3de = dc0de*(IG311+IG321*(R-3.D0)+IG331*(R-3.D0)**2.D0)         &
            & + IG301*(dc1de+dc2de*(R-3.D0)+dc3de*(R-3.D0)**2)
!
!     Derivative of the function f4
!
      a1 = (ILT/dsp)**2
      a2 = 15.D0*ds1dt**2.D0 - 4.D0*dsp*ds2dt
      a3 = ILT + 4.D0*IT1*LT1
      a4 = ILT*ds1dt/dsp
      a5 = 2.D0*ILT + 3.D0*IT1*LT1
      a6 = ILT + 5.D0*LT1*IT1
      a7 = ILT + 2.D0*LT1*IT1
      a8 = ILT*ds2dt + 15.D0*LT1*ds1dt
      a9 = ILT + LT1*IT1
      a10 = ILT*ds2dt + 6.D0*LT1*ds1dt
      a11 = 2.D0*ILT*ds2dt + 15.D0*LT1*ds1dt
      a12 = 2.D0*ILT + LT1*IT1
      a13 = ILT*ds2dt + 3.D0*LT1*ds1dt
      a14 = ILT*ds2dt + 5.D0*LT1*ds1dt
      da1 = 2.D0*ILT*dlde/dsp**2 - 2.D0*ILT**2.D0*ds1de/dsp**3
      da2 = 30.D0*ds1dt*ds2de - 4.D0*dsp*ds3de - 4.D0*ds2dt*ds1de
      da3 = dlde + 4.D0*IT1*dl1de + 4.D0*LT1*DTE
      da4 = ILT*ds2de/dsp + ds1dt*dlde/dsp - ILT*ds1dt*ds1de/dsp**2
      da5 = 2.D0*dlde + 3.D0*IT1*dl1de + 3.D0*LT1*DTE
      da6 = dlde + 5.D0*IT1*dl1de + 5.D0*LT1*DTE
!      DA6dt(ETA) =5*L1(ETA)
      da7 = dlde + 2.D0*IT1*dl1de + 2.D0*LT1*DTE
!      DA7dt(ETA) =2*L1(ETA)
      da8 = ILT*ds3de + 15.D0*LT1*ds2de + 15.D0*ds1dt*dl1de + ds2dt*dlde
      da9 = dlde + IT1*dl1de + LT1*DTE
!      DA9dt(ETA) =L1(ETA)
      da10 = ILT*ds3de + 6.D0*LT1*ds2de + 6.D0*ds1dt*dl1de + ds2dt*dlde
      da11 = 2.D0*ILT*ds3de + 15.D0*LT1*ds2de + 15.D0*ds1dt*dl1de +     &
           & 2.D0*ds2dt*dlde
      da12 = 2.D0*dlde + IT1*dlde + LT1*DTE
!     DA12dt(ETA) =L1(ETA)
      da13 = ILT*ds3de + 3.D0*LT1*ds2de + 3.D0*ds1dt*dl1de + ds2dt*dlde
      da14 = ILT*ds3de + 5.D0*LT1*ds2de + 5.D0*ds1dt*dl1de + ds2dt*dlde
!
      dd0de = 1152.D0*ETA**3*ILT*ds1de/dsp**5 -                         &
            & 288.D0*ETA**3*dlde/dsp**4 - 864.D0*ETA**2*ILT/dsp**4
      dd1 = 12.D0*LT1**2*(3.D0*dlde+2.D0*IT1*dl1de) + a1*a2*da3 +       &
          & a1*a3*da2 + a2*a3*da1 - 24.D0*a4*a5*dl1de -                 &
          & 24.D0*a4*LT1*da5 - 24.D0*a5*LT1*da4 +                       &
          & 24.D0*LT1*dl1de*(3.D0*ILT+2.D0*IT1*LT1)                     &
          & + (5.D0*IT1*ILT**3*ds1dt*(2.D0*dsp*ds3de-6.D0*ds1dt*ds2de+  &
          & 2.D0*ds2dt*ds1de))                                          &
          & /dsp**3 + (5.D0*IT1*ILT**3*ds2de*(2.D0*dsp*ds2dt-           &
          & 3.D0*ds1dt**2))                                             &
          & /dsp**3 + (15.D0*IT1*ILT**2*ds1dt*dlde*(2.D0*dsp*ds2dt-     &
          & 3.D0*ds1dt**2))                                             &
          & /dsp**3 - (15.D0*IT1*ILT**3*ds1dt*ds1de*(2.D0*dsp*ds2dt-    &
          & 3.D0*ds1dt**2))/dsp**4
!
      dd1dt = 24.D0*LT1**3 +                                            &
            & (5.D0*ILT**3*ds1dt*(2.D0*dsp*ds2dt-3.D0*ds1dt**2))/dsp**3
      dd1de = dd1 + dd1dt*DTE
      dd2 = dlde*(12.D0*a5*LT1-12.D0*a3*a4+IT1*a1*a2)                   &
          & + ILT*(12.D0*a5*dl1de-12.D0*a4*da3-12.D0*a3*da4+            &
          & 12.D0*LT1*da5+IT1*a1*da2+IT1*a2*da1)
!
      dd2dt = a1*a2*ILT
      dd2de = dd2 + dd2dt*DTE
      dd3 = ILT**2*(3.D0*da3-6.D0*IT1*da4)                              &
          & + 2.D0*ILT*dlde*(3.D0*a3-6.D0*IT1*a4)
      dd3dt = -6.D0*a4*ILT**2
      dd3de = dd3 + dd3dt*DTE
      dd4 = 3.D0*IT1*ILT**2.D0*dlde
      dd4dt = ILT**3.D0
      dd4de = dd4 + dd4dt*DTE
!
      DF4de = dd0de*(IG411+IG421*(R-4.D0)+IG431*(R-4.D0)                &
            & **2+IG441*(R-4.D0)**3)                                    &
            & + IG401*(dd1de+dd2de*(R-4.D0)+dd3de*(R-4.D0)              &
            & **2+dd4de*(R-4.D0)**3)
!
!     Derivative of function f5
!
      de0de = (3456.D0*ETA**3*ILT)/dsp**5 + (864.D0*ETA**4*dlde)        &
            & /dsp**5 - (4320.D0*ETA**4*ILT*ds1de)/dsp**6
!
      de1a = 120.D0*LT1**3*(2.D0*dlde+IT1*dl1de)                        &
           & + (10.D0*ILT**2*(45.D0*ILT*ds1dt**2*dl1de+                 &
           & 45.D0*LT1*ds1dt**2*dlde+6.D0*ILT**2*ds1dt*ds3de+           &
           & 6.D0*ILT**2*ds2dt*ds2de+2.D0*IT1*ILT*ds2dt**2*dlde+        &
           & 180.D0*IT1*LT1*ds1dt**2*dl1de+180.D0*IT1*LT1**2*ds1dt*ds2de&
           & +2.D0*IT1*ILT**2*ds2dt*ds3de+90.D0*ILT*LT1*ds1dt*ds2de+    &
           & 12.D0*ILT*ds1dt*ds2dt*dlde+30.D0*IT1*ILT*LT1*ds1dt*ds3de+  &
           & 30.D0*IT1*ILT*LT1*ds2dt*ds2de+30.D0*IT1*ILT*ds1dt*ds2dt*   &
           & dl1de+30.D0*IT1*LT1*ds1dt*ds2dt*dlde))/dsp**2
!
      de1b = 360.D0*LT1**2*dl1de*(2.D0*ILT+IT1*LT1)                     &
           & - (100.D0*ILT*dl1de*(ILT**2.D0*ds2dt+6.D0*ILT*LT1*ds1dt+   &
           & 6.D0*IT1*LT1**2*ds1dt+2.D0*IT1*ILT*LT1*ds2dt))             &
           & /dsp - (100.D0*LT1*dlde*(ILT**2*ds2dt+6.D0*ILT*LT1*ds1dt+  &
           & 6.D0*IT1*LT1**2*ds1dt+2.D0*IT1*ILT*LT1*ds2dt))             &
           & /dsp - (20*ILT**2*ds1de*(45.D0*ILT*LT1*ds1dt**2+           &
           & 6.D0*ILT**2*ds1dt*ds2dt+IT1*ILT**2*ds2dt**2+               &
           & 90.D0*IT1*LT1**2*ds1dt**2+30.D0*IT1*ILT*LT1*ds1dt*ds2dt))  &
           & /dsp**3 -                                                  &
           & (105.D0*ILT**3*ds1dt**2*(ILT*ds2de+ds1dt*dlde+IT1*ILT*ds3de&
           & +5.D0*IT1*LT1*ds2de+5.D0*IT1*ds1dt*dl1de+IT1*ds2dt*dlde))  &
           & /dsp**3
!
      de1c = -(100.D0*ILT*LT1*(ILT**2*ds3de+6.D0*ILT*LT1*ds2de+6.D0*ILT*&
           & ds1dt*dl1de+2.D0*ILT*ds2dt*dlde+6.D0*LT1*ds1dt*dlde+       &
           & 6.D0*IT1*LT1**2*ds2de+2.D0*IT1*ILT*LT1*ds3de+              &
           & 2.D0*IT1*ILT*ds2dt*dl1de+12.D0*IT1*LT1*ds1dt*dl1de+        &
           & 2.D0*IT1*LT1*ds2dt*dlde))                                  &
           & /dsp + (20.D0*ILT*dlde*(45.D0*ILT*LT1*ds1dt**2+            &
           & 6.D0*ILT**2*ds1dt*ds2dt+IT1*ILT**2*ds2dt**2+               &
           & 90.D0*IT1*LT1**2*ds1dt**2+30.D0*IT1*ILT*LT1*ds1dt*ds2dt))  &
           & /dsp**2 + (420.D0*IT1*ILT**3*ds1dt**4*dlde)                &
           & /dsp**4 + (420.D0*IT1*ILT**4*ds1dt**3*ds2de)/dsp**4
!
      de1d = -(420.D0*IT1*ILT**4*ds1dt**4*ds1de)                        &
           & /dsp**5 - (210.D0*ILT**3*ds1dt*ds2de*                      &
           & (ILT*ds1dt+IT1*ILT*ds2dt+5*IT1*LT1*ds1dt))                 &
           & /dsp**3 + (100.D0*ILT*LT1*ds1de*                           &
           & (ILT**2*ds2dt+6.D0*ILT*LT1*ds1dt+6.D0*IT1*LT1**2*ds1dt+    &
           & 2.D0*IT1*ILT*LT1*ds2dt))                                   &
           & /dsp**2 - (315.D0*ILT**2*ds1dt**2*dlde*                    &
           & (ILT*ds1dt+IT1*ILT*ds2dt+5.D0*IT1*LT1*ds1dt))              &
           & /dsp**3 + (315.D0*ILT**3*ds1dt**2*ds1de*                   &
           & (ILT*ds1dt+IT1*ILT*ds2dt+5.D0*IT1*LT1*ds1dt))/dsp**4
      de1 = de1a + de1b + de1c + de1d
!
      de1dt = 120.D0*LT1**4 + (105.D0*ILT**4*ds1dt**4)                  &
            & /dsp**4 + (10.D0*ILT**2*(ILT**2*ds2dt**2+                 &
            & 90.D0*LT1**2*ds1dt**2+30.D0*ILT*LT1*ds1dt*ds2dt))         &
            & /dsp**2 - (105.D0*ILT**3*ds1dt**2*                        &
            & (ILT*ds2dt+5.D0*LT1*ds1dt))                               &
            & /dsp**3 - (100.D0*ILT*LT1*(6.D0*LT1**2*ds1dt+             &
            & 2.D0*ILT*LT1*ds2dt))/dsp
!
      de1de = de1 + de1dt*DTE
      de2a = 240.D0*LT1**2*dlde*(ILT+IT1*LT1)                           &
           & - (20.D0*ILT**2*(ILT**2*ds3de+15.D0*ILT*LT1*ds2de+         &
           & 15.D0*ILT*ds1dt*dl1de+2.D0*ILT*ds2dt*dlde+                 &
           & 15.D0*LT1*ds1dt*dlde+30.D0*IT1*LT1**2*ds2de+               &
           & 5.D0*IT1*ILT*LT1*ds3de+5.D0*IT1*ILT*ds2dt*dl1de+           &
           & 60.D0*IT1*LT1*ds1dt*dl1de+5.D0*IT1*LT1*ds2dt*dlde))        &
           & /dsp + 240.D0*ILT*LT1**2.D0*(dlde+IT1*dl1de)               &
           & + (30.D0*ILT**3*ds2de*(3.D0*ILT*ds1dt+2.D0*IT1*ILT*ds2dt+  &
           & 15.D0*IT1*LT1*ds1dt))                                      &
           & /dsp**2 - (40.D0*ILT*dlde*(ILT**2*ds2dt+                   &
           & 15.D0*ILT*LT1*ds1dt+30.D0*IT1*LT1**2*ds1dt+                &
           & 5.D0*IT1*ILT*LT1*ds2dt))/dsp
!
      de2b = (20.D0*ILT**2*ds1de*(ILT**2*ds2dt+15.D0*ILT*LT1*ds1dt+30.D0&
           & *IT1*LT1**2*ds1dt+5.D0*IT1*ILT*LT1*ds2dt))                 &
           & /dsp**2 + 480.D0*ILT*LT1*dl1de*(ILT+IT1*LT1)               &
           & + (30.D0*ILT**3*ds1dt*(3.D0*ILT*ds2de+3.D0*ds1dt*dlde+     &
           & 2.D0*IT1*ILT*ds3de+15.D0*IT1*LT1*ds2de+                    &
           & 15.D0*IT1*ds1dt*dl1de+2.D0*IT1*ds2dt*dlde))                &
           & /dsp**2 - (420.D0*IT1*ILT**3*ds1dt**3*dlde)                &
           & /dsp**3 - (315.D0*IT1*ILT**4*ds1dt**2*ds2de)               &
           & /dsp**3 + (315.D0*IT1*ILT**4*ds1dt**3*ds1de)               &
           & /dsp**4 + (90.D0*ILT**2*ds1dt*dlde*(3.D0*ILT*ds1dt+        &
           & 2.D0*IT1*ILT*ds2dt+15.D0*IT1*LT1*ds1dt))/dsp**2
      de2c = -(60.D0*ILT**3*ds1dt*ds1de*(3.D0*ILT*ds1dt+2.D0*IT1*ILT*   &
           & ds2dt+15.D0*IT1*LT1*ds1dt))/dsp**3
      de2 = de2a + de2b + de2c
!
      de2dt = 240.D0*ILT*LT1**3 - (105.D0*ILT**4*ds1dt**3.D0)           &
            & /dsp**3 - (20.D0*ILT**2*(30.D0*LT1**2*ds1dt+              &
            & 5.D0*ILT*LT1*ds2dt))                                      &
            & /dsp + (30.D0*ILT**3*ds1dt*(2.D0*ILT*ds2dt+               &
            & 15.D0*LT1*ds1dt))/dsp**2
!
      de2de = de2 + de2dt*DTE
      de3a = 5.D0*ILT**2*(12.D0*ILT*dl1de+12.D0*LT1*dlde+               &
           & 48.D0*IT1*LT1*dl1de)                                       &
           & - (5.D0*ILT**2*ds1de*(6.D0*ILT**2*ds1dt+                   &
           & 2.D0*IT1*ILT**2*ds2dt+30.D0*IT1*ILT*LT1*ds1dt)             &
           & +5.D0*ILT**2*dsp*(6.D0*ILT**2*ds2de+12.D0*ILT*ds1dt*dlde+  &
           & 2.D0*IT1*ILT**2*ds3de+30.D0*IT1*ILT*LT1*ds2de+             &
           & 30.D0*IT1*ILT*ds1dt*dl1de+4.D0*IT1*ILT*ds2dt*dlde+         &
           & 30.D0*IT1*LT1*ds1dt*dlde)-90.D0*IT1*ILT**4*ds1dt*ds2de-    &
           & 180.D0*IT1*ILT**3*ds1dt**2.D0*dlde+                        &
           & 10.D0*ILT*dsp*dlde*(6.D0*ILT**2*ds1dt+                     &
           & 2.D0*IT1*ILT**2*ds2dt+30.D0*IT1*ILT*LT1*ds1dt))/dsp**2
      de3b = (2.D0*ds1de*(5.D0*ILT**2*dsp*(6.D0*ILT**2*ds1dt+2.D0*IT1*  &
           & ILT**2*ds2dt+30.D0*IT1*ILT*LT1*ds1dt)                      &
           & -45.D0*IT1*ILT**4*ds1dt**2))                               &
           & /dsp**3 + 10.D0*ILT*dlde*(12.D0*ILT*LT1+24.D0*IT1*LT1**2)
      de3 = de3a + de3b
!
      de3dt = 120.D0*ILT**2*LT1**2.D0 +                                 &
            & (45.D0*ILT**4*ds1dt**2-5.D0*ILT**2*dsp*(2.D0*ILT**2*ds2dt+&
            & 30.D0*ILT*LT1*ds1dt))/dsp**2
!
      de3de = de3 + de3dt*DTE
      de4 = 4.D0*ILT**3*(dlde+5.D0*IT1*dl1de)                           &
          & + 12.D0*ILT**2*dlde*(ILT+5.D0*IT1*LT1)                      &
          & - (10.D0*IT1*ILT**4*ds2de)                                  &
          & /dsp - (40.D0*IT1*ILT**3*ds1dt*dlde)                        &
          & /dsp + (10.D0*IT1*ILT**4*ds1dt*ds1de)/dsp**2
!
      de4dt = 20.D0*ILT**3*LT1 - (10.D0*ILT**4*ds1dt)/dsp
      de4de = de4 + de4dt*DTE
      de5 = 4.D0*IT1*ILT**3*dlde
      de5dt = ILT**4
      de5de = de5 + de5dt*DTE
!
      DF5de = de0de*(IG511+IG521*(R-5.D0)+IG531*(R-5.D0)                &
            & **2.D0+IG541*(R-5.D0)**3+IG551*(R-5.D0)**4)               &
            & + IG501*(de1de+de2de*(R-5.D0)+de3de*(R-5.D0)              &
            & **2+de4de*(R-5.D0)**3+de5de*(R-5.D0)**4)
!
!     !%-------------------------SECOND ROOT--------------------------%!
!     Derivatives of the GR(R) shells with respect to ETA for R >=1.D0
!     Note: in subroutine GRPRLM ISP = 3(1-ETA)^2*S1, whereas in
!     subroutine GRPRLD,DSP = S1,The "D" represents S1 being used for
!     the derivs
!     The following derivatives are simiar to the above derivatives but
!     in this case the second root of ti (IT2) is used throughout.
!     Note: in subroutine GRPRLM ISPC = 2(1-ETA)^2*S1, whereas in
!     subroutine GRPRLD,DISP = S1 (evaluated at IT2).
      disp = XF1*ISPc
      disde = -2.D0*ETAm*IT2**3 + 6.D0*IT2**2*(1.D0-2.D0*ETA)           &
            & + 36.D0*ETA*IT2 - 48.D0*ETA - 12.D0
!     The following DITE is used frequently and is the derivative of the
!     complex root ti (IT2) with respect to density (Eta)
      DITe = -(disde/disp)
      dis1 = 36.D0*ETA - 12.D0*ETA*IT2 + 3.D0*IT2**2*(2.D0*ETA-2.D0)    &
           & - 12.D0*IT2*(ETA-1.D0)
      dis1dt = 6.D0*IT2*(ETA-1.D0)**2 + 12.D0*ETA*ETAm
      dis1de = dis1 + dis1dt*DITe
!     Derivatives of S2:
      dis2 = 6.D0*(2.D0*ETA-2.D0)*IT2 - 12.D0*(2.D0*ETA-1.D0)
      dis2dt = 6.D0*(1.D0-ETA)**2
      dis2de = dis2 + 6.D0*(ETA-1.D0)**2*DITe
!     Derivative of S3
      dis3de = 12.D0*(ETA-1.D0)
!     Derivative of L with respect to density (ETA) and dL1/dETA
      dilde = IT2/2.D0 + DITe*LT1 + 2.D0
      dil1de = 1.D0/2.D0
!     Derivatives of the function f1 ("0" = zero, not "oh")
      dia0dt = ILTc/disp
      dia0 = IT2/disp*dilde - IT2*ILTc/(disp**2.D0)*dis1de
      DIF1de = dia0 + dia0dt*DITe
!     Derivatives of the function f2
      dib0de = 12.D0/(disp**2.D0)                                       &
             & *(-ILTc-ETA*dilde+2.D0*ETA*ILTc*dis1de/(disp))
      dib1 = dilde*(1.D0-IT2*dis1dt/disp)                               &
           & + IT2*ILTc/disp*(-dis2de+dis1de*dis1dt/disp)               &
           & + 2.D0*IT2*dl1de
      dib1dt = -dis1dt*ILTc/disp + 2.D0*LT1
      dib1de = dib1 + dib1dt*DITe
      dib2 = IT2*dilde
      dib2de = dib2 + ILTc*DITe
      DIF2de = dib0de*(IG212+IG222*(R)) + IG202*(dib1de+dib2de*(R-2.D0))
!     Derivatives of the function f3
      dia1de = 2.D0*ILTc*dilde/disp**2 - (2.D0*ILTc**2*dis1de)/disp**3
!
      dic0de = 72.D0*ETA/(disp**3)                                      &
             & *(2.D0*ILTc+ETA*dilde-3.D0*ETA*ILTc*dis1de/disp)
      siinv = 1.D0/disp
      dc12 = 6.D0*dl1de*(ILTc+IT2*LT1) + 6.D0*LT1*(dilde+IT2*dl1de)     &
           & - siinv*(3.D0*ILTc*dis2de*(ILTc+3.D0*IT2*LT1))             &
           & - siinv*(3.D0*dis1dt*dilde*(ILTc+3.D0*IT2*LT1))            &
           & - siinv*(3.D0*ILTc*dis1dt*(dilde+3.D0*IT2*dl1de))          &
           & - siinv**2.D0*                                             &
           & (IT2*ILTc**2*(disp*dis3de-6.D0*dis1dt*dis2de+dis2dt*dis1de)&
           & )                                                          &
           & - siinv**2*(2.D0*IT2*ILTc*dilde*(disp*dis2dt-3.D0*dis1dt**2&
           & )) + siinv**2*(3.D0*ILTc*dis1dt*dis1de*(ILTc+3.D0*IT2*LT1))&
           & + siinv**3*(2.D0*IT2*ILTc**2.D0*dis1de*(disp*dis2dt-3.D0*  &
           & dis1dt**2))
!
      dic1dt = 6.D0*LT1**2.D0 -                                         &
             & siinv**2*(ILTc**2*(disp*dis2dt-3.D0*dis1dt**2))          &
             & - 9.D0*ILTc*LT1*dis1dt/disp
      dic1de = dc12 + dic1dt*DITe
      dic2 = dilde*(6.D0*IT2*LT1-ILTc*((3.D0*IT2*dis1dt)/disp-2.D0))    &
           & - ILTc*                                                    &
           & (ILTc*((3.D0*IT2*dis2de)/disp-(3.D0*IT2*dis1dt*dis1de)     &
           & /disp**2)-6.D0*IT2*dil1de+                                 &
           & dilde*((3.D0*IT2*dis1dt)/disp-2.D0))
      dic2dt = ILTc*(6.D0*LT1-3.D0*ILTc*dis1dt/disp)
      dic2de = dic2 + dic2dt*DITe
      dic3 = 2.D0*IT2*ILTc*dilde
      dic3dt = ILTc**2
      dic3de = dic3 + dic3dt*DITe
      DIF3de = dic0de*(IG312+IG322*(R-3.D0)+IG332*(R-3.D0)**2)          &
             & + IG302*(dic1de+dic2de*(R-3.D0)+dic3de*(R-3.D0)**2)
!
!    Derivatives of function f4
!
      ai1 = (ILTc/disp)**2
      ai2 = 15.D0*dis1dt**2 - 4.D0*disp*dis2dt
      ai3 = ILTc + 4.D0*IT2*LT1
      ai4 = ILTc*dis1dt/disp
      ai5 = 2.D0*ILTc + 3.D0*IT2*LT1
      dai1 = 2.D0*ILTc*dilde/disp**2 - 2.D0*ILTc**2*dis1de/disp**3
      dai2 = 30.D0*dis1dt*dis2de - 4.D0*disp*dis3de - 4.D0*dis2dt*dis1de
      dai3 = dilde + 4.D0*IT2*dil1de + 4.D0*LT1*DITe
      dai4 = ILTc*dis2de/disp + dis1dt*dilde/disp -                     &
           & ILTc*dis1dt*dis1de/disp**2
      dai5 = 2.D0*dilde + 3.D0*IT2*dil1de + 3.D0*LT1*DITe
!
      did0de = 1152.D0*ETA**3*ILTc*dis1de/disp**5 -                     &
             & 288.D0*ETA**3*dilde/disp**4 - 864.D0*ETA**2*ILTc/disp**4
      did1 = 12.D0*LT1**2*(3.D0*dilde+2.D0*IT2*dil1de) + ai1*ai2*dai3 + &
           & ai1*ai3*dai2 + ai2*ai3*dai1 - 24.D0*ai4*ai5*dil1de -       &
           & 24.D0*ai4*LT1*dai5 - 24.D0*ai5*LT1*dai4 +                  &
           & 24.D0*LT1*dil1de*(3.D0*ILTc+2.D0*IT2*LT1)                  &
           & + (5.D0*IT2*ILTc**3*dis1dt*(2.D0*disp*dis3de-              &
           & 6.D0*dis1dt*dis2de+2.D0*dis2dt*dis1de))/disp**3 +          &
           & (5.D0*IT2*ILTc**3*dis2de*(2.D0*disp*dis2dt-3.D0*dis1dt**2))&
           & /disp**3.D0 +                                              &
           & (15.D0*IT2*ILTc**2*dis1dt*dilde*(2.D0*disp*dis2dt-         &
           & 3.D0*dis1dt**2))/disp**3 -                                 &
           & (15.D0*IT2*ILTc**3*dis1dt*dis1de*(2.D0*disp*dis2dt-        &
           & 3.D0*dis1dt**2))/disp**4
!
      did1dt = 24.D0*LT1**3 +                                           &
             & (5.D0*ILTc**3*dis1dt*(2.D0*disp*dis2dt-3.D0*dis1dt**2))  &
             & /disp**3
      did1de = did1 + did1dt*DITe
!
      did2 = dilde*(12.D0*ai5*LT1-12.D0*ai3*ai4+IT2*ai1*ai2)            &
           & + ILTc*(12.D0*ai5*dil1de-12.D0*ai4*dai3-12.D0*ai3*dai4+    &
           & 12.D0*LT1*dai5+IT2*ai1*dai2+IT2*ai2*dai1)
!
      did2dt = ai1*ai2*ILTc
      did2de = did2 + did2dt*DITe
      did3 = ILTc**2*(3.D0*dai3-6.D0*IT2*dai4)                          &
           & + 2.D0*ILTc*dilde*(3.D0*ai3-6.D0*IT2*ai4)
      did3dt = -6.D0*ai4*ILTc**2
      did3de = did3 + did3dt*DITe
      did4 = 3.D0*IT2*ILTc**2*dilde
      did4dt = ILTc**3
      did4de = did4 + did4dt*DITe
!
      DIF4de = did0de*(IG412+IG422*(R-4.D0)+IG432*(R-4.D0)              &
             & **2+IG442*(R-4.D0)**3)                                   &
             & + IG402*(did1de+did2de*(R-4.D0)+did3de*(R-4.D0)          &
             & **2+did4de*(R-4.D0)**3.D0)
!
!     Derivative of function f5
!
      die0de = (3456.D0*ETA**3*ILTc)/disp**5.D0 + (864.D0*ETA**4*dilde) &
             & /disp**5 - (4320.D0*ETA**4*ILTc*dis1de)/disp**6
!
      die1a = 120.D0*LT1**3*(2.D0*dilde+IT2*dil1de)                     &
            & + (10.D0*ILTc**2*(45.D0*ILTc*dis1dt**2*dil1de+            &
            & 45.D0*LT1*dis1dt**2*dilde+6.D0*ILTc**2*dis1dt*dis3de+     &
            & 6.D0*ILTc**2.D0*dis2dt*dis2de+                            &
            & 2.D0*IT2*ILTc*dis2dt**2*dilde+                            &
            & 180.D0*IT2*LT1*dis1dt**2*dil1de+                          &
            & 180.D0*IT2*LT1**2*dis1dt*dis2de+                          &
            & 2.D0*IT2*ILTc**2*dis2dt*dis3de+                           &
            & 90.D0*ILTc*LT1*dis1dt*dis2de+                             &
            & 12.D0*ILTc*dis1dt*dis2dt*dilde+30.D0*IT2*ILTc*LT1*dis1dt* &
            & dis3de+30.D0*IT2*ILTc*LT1*dis2dt*dis2de+                  &
            & 30.D0*IT2*ILTc*dis1dt*dis2dt*dil1de+                      &
            & 30.D0*IT2*LT1*dis1dt*dis2dt*dilde))/disp**2
!
      die1b = 360.D0*LT1**2*dil1de*(2.D0*ILTc+IT2*LT1)                  &
            & - (100.D0*ILTc*dil1de*(ILTc**2*dis2dt+                    &
            & 6.D0*ILTc*LT1*dis1dt+6.D0*IT2*LT1**2*dis1dt+              &
            & 2.D0*IT2*ILTc*LT1*dis2dt))                                &
            & /disp - (100.D0*LT1*dilde*(ILTc**2*dis2dt+                &
            & 6.D0*ILTc*LT1*dis1dt+6.D0*IT2*LT1**2*dis1dt+              &
            & 2.D0*IT2*ILTc*LT1*dis2dt))                                &
            & /disp - (20*ILTc**2*dis1de*(45.D0*ILTc*LT1*dis1dt**2+     &
            & 6.D0*ILTc**2*dis1dt*dis2dt+IT2*ILTc**2*dis2dt**2+         &
            & 90.D0*IT2*LT1**2*dis1dt**2.D0+30.D0*IT2*ILTc*LT1*dis1dt*  &
            & dis2dt))/disp**3 -                                        &
            & (105.D0*ILTc**3*dis1dt**2*(ILTc*dis2de+dis1dt*dilde+      &
            & IT2*ILTc*dis3de+5.D0*IT2*LT1*dis2de+                      &
            & 5.D0*IT2*dis1dt*dil1de+IT2*dis2dt*dilde))/disp**3
!
      die1c = -                                                         &
            & (100.D0*ILTc*LT1*(ILTc**2*dis3de+6.D0*ILTc*LT1*dis2de+6.D0&
            & *ILTc*dis1dt*dil1de+2.D0*ILTc*dis2dt*dilde+               &
            & 6.D0*LT1*dis1dt*dilde+6.D0*IT2*LT1**2*dis2de+             &
            & 2.D0*IT2*ILTc*LT1*dis3de+2.D0*IT2*ILTc*dis2dt*dil1de+     &
            & 12.D0*IT2*LT1*dis1dt*dil1de+2.D0*IT2*LT1*dis2dt*dilde))   &
            & /disp +                                                   &
            & (20.D0*ILTc*dilde*(45.D0*ILTc*LT1*dis1dt**2+6.D0*ILTc**2* &
            & dis1dt*dis2dt+IT2*ILTc**2*dis2dt**2+                      &
            & 90.D0*IT2*LT1**2*dis1dt**2+30.D0*IT2*ILTc*LT1*dis1dt*     &
            & dis2dt))/disp**2 + (420.D0*IT2*ILTc**3*dis1dt**4*dilde)   &
            & /disp**4 + (420.D0*IT2*ILTc**4*dis1dt**3*dis2de)/disp**4
!
      die1d = -(420.D0*IT2*ILTc**4*dis1dt**4.D0*dis1de)/disp**5 -       &
            & (210.D0*ILTc**3*dis1dt*dis2de*                            &
            & (ILTc*dis1dt+IT2*ILTc*dis2dt+5.D0*IT2*LT1*dis1dt))        &
            & /disp**3 +                                                &
            & (100.D0*ILTc*LT1*dis1de*(ILTc**2*dis2dt+6.D0*ILTc*LT1*    &
            & dis1dt+6.D0*IT2*LT1**2*dis1dt+2.D0*IT2*ILTc*LT1*dis2dt))  &
            & /disp**2 -                                                &
            & (315.D0*ILTc**2*dis1dt**2*dilde*(ILTc*dis1dt+IT2*ILTc*    &
            & dis2dt+5.D0*IT2*LT1*dis1dt))/disp**3 +                    &
            & (315.D0*ILTc**3*dis1dt**2*dis1de*                         &
            & (ILTc*dis1dt+IT2*ILTc*dis2dt+5.D0*IT2*LT1*dis1dt))/disp**4
      die1 = die1a + die1b + die1c + die1d
!
      die1dt = 120.D0*LT1**4 + (105.D0*ILTc**4*dis1dt**4)/disp**4 +     &
             & (10.D0*ILTc**2*(ILTc**2*dis2dt**2+90.D0*LT1**2*dis1dt**2+&
             & 30.D0*ILTc*LT1*dis1dt*dis2dt))/disp**2 -                 &
             & (105.D0*ILTc**3*dis1dt**2.D0*                            &
             & (ILTc*dis2dt+5.D0*LT1*dis1dt))/disp**3 -                 &
             & (100.D0*ILTc*LT1*(6.D0*LT1**2*dis1dt+                    &
             & 2.D0*ILTc*LT1*dis2dt))/disp
      die1de = die1 + die1dt*DITe
!
      die2a = 240.D0*LT1**2*dilde*(ILTc+IT2*LT1)                        &
            & - (20.D0*ILTc**2*(ILTc**2*dis3de+15.D0*ILTc*LT1*dis2de+   &
            & 15.D0*ILTc*dis1dt*dil1de+2.D0*ILTc*dis2dt*dilde+          &
            & 15.D0*LT1*dis1dt*dilde+30.D0*IT2*LT1**2*dis2de+           &
            & 5.D0*IT2*ILTc*LT1*dis3de+5.D0*IT2*ILTc*dis2dt*dil1de+     &
            & 60.D0*IT2*LT1*dis1dt*dil1de+5.D0*IT2*LT1*dis2dt*dilde))   &
            & /disp + 240.D0*ILTc*LT1**2*(dilde+IT2*dil1de)             &
            & + (30.D0*ILTc**3*dis2de*(3.D0*ILTc*dis1dt+                &
            & 2.D0*IT2*ILTc*dis2dt+15.D0*IT2*LT1*dis1dt))/disp**2
!
      die2b = -(40.D0*ILTc*dilde*(ILTc**2*dis2dt+15.D0*ILTc*LT1*dis1dt+ &
            & 30.D0*IT2*LT1**2*dis1dt+5.D0*IT2*ILTc*LT1*dis2dt))        &
            & /disp + (20.D0*ILTc**2*dis1de*(ILTc**2*dis2dt+            &
            & 15.D0*ILTc*LT1*dis1dt+30.D0*IT2*LT1**2*dis1dt+            &
            & 5.D0*IT2*ILTc*LT1*dis2dt))/disp**2.D0 +                   &
            & 480.D0*ILTc*LT1*dil1de*(ILTc+IT2*LT1)                     &
            & + (30.D0*ILTc**3*dis1dt*(3.D0*ILTc*dis2de+                &
            & 3.D0*dis1dt*dilde+2.D0*IT2*ILTc*dis3de+                   &
            & 15.D0*IT2*LT1*dis2de+15.D0*IT2*dis1dt*dil1de+             &
            & 2.D0*IT2*dis2dt*dilde))/disp**2 -                         &
            & (420.D0*IT2*ILTc**3*dis1dt**3*dilde)/disp**3 -            &
            & (315.D0*IT2*ILTc**4*dis1dt**2*dis2de)/disp**3
!
      die2c = (315.D0*IT2*ILTc**4*dis1dt**3*dis1de)/disp**4 +           &
            & (90.D0*ILTc**2*dis1dt*dilde*(3.D0*ILTc*dis1dt+            &
            & 2.D0*IT2*ILTc*dis2dt+15.D0*IT2*LT1*dis1dt))/disp**2 -     &
            & (60.D0*ILTc**3*dis1dt*dis1de*(3.D0*ILTc*dis1dt+           &
            & 2.D0*IT2*ILTc*dis2dt+15.D0*IT2*LT1*dis1dt))/disp**3
      die2 = die2a + die2b + die2c
!
      die2dt = 240.D0*ILTc*LT1**3 - (105.D0*ILTc**4*dis1dt**3)/disp**3 -&
             & (20.D0*ILTc**2*(30.D0*LT1**2*dis1dt+5.D0*ILTc*LT1*dis2dt)&
             & )                                                        &
             & /disp + (30.D0*ILTc**3*dis1dt*(2.D0*ILTc*dis2dt+15.D0*LT1&
             & *dis1dt))/disp**2
      die2de = die2 + die2dt*DITe
!
      die3a = 5.D0*ILTc**2*(12.D0*ILTc*dil1de+12.D0*LT1*dilde+          &
            & 48.D0*IT2*LT1*dil1de)                                     &
            & - (5.D0*ILTc**2*dis1de*(6.D0*ILTc**2*dis1dt+              &
            & 2.D0*IT2*ILTc**2*dis2dt+30.D0*IT2*ILTc*LT1*dis1dt)        &
            & +5.D0*ILTc**2*disp*(6.D0*ILTc**2*dis2de+                  &
            & 12.D0*ILTc*dis1dt*dilde+2.D0*IT2*ILTc**2*dis3de+          &
            & 30.D0*IT2*ILTc*LT1*dis2de+30.D0*IT2*ILTc*dis1dt*dil1de+   &
            & 4.D0*IT2*ILTc*dis2dt*dilde+30.D0*IT2*LT1*dis1dt*dilde)    &
            & -90.D0*IT2*ILTc**4*dis1dt*dis2de-                         &
            & 180.D0*IT2*ILTc**3.D0*dis1dt**2.D0*dilde+                 &
            & 10.D0*ILTc*disp*dilde*(6.D0*ILTc**2*dis1dt+               &
            & 2.D0*IT2*ILTc**2*dis2dt+30.D0*IT2*ILTc*LT1*dis1dt))       &
            & /disp**2
!
      die3b = (2.D0*dis1de*(5.D0*ILTc**2*disp*(6.D0*ILTc**2*dis1dt+2.D0*&
            & IT2*ILTc**2*dis2dt+30.D0*IT2*ILTc*LT1*dis1dt)             &
            & -45.D0*IT2*ILTc**4*dis1dt**2))/disp**3 +                  &
            & 10.D0*ILTc*dilde*(12.D0*ILTc*LT1+24.D0*IT2*LT1**2)
      die3 = die3a + die3b
!
      die3dt = 120.D0*ILTc**2*LT1**2.D0 +                               &
             & (45.D0*ILTc**4*dis1dt**2-5.D0*ILTc**2*disp*              &
             & (2.D0*ILTc**2*dis2dt+30.D0*ILTc*LT1*dis1dt))/disp**2
      die3de = die3 + die3dt*DITe
!
      die4 = 4.D0*ILTc**3*(dilde+5.D0*IT2*dil1de)                       &
           & + 12.D0*ILTc**2*dilde*(ILTc+5.D0*IT2*LT1)                  &
           & - (10.D0*IT2*ILTc**4*dis2de)                               &
           & /disp - (40.D0*IT2*ILTc**3*dis1dt*dilde)                   &
           & /disp + (10.D0*IT2*ILTc**4*dis1dt*dis1de)/disp**2
      die4dt = 20.D0*ILTc**3*LT1 - (10.D0*ILTc**4*dis1dt)/disp
      die4de = die4 + die4dt*DITe
!
      die5 = 4.D0*IT2*ILTc**3*dilde
      die5dt = ILTc**4
      die5de = die5 + die5dt*DITe
!
      DIF5de = die0de*(IG512+IG522*(R-5.D0)+IG532*(R-5.D0)              &
             & **2+IG542*(R-5.D0)**3+IG552*(R-5.D0)**4)                 &
             & + IG502*(die1de+die2de*(R-5.D0)+die3de*(R-5.D0)          &
             & **2+die4de*(R-5.D0)**3+die5de*(R-5.D0)**4)
      END
!*==GR.spg  processed by SPAG 6.72Dc at 03:26 on  5 Apr 2020
!
!     !%--------------------------------------------------------------%!
!
      DOUBLE PRECISION FUNCTION GR(R1)
!
!     ORIGINALLY CALCULATED R TIMES RADIAL DISTRIBUTION FUNCTION SO I
!     DIVIDED BY R TO GIVE Y(R)=G(R) FOR R>D AND HAVE ADDED Y(R) FOR
!     R<D.  IF YOU WANT ONLY G(R) THEN IGNORE R<D.
!
      IMPLICIT NONE
!*--GR907
!*** Start of declarations inserted by SPAG
      REAL*8 a , b , c , C1 , C2 , C3 , d , ETA , ex1 , EX11 , hr ,     &
           & IG101 , IG201 , IG211 , IG221 , IG301 , IG311 , IG321 ,    &
           & IG331 , IG401
      REAL*8 IG411 , IG421 , IG431 , IG441 , IG501 , IG511 , IG521 ,    &
           & IG531 , IG541 , IG551 , IT1 , PI , r , R1 , r3 , rho ,     &
           & rrg , rx
!*** End of declarations inserted by SPAG
      COMMON /PIE   / PI
      COMMON /BH    / C1 , C2 , C3
      COMMON /RGR   / IG102 , IG202 , IG212 , IG222 , IG302 , IG312 ,   &
                    & IG322 , IG332 , IG402 , IG412 , IG422 , IG432 ,   &
                    & IG442 , IG502 , IG512 , IG522 , IG532 , IG542 ,   &
                    & IG552 , IT2 , EX12 , IG101 , IG201 , IG211 ,      &
                    & IG221 , IG301 , IG311 , IG321 , IG331 , IG401 ,   &
                    & IG411 , IG421 , IG431 , IG441 , IG501 , IG511 ,   &
                    & IG521 , IG531 , IG541 , IG551 , IT1 , EX11 , ETA
      COMPLEX*16 EX12 , IG102 , IG202 , IG212 , IG222 , IG302 , IG312 , &
               & IG322 , IG332 , IG402 , IG412 , IG422 , IG432 , IG442 ,&
               & IG502 , IG512 , IG522 , IG532 , IG542 , IG552 , IT2
      COMPLEX*16 CDEXP , ex2
!
      r = R1
      IF ( r.GE..99999D0 .AND. r.LE.6.D0 ) THEN
!     The following calculates G(R) analytically for 1<=R<=6.
         rrg = 0.D0
!
!         FIRST SHELL CONTRIBUTION
!
         ex1 = EX11*DEXP(r*IT1)
         ex2 = EX12*CDEXP(r*IT2)
         rrg = rrg + IG101*ex1 + 2.D0*IG102*ex2
         IF ( r.GT.2.D0 ) THEN
!
!         SECOND SHELL CONTRIBUTION
!
            ex1 = ex1*EX11
            ex2 = ex2*EX12
            rrg = rrg + IG201*ex1*(IG211+IG221*r)                       &
                & + 2.D0*IG202*ex2*(IG212+IG222*r)
            IF ( r.GT.3.D0 ) THEN
!
!         THIRD SHELL CONTRIBUTION
!
               rx = r - 3.D0
               ex1 = ex1*EX11
               ex2 = ex2*EX12
               rrg = rrg + IG301*ex1*((IG331*rx+IG321)*rx+IG311)        &
                   & + 2.D0*IG302*ex2*((IG332*rx+IG322)*rx+IG312)
               IF ( r.GT.4.D0 ) THEN
!
!         FOURTH SHELL CONTRIBUTION
!
                  rx = r - 4.D0
                  ex1 = ex1*EX11
                  ex2 = ex2*EX12
                  rrg = rrg + IG401*ex1*(((IG441*rx+IG431)*rx+IG421)    &
                      & *rx+IG411)                                      &
                      & + 2.D0*ex2*IG402*(((IG442*rx+IG432)*rx+IG422)   &
                      & *rx+IG412)
                  IF ( r.GT.5.D0 ) THEN
!
!         FIFTH SHELL CONTRIBUTION
!
                     rx = r - 5.D0
                     ex1 = ex1*EX11
                     ex2 = ex2*EX12
                     rrg = rrg +                                        &
                         & IG501*ex1*((((IG551*rx+IG541)*rx+IG531)*rx+  &
                         & IG521)*rx+IG511)                             &
                         & + 2.D0*ex2*IG502*((((IG552*rx+IG542)         &
                         & *rx+IG532)*rx+IG522)*rx+IG512)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         GR = rrg/r
      ELSEIF ( r.GT.6.D0 ) THEN
!     The following extrapolates G(R) past r=6 empirically.
         rho = ETA*6.D0/PI
         a = 3.64779052 - 4.5326958*rho + 2.8254675*rho**1.5D0 -        &
           & 1.1230973*DLOG(rho)/rho + 0.02371867*DLOG(rho)/rho**2
         b = 7.0399098 - 13.686034*rho - 0.85088666*rho**3 +            &
           & 2.9323443*DEXP(rho) + 6.6465928*DSQRT(rho)*DLOG(rho)
         c = 5.9382949 + 3.7429435*rho + 1.108728*rho**2*DLOG(rho)      &
           & - 0.96469546*DEXP(rho) + 0.19983335*DLOG(rho)
         d = 1.1350252 - 1.9259676*rho - 1.8835841*rho**2*DLOG(rho)     &
           & + 1.0381572*rho**3 + 0.00007112*DLOG(rho)
         hr = a*DEXP(-b*r)*DSIN(c*r+d)/r
         GR = hr + 1.D0
         RETURN
      ELSE
!     The following calculates Y(R) for r<1.
         r3 = r*r*r
         GR = C1 + C2*r + C3*r3
         RETURN
      ENDIF
      END
!*==DRRG.spg  processed by SPAG 6.72Dc at 03:26 on  5 Apr 2020
!
!     !%--------------------------------------------------------------%!
!
      DOUBLE PRECISION FUNCTION DRRG(R1)
99001 FORMAT (2X,2F16.6)
!
!     Calculates the derivative of the Radial Distribution Function
!
      IMPLICIT NONE
!*--DRRG1016
!*** Start of declarations inserted by SPAG
      REAL*8 a , b , c , C1 , C2 , C3 , d , da , db , dc , dd , DEC1 ,  &
           & DEC2 , DEC3 , DF1de , DF2de , DF3de , DF4de , DF5de , DTE
      REAL*8 ETA , ex1 , EX11 , fn11 , fn21 , fn31 , fn41 , fn51 ,      &
           & IG101 , IG201 , IG211 , IG221 , IG301 , IG311 , IG321 ,    &
           & IG331 , IG401 , IG411 , IG421 , IG431
      REAL*8 IG441 , IG501 , IG511 , IG521 , IG531 , IG541 , IG551 ,    &
           & IT1 , PI , r , R1 , r3 , rho
!*** End of declarations inserted by SPAG
      COMMON /PIE   / PI
      COMMON /BH    / C1 , C2 , C3
      COMMON /DBH   / DEC1 , DEC2 , DEC3
      COMMON /RGR   / IG102 , IG202 , IG212 , IG222 , IG302 , IG312 ,   &
                    & IG322 , IG332 , IG402 , IG412 , IG422 , IG432 ,   &
                    & IG442 , IG502 , IG512 , IG522 , IG532 , IG542 ,   &
                    & IG552 , IT2 , EX12 , IG101 , IG201 , IG211 ,      &
                    & IG221 , IG301 , IG311 , IG321 , IG331 , IG401 ,   &
                    & IG411 , IG421 , IG431 , IG441 , IG501 , IG511 ,   &
                    & IG521 , IG531 , IG541 , IG551 , IT1 , EX11 , ETA
!     Common block for derivatives
      COMMON /DRGR1 / DF1de , DIF1de , DTE , DITe , DF2de , DIF2de ,    &
                    & DF3de , DIF3de , DF4de , DIF4de , DF5de , DIF5de
      COMPLEX*16 EX12 , IG102 , IG202 , IG212 , IG222 , IG302 , IG312 , &
               & IG322 , IG332 , IG402 , IG412 , IG422 , IG432 , IG442 ,&
               & IG502 , IG512 , IG522 , IG532 , IG542 , IG552 , IT2
      COMPLEX*16 CDEXP , fn12 , fn22 , fn32 , fn42 , ex2
      COMPLEX*16 DITe , DIF1de , DIF2de , DIF3de , DIF4de , DIF5de ,    &
               & fn52
      r = R1
      IF ( r.GE..99999D0 .AND. r.LE.6.0D0 ) THEN
!     The following calculates the density derivative of G(R) analytically
!     for 1<=R<=6.
         DRRG = 0.D0
!
!         FIRST SHELL CONTRIBUTION
!
         ex1 = EX11*DEXP(r*IT1)
         ex2 = EX12*CDEXP(r*IT2)
         fn11 = IG101
         fn12 = IG102
         DRRG = DRRG + ex1*(DF1de+(r-1.D0)*fn11*DTE)                    &
              & + 2.D0*ex2*(DIF1de+(r-1.D0)*fn12*DITe)
         IF ( r.GT.2.D0 ) THEN
!
!         SECOND SHELL CONTRIBUTION
!
            ex1 = ex1*EX11
            ex2 = ex2*EX12
            fn21 = IG201*(IG211+IG221*(r))
            fn22 = IG202*(IG212+IG222*(r))
            DRRG = DRRG + ex1*(DF2de+(r-2.D0)*fn21*DTE)                 &
                 & + 2.D0*ex2*(DIF2de+(r-2.D0)*fn22*DITe)
            IF ( r.GT.3.D0 ) THEN
!
!         THIRD SHELL CONTRIBUTION
!
               ex1 = ex1*EX11
               ex2 = ex2*EX12
               fn31 = IG301*(IG311+IG321*(r-3.D0)+IG331*(r-3.D0)**2.D0)
               fn32 = IG302*(IG312+IG322*(r-3.D0)+IG332*(r-3.D0)**2.D0)
               DRRG = DRRG + ex1*(DF3de+(r-3.D0)*fn31*DTE)              &
                    & + 2.D0*ex2*(DIF3de+(r-3.D0)*fn32*DITe)
               IF ( r.GT.4.D0 ) THEN
!
!         FOURTH SHELL CONTRIBUTION
!
                  ex1 = ex1*EX11
                  ex2 = ex2*EX12
                  fn41 = IG401*(IG411+IG421*(r-4.D0)+IG431*(r-4.D0)     &
                       & **2.D0+IG441*(r-4.D0)**3.D0)
                  fn42 = IG402*(IG412+IG422*(r-4.D0)+IG432*(r-4.D0)     &
                       & **2.D0+IG442*(r-4.D0)**3.D0)
                  DRRG = DRRG + ex1*(DF4de+(r-4.D0)*fn41*DTE)           &
                       & + 2.D0*ex2*(DIF4de+(r-4.D0)*fn42*DITe)
                  IF ( r.GT.5.D0 ) THEN
!
!         FIFTH SHELL CONTRIBUTION
!
                     ex1 = ex1*EX11
                     ex2 = ex2*EX12
                     fn51 = IG501*(IG511+IG521*(r-5.D0)+IG531*(r-5.D0)  &
                          & **2.D0+IG541*(r-5.D0)**3.D0+IG551*(r-5.D0)  &
                          & **4.D0)
                     fn52 = IG502*(IG512+IG522*(r-5.D0)+IG532*(r-5.D0)  &
                          & **2.D0+IG542*(r-5.D0)**3.D0+IG552*(r-5.D0)  &
                          & **4.D0)
                     DRRG = DRRG + ex1*(DF5de+(r-5.D0)*fn51*DTE)        &
                          & + 2.D0*ex2*(DIF5de+(r-5.D0)*fn52*DITe)
                     IF ( r.LE.6.D0 ) THEN
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         DRRG = PI*DRRG/(6.D0*r)
      ELSEIF ( r.GT.6.0D0 ) THEN
!     The following extrapolates the density derivative of G(R) past r=6
!     empirically.
         rho = ETA*6.D0/PI
         a = 3.64779052 - 4.5326958*rho + 2.8254675*rho**1.5D0 -        &
           & 1.1230973*DLOG(rho)/rho + 0.02371867*DLOG(rho)/rho**2
         b = 7.0399098 - 13.686034*rho - 0.85088666*rho**3 +            &
           & 2.9323443*DEXP(rho) + 6.6465928*DSQRT(rho)*DLOG(rho)
         c = 5.9382949 + 3.7429435*rho + 1.108728*rho**2*DLOG(rho)      &
           & - 0.96469546*DEXP(rho) + 0.19983335*DLOG(rho)
         d = 1.1350252 - 1.9259676*rho - 1.8835841*rho**2*DLOG(rho)     &
           & + 1.0381572*rho**3 + 0.00007112*DLOG(rho)
         da = 2.2578288 - 2.6419506/rho - 0.57462084/rho**2 -           &
            & 0.49861017/rho**3 + 0.03826629/rho**4
         db = -1.1150611 + 1.5019557*rho - 1.5365115*rho**1.5D0 +       &
            & 0.79801674*rho**3 - 1.2700068/rho
         dc = 2.5339359 - 5.0383713*rho + 6.0947248*rho**1.5D0 -        &
            & 1.3679322*rho**2.5D0 + 0.20489549/rho
         dd = -18.494398 + 17.445882*rho - 16.716728*DSQRT(rho)         &
            & *DLOG(rho) + 0.45062181/rho**1.5D0 - 0.12072079/rho**2
         DRRG = DEXP(-b*r)                                              &
              & *(da*DSIN(c*r+d)+DCOS(c*r+d)*a*(dc*r+dd)-r*db*DSIN      &
              & (c*r+d)*a)/r
         RETURN
      ELSE
!     The following calculates the density derivative of Y(R) for r<1.
         r3 = r*r*r
         DRRG = (DEC1+DEC2*r+DEC3*r3)*PI/(6.D0)
         RETURN
      ENDIF
      END
