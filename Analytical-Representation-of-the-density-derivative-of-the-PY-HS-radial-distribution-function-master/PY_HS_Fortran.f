      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/PIE/PI
 1003 FORMAT(F16.2,F16.6)
 1004 FORMAT(F16.2,2F16.6)
 1005 FORMAT("Density (RHO)",F16.1)
C     If you find this useful, please cite this as:
C     Braden D. Kelly and William R. Smith and Douglas Henderson, 
C     Analytical representation of the density derivative of the 
C     Percus–Yevick hard-sphere radial distribution function,
C     Molecular Physics, 114,16-17,2446-2450,2016
C     https://doi.org/10.1080/00268976.2016.1164908

C     THIS PROGRAM CALCULATES BOTH Y(R), WHICH IS EQUAL TO G(R) FOR R>1.
C     AS WELL AS THE DENSITY DERIVATIVE OF Y(R). THE DISTANCES ARE 
C     EXPECTED TO BE IN UNITS OF R/D AND THE DENSITY IS EXPECTED TO BE 
C     IN UNITS OF RHO*(D**3) WHERE D IS THE HARD SPHERE DIAMETER. YOU 
C     PROCEED IN THE FOLLOWING FASHION.          
      ONE=1.D0
      PI=4.D0*DATAN(ONE)                                                                                
      D=1.D0 
C     YOU CAN USE ANY UNITS YOU WANT BUT IF D=1 DISTANCES ARE IN UNITS OF R/D 
      DR=0.1D0
      NR=1
      RHO=1.0D0
      RHOD=RHO*(D**3)
C     GRPRLM GIVES THE PY APPROXIMATION; 
      CALL GRPRLM(RHOD)                                                    
C
      PRINT 1005,RHO                                                                     
C
C     YOU HAVE NOW SET UP THE CONSTANTS WHICH DEPEND UPON RHO BUT ARE 
c     INDEPENDANT OF R.
C                                                                                
C     WHEN YOU COME TO THE PART OF THE PROGRAM WHICH NEEDS Y(R) OR G(R)
C     DO SOMETHING LIKE THE FOLLOWING
      FLAG=1.D0
C     IF FLAG >=1 THEN G(R) AS WELL AS THE DERIVATIVE OF G(R) WILL BE 
C     CALCULATED. OTHERWISE ONLY G(R) WILL BE CALCULATED.
C                                        
      IF(FLAG.GE.1.D0) THEN                                                                          
      DO 50 J=1,101                                                             
      R=(J-1)*DR                                                                      
      RD=R/D                                                                         
      Y=GR(RD)
      CALL GRPRLD(RHOD,RD)                                                                       
      DY=DRRG(RD)
C      IF(RD.LT.1.D0)Y=0.D0
C      IF(RD.LT.1.D0)DY=0.D0
      PRINT 1004,RD,Y,DY   
   50 CONTINUE
      ELSE
      DO 100 J=1,101                                                              
      R=(J-1)*DR                                                                      
      RD=R/D                
      Y=GR(RD)
C      IF(RD.LT.1.D0)Y=0.D0
      PRINT 1003,RD,Y
  100 CONTINUE
      END IF
      END
C
C     !%--------------------------------------------------------------%!      
C                                                       
      SUBROUTINE GRPRLM(RHO)                                                    
      IMPLICIT REAL*8(A-Z)  
      COMMON/PIE/PI                                                    
      COMMON/BH/C1,C2,C3                                                        
      COMMON /RGR/ IG102,IG202,IG212,IG222,IG302,IG312,IG322,IG332        
     1,IG402,IG412,IG422,IG432,IG442,IG502,IG512,IG522,IG532,IG542
     2,IG552,IT2,EX12,IG101,IG201,IG211,IG221,IG301,IG311,IG321,IG331
     3,IG401,IG411,IG421,IG431,IG441,IG501,IG511,IG521,IG531,IG541
     4,IG551,IT1,EX11,ETA                                                                      
C     Common variables to pass to derivative subroutine GRPRLD                                     
      COMMON/DERIV/ ISP,XF1,XF2,XF3,XF4,XF5,ETAM,LT1,LT2,ILT,ISPC,ILTC
     1,SDP,ISDP,IXL1,IXL2,IXL3,IXL4,IXL5      
      COMPLEX*16 EX12,IG102,IG202,IG212,IG222,IG302,IG312,IG322,IG332        
     1,IG402,IG412,IG422,IG432,IG442,IG502,IG512,IG522,IG532,IG542
     2,IG552,IT2                                           
      COMPLEX*16 DCMPLX,CDEXP,IXL1C,IXL2C,IXL3C,IXL4C,IXL5C,ILTC,JX        
     1,ISPC,ISDPC,IXX1C
C           
      ETA=PI*RHO/6.D0                                                          
      ETAM =1.D0-ETA                                                            
      BOT1=ETAM**4                                                              
      TOP=1.D0+2.D0*ETA                                                         
      TOP1=TOP*TOP                                                              
      C1=TOP1/BOT1                                                              
      BOT2=4.D0*BOT1                                                            
      TOPP=2.D0+ETA                                                             
      TOP2=TOPP*TOPP                                                            
      C2=-TOP2/BOT2                                                             
      C2=6.D0*ETA*C2                                                          
      C3=ETA*C1/2.D0                                                            
      ETA12=12.D0*ETA                                                           
      LT1=1.D0+ETA/2.D0                                                         
      LT2=1.D0+2.D0*ETA                                                         
      LP=LT1                                                                    
      XETA=1.D0/(1.D0-ETA)                                                      
      ETA3=(1.D0-ETA)**2                                                        
      SP1=1.D0                                                                  
      SP2=4.D0*ETA*XETA                                                         
      SP3=6.D0*(ETA*XETA)**2                                                    
      SDP1=2.D0                                                                 
      SDP2=SP2                                                                  
      STP=2.D0                                                                  
      FF=(-ETA+3.D0)*ETA+3.D0                                                   
      PAR=DSQRT(1.D0+2.D0*(ETA**2/FF)**2)                                       
      YPLUS=(1.D0+PAR)**(1.D0/3.D0)                                             
      YMINUS=-(-1.D0+PAR)**(1.D0/3.D0)                                          
      JX=CDEXP(DCMPLX(0.D0,2.D0*PI/3.D0))                                      
      PAR=(2.D0*ETA*FF)**(1.D0/3.D0)                                            
      IT1=XETA*(-2.D0*ETA+PAR*(YPLUS+YMINUS))                                 
      IT2=XETA*(-2.D0*ETA+PAR*(YPLUS*JX+YMINUS/JX))                          
      XF1=3.D0*ETA3                                                             
      PAR=ETA3**2/ETA                                                           
      XF2=PAR*3.D0/4.D0                                                         
      PAR=PAR*ETA3/ETA                                                          
      XF3=PAR*3.D0/8.D0                                                         
      PAR=PAR*ETA3/ETA                                                          
      XF4=PAR*9.D0/32.D0
C     !%----------------- Variables for First Root -------------------%!                                                       
      ILT=LT1*IT1+LT2                                                           
      ISDP=SDP1*IT1+SDP2                                                        
      ISP=(SP1*IT1+SP2)*IT1+SP3                                                 
      EX11=DEXP(-IT1)                                                           
      IXL1=(ILT/ISP)**2                                                         
      IXL2=15.D0*ISDP**2-4.D0*ISP*STP                                           
      IXL3=ILT+4.D0*IT1*LP                                                      
      IXL4=ILT*ISDP/ISP                                                         
      IXL5=2.D0*ILT+3.D0*IT1*LP                                                 
      IG101=IT1*ILT/XF1/ISP                                                     
      IG201=-ILT/XF2/ISP**2                                                     
      IG211=ILT*(1.D0-IT1*(2.D0+ISDP/ISP))+2.D0*LP*IT1                          
      IG221=ILT*IT1                                                             
      IG301=ILT/XF3/ISP**3                                                      
      IG311=ILT**2*IT1/ISP**2*(3.D0*ISDP**2-ISP*STP)-3.D0*ILT*ISDP/ISP*(        
     1ILT+3.D0*IT1*LP)+6.D0*LP*(ILT+LP*IT1)                                     
      IG321=ILT*(6.D0*LP*IT1+ILT*(2.D0-3.D0*ISDP*IT1/ISP))                      
      IG331=IG221*ILT                                                           
      IG401=-ILT/XF4/ISP**4                                                     
      IG411=5.D0*IT1*IXL1*ILT/ISP*ISDP*(2.D0*ISP*STP-3.D0*ISDP**2)+IXL1*        
     1IXL2*IXL3-24.D0*LP*IXL4*IXL5+12.D0*LP**2*(3.D0*ILT+2.D0*IT1*LP)           
      IG421=(IT1*IXL1*IXL2-12.D0*(IXL4*IXL3-LP*IXL5))*ILT                       
      IG431=(-6.D0*IT1*IXL4+3.D0*IXL3)*ILT**2                                   
      IG441=IG331*ILT
C
      ISP=ISP*XF1
      ISDP=ISDP*XF1
      STP=STP*XF1
      IXX1=ILT+5.D0*IT1*LP
      IG501=864.D0*ETA**4.D0*ILT/(ISP**5)
C
      IG511=120.D0*LT1**3*(2.D0*ILT+LT1*IT1)+(10.D0*ILT**2*(6.D0*ILT**2
     1*ISDP*STP+ IT1*ILT**2*STP**2+45.D0*ILT*LT1*ISDP**2+30.D0*IT1*ILT
     2*LT1*ISDP*STP+90.D0*IT1*LT1**2.D0*ISDP**2))/ISP**2-(100.D0*ILT*LT1
     3*(ILT**2*STP+ 6.D0*ILT*LT1*ISDP+6.D0*LT1**2*ISDP*IT1+2.D0*ILT*LT1
     4*STP*IT1))/ISP-(105.D0*ILT**3*ISDP**2*(ILT*ISDP+ILT*STP*IT1+5.D0
     5*LT1*ISDP*IT1))/ISP**3+(105.D0*ILT**4*ISDP**4*IT1)/ISP**4 
C
      IG521=240.D0*ILT*LT1**2*(ILT+LT1*IT1)-(20.D0*ILT**2*(ILT**2.D0*STP
     1+15.D0*ILT*LT1*ISDP+30.D0*LT1**2*ISDP*IT1+5.D0*ILT*LT1*STP*IT1))
     2/ISP+(30.D0*ILT**3*ISDP*(3.D0*ILT*ISDP+2.D0*ILT*STP*IT1+15.D0*LT1
     3*ISDP*IT1))/ISP**2-(105.D0*ILT**4*ISDP**3*IT1)/ISP**3
C
      IG531=5.D0*ILT**2*(24.D0*IT1*LT1**2+12.D0*ILT*LT1)+(45.D0*ILT**4
     1*ISDP**2*IT1- 5.D0*ILT**2*ISP*(6.D0*ILT**2*ISDP+2.D0*ILT**2*STP
     2*IT1+30.D0*ILT*LT1*ISDP*IT1))/ISP**2
C
      IG541=4.D0*ILT**3*IXX1-10.D0*IT1*ILT**3*IXL4
      IG551=IT1*ILT**4
C
      ISP=ISP/XF1
      ISDP=ISDP/XF1
      STP=STP/XF1
C     !%--------------- Variables for Second Root --------------------%!                                                     
      ILTC=LT1*IT2+LT2                                                          
      ISDPC=SDP1*IT2+SDP2                                                       
      ISPC=(SP1*IT2+SP2)*IT2+SP3                                                
      EX12 = CDEXP((0.D0,0.D0) - IT2)                                           
      IXL1C=(ILTC/ISPC)**2                                                      
      IXL2C=15.D0*ISDPC**2-4.D0*ISPC*STP                                        
      IXL3C=ILTC+4.D0*IT2*LP                                                    
      IXL4C=ILTC*ISDPC/ISPC                                                     
      IXL5C=2.D0*ILTC+3.D0*IT2*LP                                               
      IG102=IT2*ILTC/XF1/ISPC                                                   
      IG202=-ILTC/XF2/ISPC**2                                                   
      IG212=ILTC*(1.D0-IT2*(2.D0+ISDPC/ISPC))+2.D0*LP*IT2                       
      IG222=ILTC*IT2                                                            
      IG302=ILTC/XF3/ISPC**3                                                    
      IG312=ILTC**2*IT2/ISPC**2*(3.D0*ISDPC**2-ISPC*STP)-3.D0*ILTC*ISDPC        
     1/ISPC*(ILTC+3.D0*IT2*LP)+6.D0*LP*(ILTC+LP*IT2)                            
      IG322=ILTC*(6.D0*LP*IT2+ILTC*(2.D0-3.D0*ISDPC*IT2/ISPC))                  
      IG332=IG222*ILTC                                                          
      IG402=-ILTC/XF4/ISPC**4                                                   
      IG412=5.D0*IT2*IXL1C*ILTC/ISPC*ISDPC*(2.D0*ISPC*STP-3.D0*ISDPC**2)        
     1+IXL1C*IXL2C*IXL3C-24.D0*LP*IXL4C*IXL5C+12.D0*LP**2*(3.D0*ILTC+2.D        
     20*IT2*LP)                                                                 
      IG422=(IT2*IXL1C*IXL2C-12.D0*(IXL4C*IXL3C-LP*IXL5C))*ILTC                 
      IG432=(-6.D0*IT2*IXL4C+3.D0*IXL3C)*ILTC**2                                
      IG442=IG332*ILTC 
C
      ISPC=ISPC*XF1
      ISDPC=ISDPC*XF1
      STP=STP*XF1
      IXX1C=ILTC+5.D0*IT2*LP
      IG502=864.D0*ETA**4*ILTC/ISPC**5
      IG512=120.D0*LT1**3*(2.D0*ILTC+LT1*IT2)+(10.D0*ILTC**2*(6.D0*ILTC
     1**2*ISDPC*STP+IT2*ILTC**2*STP**2.D0+45.D0*ILTC*LT1*ISDPC**2+30.D0
     2*IT2*ILTC*LT1*ISDPC*STP+90.D0*IT2*LT1**2*ISDPC**2))/ISPC**2-(100
     3.D0*ILTC*LT1*(ILTC**2*STP+6.D0*ILTC*LT1*ISDPC+6.D0*LT1**2.D0*ISDPC
     4*IT2+2.D0*ILTC*LT1*STP*IT2))/ISPC-(105.D0*ILTC**3*ISDPC**2*(ILTC
     5*ISDPC+ILTC*STP*IT2+5.D0*LT1*ISDPC*IT2))/ISPC**3+(105.D0*ILTC**4
     6*ISDPC**4*IT2)/ISPC**4 
      IG522=240.D0*ILTC*LT1**2*(ILTC+LT1*IT2)-(20.D0*ILTC**2*(ILTC**2
     1*STP+15.D0*ILTC*LT1*ISDPC+30.D0*LT1**2*ISDPC*IT2+5.D0*ILTC*LT1*STP
     2*IT2))/ISPC+(30.D0*ILTC**3*ISDPC*(3.D0*ILTC*ISDPC+2.D0*ILTC*STP
     3*IT2+15.D0*LT1*ISDPC*IT2))/ISPC**2-(105.D0*ILTC**4*ISDPC**3*IT2)
     4/ISPC**3
      IG532=5.D0*ILTC**2*(24.D0*IT2*LT1**2+12.D0*ILTC*LT1)+(45.D0*ILTC**
     14*ISDPC**2*IT2-5.D0*ILTC**2*ISPC*(6.D0*ILTC**2*ISDPC+2.D0*ILTC**2
     2*STP*IT2+30.D0*ILTC*LT1*ISDPC*IT2))/ISPC**2
      IG542=4.D0*ILTC**3*IXX1C-10.D0*IT2*ILTC**3*IXL4C
      IG552=IT2*ILTC**4
C
      ISPC=ISPC/XF1
      ISDPC=ISDPC/XF1
      STP=STP/XF1
      RETURN
      END
C
C     !%--------------------------------------------------------------%!
C
C     Calculates the analytical derivatives of the G(R) components
C
      SUBROUTINE GRPRLD(RHO,R)
 1003 FORMAT(2X,2F16.6)                                                    
      IMPLICIT REAL*8(A-Z)  
      COMMON/PIE/PI                                                   
      COMMON/BH/C1,C2,C3
      COMMON/DBH/DEC1,DEC2,DEC3                                                        
      COMMON /RGR/ IG102,IG202,IG212,IG222,IG302,IG312,IG322,IG332        
     1,IG402,IG412,IG422,IG432,IG442,IG502,IG512,IG522,IG532,IG542
     2,IG552,IT2,EX12,IG101,IG201,IG211,IG221,IG301,IG311,IG321,IG331
     3,IG401,IG411,IG421,IG431,IG441,IG501,IG511,IG521,IG531,IG541
     4,IG551,IT1,EX11,ETA                                                                    
C     Passed from subroutine GRPRLM to subroutine GRPRLD
      COMMON/DERIV/ ISP,XF1,XF2,XF3,XF4,XF5,ETAM,LT1,LT2,ILT,ISPC,ILTC
     1,SDP,ISDP,IXL1,IXL2,IXL3,IXL4,IXL5
C     Common block for derivatives                                   
      COMMON/DRGR1/ DF1DE,DIF1DE,DTE,DITE,DF2DE,DIF2DE,DF3DE,DIF3DE
     1,DF4DE,DIF4DE,DF5DE,DIF5DE
      COMPLEX*16 EX12,IG102,IG202,IG212,IG222,IG302,IG312,IG322,IG332        
     1,IG402,IG412,IG422,IG432,IG442,IG502,IG512,IG522,IG532,IG542
     2,IG552,IT2                                           
      COMPLEX*16 ILTC,ISPC,SIINV,DIS3DE    
C
C     Derivative Variables
C
C     Complex variables Shells One to Three
      COMPLEX*16 DISP,DISDE,DITE,DIS1,DIS1DT,DIS1DE,DILDE,DIA0DT,DIA0
     1,DIF1DE,DIB0DE,DIB1,DIB1DT,DIB1DE,DIB2,DIB2DE,DIF2DE,DIS2DE,DIS2
      COMPLEX*16 DIA1DE,DIC0DE,DIC1DT,DIC1DE,DIC2,DIC2DT
     1,DIC2DE,DIC3,DIC3DT,DIC3DE,DIF3DE
C         --Fourth Shell
      COMPLEX*16 DIF4DE,DID0DE,DID1,DID1DT,DID1DE,DID2,DID2DT
     1,DID2DE,DID3,DID3DT,DID3DE,DID4,DID4DT,DID4DE
      COMPLEX*16 AI1,AI2,AI3,AI4,AI5,DAI1,DAI2,DAI3,DAI4,DAI5
C        --Fifth Shell
      COMPLEX*16 DIE0DE,DIE1,DIE1DT,DIE1DE,DIE2,DIE2DT,DIE2DE,DIE3
     1,DIE3DT,DIE3DE,DIE4,DIE4DT,DIE4DE,DIE5,DIE5DT,DIE5DE
      COMPLEX*16 DIE1a,DIE1b,DIE1c,DIE1d,DIE2a,DIE2b,DIE2c,DIE3a,DIE3b
     2,DIF5DE
C
C     !%------------------------DERIVATIVES---------------------------%!
C      Derivatives of Y(R) for R<=1.D0
      DTER=(1.D0-ETA)**5
      DEC1=(4.D0*(1.D0+2.D0*ETA)*(2.D0+ETA))/DTER
      DEC2=-(3.D0*(2.D0+ETA)*(ETA*ETA+9.D0*ETA+2.D0))/(2.D0*DTER)
      DEC3=(1.D0+2*ETA)*(2.D0*ETA*ETA+9.D0*ETA+1.D0)/(2.D0*DTER)
C     !%--------------------------------------------------------------%!
C     Derivatives of the GR(R) shells with respect to ETA for R >=1.D0
C     Note: in subroutine GRPRLM ISP = 3(1-ETA)^2*S1, whereas in 
C     subroutine GRPRLD,DSP = S1,The "D" represents S1 being used for 
C     the derivs
      DSP=ISP*XF1
      DSDE=-2.D0*ETAM*IT1**3+6.D0*IT1**2*(1.D0-2.D0*ETA)+36.D0*
     1ETA*IT1-48.D0*ETA-12.D0
      DTE=-(DSDE/DSP)
C     Derivatives of S1:
      DS1=36.D0*ETA-12.D0*ETA*IT1+3.D0*IT1**2*(2.D0*ETA-2.D0)-12.D0*IT1
     1*(ETA-1.D0)
      DS1DT=6.D0*IT1*(ETA-1.D0)**2+12.D0*ETA*(1.D0-ETA)
      DS1DE=DS1+DS1DT*DTE
C     Derivatives of S2:
      DS2=6.D0*(2.D0*ETA-2.D0)*IT1-12.D0*(2.D0*ETA-1.D0)
      DS2DT=6.D0*(1.D0-ETA)**2
      DS2DE=DS2+6.D0*(ETA-1.D0)**2*DTE
C     Derivative of S3
      DS3DE=12.D0*(ETA-1.D0)
C     Derivative of L with respect to density (ETA) and dL1/dETA
      DLDE=IT1/2.D0+DTE*LT1+2.D0
      DL1DE=1.D0/2.D0
C     Derivatives of the function f1 ("0" = zero, not "oh")
      DA0DT=ILT/DSP
      DA0=IT1/DSP*DLDE-IT1*ILT/(DSP**2)*DS1DE
      DF1DE=DA0+DA0DT*DTE
C     Derivatives of the function f2
      DB0DE=12.D0/(DSP**2)*(-ILT-ETA*DLDE+2.D0*ETA*ILT*DS1DE/DSP)
      DB1=DLDE*(1.D0-IT1*DS1DT/DSP)+IT1*ILT/DSP*(-DS2DE+DS1DE*DS1DT/DSP
     1)+2.D0*IT1*DL1DE 
      DB1DT=-DS1DT*ILT/DSP+2.D0*LT1
      DB1DE=DB1+DB1DT*DTE
      DB2=IT1*DLDE
      DB2DE=DB2+ILT*DTE
      DF2DE=DB0DE*(IG211+IG221*(R))+IG201*(DB1DE+DB2DE*(R-2.D0))
C
C     Derivatives of the function f3
C
      DA1DE=2.D0*ILT*DLDE/(DSP**2)-(2.D0*(ILT**2)*DS1DE)/(DSP**3)
      DC0DE=72.D0*ETA/(DSP**3)*(2.D0*ILT+ETA*DLDE-3.D0*ETA*ILT*DS1DE/
     1DSP)
      A1CONV=9.D0*(1.D0-ETA)**4
      SINV=1.D0/DSP   
      DC11=6.D0*DL1DE*(ILT+IT1*LT1)+6.D0*LT1*(DLDE+IT1*DL1DE)-SINV*(3.D0
     1*ILT*DS2DE*(ILT+3.D0*IT1*LT1))-SINV*(3.D0*DS1DT*DLDE*(ILT+3.D0*IT1
     2*LT1))-SINV*(3.D0*ILT*DS1DT*(DLDE+3.D0*IT1*DL1DE))-SINV**2*(IT1*IL
     3T**2*(DSP*DS3DE-6.D0*DS1DT*DS2DE+DS2DT*DS1DE))-SINV**2*(2.D0*IT1*I
     4LT*DLDE*(DSP*DS2DT-3.D0*DS1DT**2))+SINV**2*(3.D0*ILT*DS1DT*DS1DE*(
     5ILT+3.D0*IT1*LT1))+SINV**3*(2.D0*IT1*ILT**2*DS1DE*(DSP*DS2DT-3.D0
     6*DS1DT**2))
      DC1DT=6.D0*LT1**2-SINV**2.D0*(ILT**2*(DSP*DS2DT-3.D0*
     1DS1DT**2))-9.D0*ILT*LT1*DS1DT/DSP
      DC1DE=DC11+DC1DT*DTE
      DC2=DLDE*(6.D0*IT1*LT1-ILT*((3.D0*IT1*DS1DT)/DSP-2.D0))-ILT*(ILT*
     1((3.D0*IT1*DS2DE)/DSP-(3.D0*IT1*DS1DT*DS1DE)/DSP**2)-6.D0*IT1*
     2DL1DE+DLDE*((3.D0*IT1*DS1DT)/DSP-2.D0))
      DC2DT=ILT*(6.D0*LT1-3.D0*ILT*DS1DT/DSP)
      DC2DE=DC2+DC2DT*DTE
      DC3=2.D0*IT1*ILT*DLDE
      DC3DT=ILT**2
      DC3DE=DC3+DC3DT*DTE
      DF3DE=DC0DE*(IG311+IG321*(R-3.D0)+IG331*(R-3.D0)**2.D0)+IG301*(DC
     11DE+DC2DE*(R-3.D0)+DC3DE*(R-3.D0)**2)
C
C     Derivative of the function f4
C
      A1=(ILT/DSP)**2
      A2=15.D0*DS1DT**2.D0-4.D0*DSP*DS2DT
      A3=ILT+4.D0*IT1*LT1
      A4=ILT*DS1DT/DSP
      A5=2.D0*ILT+3.D0*IT1*LT1
      A6=ILT+5.D0*LT1*IT1;
      A7=ILT+2.D0*LT1*IT1;
      A8=ILT*DS2DT+15.D0*LT1*DS1DT;
      A9=ILT+LT1*IT1;
      A10=ILT*DS2DT+6.D0*LT1*DS1DT;
      A11=2.D0*ILT*DS2DT+15.D0*LT1*DS1DT;
      A12=2.D0*ILT+LT1*IT1;
      A13=ILT*DS2DT+3.D0*LT1*DS1DT;
      A14=ILT*DS2DT+5.D0*LT1*DS1DT;
      DA1=2.D0*ILT*DLDE/DSP**2-2.D0*ILT**2.D0*DS1DE/DSP**3
      DA2=30.D0*DS1DT*DS2DE-4.D0*DSP*DS3DE-4.D0*DS2DT*DS1DE
      DA3=DLDE+4.D0*IT1*DL1DE+4.D0*LT1*DTE
      DA4=ILT*DS2DE/DSP+DS1DT*DLDE/DSP-ILT*DS1DT*DS1DE/DSP**2
      DA5=2.D0*DLDE+3.D0*IT1*DL1DE+3.D0*LT1*DTE
      DA6=DLDE+5.D0*IT1*DL1DE +5.D0*LT1*DTE
C      DA6dt(ETA) =5*L1(ETA)
      DA7=DLDE+2.D0*IT1*DL1DE+2.D0*LT1*DTE
C      DA7dt(ETA) =2*L1(ETA)
      DA8=ILT*DS3DE+15.D0*LT1*DS2DE+15.D0*DS1DT*DL1DE+DS2DT*DLDE 
      DA9=DLDE+IT1*DL1DE+LT1*DTE 
C      DA9dt(ETA) =L1(ETA)
      DA10=ILT*DS3DE+6.D0*LT1*DS2DE+6.D0*DS1DT*DL1DE+DS2DT*DLDE 
      DA11=2.D0*ILT*DS3DE+15.D0*LT1*DS2DE+15.D0*DS1DT*DL1DE+2.D0*DS2DT
     1*DLDE 
      DA12=2.D0*DLDE+IT1*DLDE+LT1*DTE 
C     DA12dt(ETA) =L1(ETA)
      DA13=ILT*DS3DE+3.D0*LT1*DS2DE+3.D0*DS1DT*DL1DE+DS2DT*DLDE 
      DA14=ILT*DS3DE+5.D0*LT1*DS2DE+5.D0*DS1DT*DL1DE+DS2DT*DLDE
C
      DD0DE=1152.D0*ETA**3*ILT*DS1DE/DSP**5-288.D0*ETA**3*DLDE/DSP**4
     1-864.D0*ETA**2*ILT/DSP**4
      DD1=12.D0*LT1**2*(3.D0*DLDE+2.D0*IT1*DL1DE)+A1*A2*DA3+A1*A3*DA2
     1+A2*A3*DA1-24.D0*A4*A5*DL1DE-24.D0*A4*LT1*DA5-24.D0*A5*LT1*DA4+24.
     2D0*LT1*DL1DE*(3.D0*ILT+2.D0*IT1*LT1)+(5.D0*IT1*ILT**3*DS1DT*(2.
     3D0*DSP*DS3DE-6.D0*DS1DT*DS2DE+2.D0*DS2DT*DS1DE))/DSP**3+(5.D0*IT1
     4*ILT**3*DS2DE*(2.D0*DSP*DS2DT-3.D0*DS1DT**2))/DSP**3+(15.D0*IT1
     5*ILT**2*DS1DT*DLDE*(2.D0*DSP*DS2DT-3.D0*DS1DT**2))/DSP**3-(15.D0
     6*IT1*ILT**3*DS1DT*DS1DE*(2.D0*DSP*DS2DT-3.D0*DS1DT**2))/DSP**4
C
      DD1DT=24.D0*LT1**3+(5.D0*ILT**3*DS1DT*(2.D0*DSP*DS2DT-3.D0*DS1DT
     1**2))/DSP**3
      DD1DE=DD1+DD1DT*DTE
      DD2=DLDE*(12.D0*A5*LT1-12.D0*A3*A4+IT1*A1*A2)+ILT*(12.D0*A5*DL1DE
     1-12.D0*A4*DA3-12.D0*A3*DA4+12.D0*LT1*DA5+IT1*A1*DA2+IT1*A2*DA1)
C
      DD2DT=A1*A2*ILT
      DD2DE=DD2+DD2DT*DTE
      DD3=ILT**2*(3.D0*DA3-6.D0*IT1*DA4)+2.D0*ILT*DLDE*(3.D0*A3-6.D0
     1*IT1*A4)
      DD3DT=-6.D0*A4*ILT**2
      DD3DE=DD3+DD3DT*DTE
      DD4=3.D0*IT1*ILT**2.D0*DLDE
      DD4DT=ILT**3.D0
      DD4DE=DD4+DD4DT*DTE
C
      DF4DE=DD0DE*(IG411+IG421*(R-4.D0)+IG431*(R-4.D0)**2+IG441*(R-4.D0)
     1**3)+IG401*(DD1DE+DD2DE*(R-4.D0)+DD3DE*(R-4.D0)**2+DD4DE*(R-4.D0)
     2**3)
C
C     Derivative of function f5
C 
      DE0DE=(3456.D0*ETA**3*ILT)/DSP**5+(864.D0*ETA**4*DLDE)/DSP**5-
     1(4320.D0*ETA**4*ILT*DS1DE)/DSP**6
C
      DE1A=120.D0*LT1**3*(2.D0*DLDE+IT1*DL1DE)+(10.D0*ILT**2*(45.D0*ILT
     1*DS1DT**2*DL1DE+45.D0*LT1*DS1DT**2*DLDE+6.D0*ILT**2*DS1DT*DS3DE+
     26.D0*ILT**2*DS2DT*DS2DE+2.D0*IT1*ILT*DS2DT**2*DLDE+180.D0*IT1*LT1
     3*DS1DT**2*DL1DE+180.D0*IT1*LT1**2*DS1DT*DS2DE+2.D0*IT1*ILT**2
     4*DS2DT*DS3DE+90.D0*ILT*LT1*DS1DT*DS2DE+12.D0*ILT*DS1DT*DS2DT*DLDE
     5+30.D0*IT1*ILT*LT1*DS1DT*DS3DE+30.D0*IT1*ILT*LT1*DS2DT*DS2DE+30.D0
     6*IT1*ILT*DS1DT*DS2DT*DL1DE+30.D0*IT1*LT1*DS1DT*DS2DT*DLDE))/DSP**2
C
      DE1B=360.D0*LT1**2*DL1DE*(2.D0*ILT+IT1*LT1)-(100.D0*ILT*DL1DE*(ILT
     1**2.D0*DS2DT+6.D0*ILT*LT1*DS1DT+6.D0*IT1*LT1**2*DS1DT+2.D0*IT1*ILT
     2*LT1*DS2DT))/DSP-(100.D0*LT1*DLDE*(ILT**2*DS2DT+6.D0*ILT*LT1*DS1DT
     3+6.D0*IT1*LT1**2*DS1DT+2.D0*IT1*ILT*LT1*DS2DT))/DSP-(20*ILT**2
     4*DS1DE*(45.D0*ILT*LT1*DS1DT**2+6.D0*ILT**2*DS1DT*DS2DT+IT1*ILT**2
     5*DS2DT**2+90.D0*IT1*LT1**2*DS1DT**2+30.D0*IT1*ILT*LT1*DS1DT*DS2DT)
     6)/DSP**3-(105.D0*ILT**3*DS1DT**2*(ILT*DS2DE+DS1DT*DLDE+IT1*ILT
     7*DS3DE+5.D0*IT1*LT1*DS2DE+5.D0*IT1*DS1DT*DL1DE+IT1*DS2DT*DLDE))
     8/DSP**3
C
      DE1C=-(100.D0*ILT*LT1*(ILT**2*DS3DE+6.D0*ILT*LT1*DS2DE+6.D0*ILT
     1*DS1DT*DL1DE+2.D0*ILT*DS2DT*DLDE+6.D0*LT1*DS1DT*DLDE+6.D0*IT1*LT1
     2**2*DS2DE+2.D0*IT1*ILT*LT1*DS3DE+2.D0*IT1*ILT*DS2DT*DL1DE+12.D0
     3*IT1*LT1*DS1DT*DL1DE+2.D0*IT1*LT1*DS2DT*DLDE))/DSP+(20.D0*ILT*DLDE
     4*(45.D0*ILT*LT1*DS1DT**2+6.D0*ILT**2*DS1DT*DS2DT+IT1*ILT**2*DS2DT
     5**2+90.D0*IT1*LT1**2*DS1DT**2+30.D0*IT1*ILT*LT1*DS1DT*DS2DT))/DSP
     6**2+(420.D0*IT1*ILT**3*DS1DT**4*DLDE)/DSP**4+(420.D0*IT1*ILT**4
     7*DS1DT**3*DS2DE)/DSP**4
C
      DE1D=-(420.D0*IT1*ILT**4*DS1DT**4*DS1DE)/DSP**5-(210.D0*ILT**3*DS1
     1DT*DS2DE*(ILT*DS1DT+IT1*ILT*DS2DT+5*IT1*LT1*DS1DT))/DSP**3+(100.D0
     2*ILT*LT1*DS1DE*(ILT**2*DS2DT+6.D0*ILT*LT1*DS1DT+6.D0*IT1*LT1**2
     3*DS1DT+2.D0*IT1*ILT*LT1*DS2DT))/DSP**2-(315.D0*ILT**2*DS1DT**2*DLD
     4E*(ILT*DS1DT+IT1*ILT*DS2DT+5.D0*IT1*LT1*DS1DT))/DSP**3+(315.D0*ILT
     5**3*DS1DT**2*DS1DE*(ILT*DS1DT+IT1*ILT*DS2DT+5.D0*IT1*LT1*DS1DT))
     6/DSP**4
      DE1=DE1A+DE1B+DE1C+DE1D
C
      DE1DT=120.D0*LT1**4+(105.D0*ILT**4*DS1DT**4)/DSP**4+(10.D0*ILT**2
     1*(ILT**2*DS2DT**2+90.D0*LT1**2*DS1DT**2+30.D0*ILT*LT1*DS1DT*DS2DT)
     2)/DSP**2-(105.D0*ILT**3*DS1DT**2*(ILT*DS2DT+5.D0*LT1*DS1DT))/DSP
     3**3-(100.D0*ILT*LT1*(6.D0*LT1**2*DS1DT+2.D0*ILT*LT1*DS2DT))/DSP
C 
      DE1DE=DE1+DE1DT*DTE
      DE2A=240.D0*LT1**2*DLDE*(ILT+IT1*LT1)-(20.D0*ILT**2*(ILT**2*DS3DE
     1+15.D0*ILT*LT1*DS2DE+15.D0*ILT*DS1DT*DL1DE+2.D0*ILT*DS2DT*DLDE+15.
     2D0*LT1*DS1DT*DLDE+30.D0*IT1*LT1**2*DS2DE+5.D0*IT1*ILT*LT1*DS3DE+5.
     3D0*IT1*ILT*DS2DT*DL1DE+60.D0*IT1*LT1*DS1DT*DL1DE+5.D0*IT1*LT1*DS2D
     4T*DLDE))/DSP+240.D0*ILT*LT1**2.D0*(DLDE+IT1*DL1DE)+(30.D0*ILT**3*D
     5S2DE*(3.D0*ILT*DS1DT+2.D0*IT1*ILT*DS2DT+15.D0*IT1*LT1*DS1DT))/DSP
     6**2-(40.D0*ILT*DLDE*(ILT**2*DS2DT+15.D0*ILT*LT1*DS1DT+30.D0*IT1
     7*LT1**2*DS1DT+5.D0*IT1*ILT*LT1*DS2DT))/DSP
C
      DE2B=(20.D0*ILT**2*DS1DE*(ILT**2*DS2DT+15.D0*ILT*LT1*DS1DT+30.D0
     1*IT1*LT1**2*DS1DT+5.D0*IT1*ILT*LT1*DS2DT))/DSP**2+480.D0*ILT*LT1
     2*DL1DE*(ILT+IT1*LT1)+(30.D0*ILT**3*DS1DT*(3.D0*ILT*DS2DE+3.D0*DS1D
     3T*DLDE+2.D0*IT1*ILT*DS3DE+15.D0*IT1*LT1*DS2DE+15.D0*IT1*DS1DT*DL1D
     4E+2.D0*IT1*DS2DT*DLDE))/DSP**2-(420.D0*IT1*ILT**3*DS1DT**3*DLDE)
     5/DSP**3-(315.D0*IT1*ILT**4*DS1DT**2*DS2DE)/DSP**3+(315.D0*IT1*ILT
     6**4*DS1DT**3*DS1DE)/DSP**4+(90.D0*ILT**2*DS1DT*DLDE*(3.D0*ILT*DS1D
     7T+2.D0*IT1*ILT*DS2DT+15.D0*IT1*LT1*DS1DT))/DSP**2
      DE2C=-(60.D0*ILT**3*DS1DT*DS1DE*(3.D0*ILT*DS1DT+2.D0*IT1*ILT
     1*DS2DT+15.D0*IT1*LT1*DS1DT))/DSP**3
      DE2=DE2A+DE2B+DE2C
C 
      DE2DT=240.D0*ILT*LT1**3-(105.D0*ILT**4*DS1DT**3.D0)/DSP**3-(20.D0
     1*ILT**2*(30.D0*LT1**2*DS1DT+5.D0*ILT*LT1*DS2DT))/DSP+(30.D0*ILT**3
     2*DS1DT*(2.D0*ILT*DS2DT+15.D0*LT1*DS1DT))/DSP**2 
C 
      DE2DE=DE2+DE2DT*DTE
      DE3A=5.D0*ILT**2*(12.D0*ILT*DL1DE+12.D0*LT1*DLDE+48.D0*IT1*LT1*DL1
     1DE)-(5.D0*ILT**2*DS1DE*(6.D0*ILT**2*DS1DT+2.D0*IT1*ILT**2*DS2DT
     2+30.D0*IT1*ILT*LT1*DS1DT)+5.D0*ILT**2*DSP*(6.D0*ILT**2*DS2DE+12.D0
     3*ILT*DS1DT*DLDE+2.D0*IT1*ILT**2*DS3DE+30.D0*IT1*ILT*LT1*DS2DE+30.D
     40*IT1*ILT*DS1DT*DL1DE+4.D0*IT1*ILT*DS2DT*DLDE+30.D0*IT1*LT1*DS1DT
     5*DLDE)-90.D0*IT1*ILT**4*DS1DT*DS2DE-180.D0*IT1*ILT**3*DS1DT**2.D0
     6*DLDE+10.D0*ILT*DSP*DLDE*(6.D0*ILT**2*DS1DT+2.D0*IT1*ILT**2*DS2DT
     7+30.D0*IT1*ILT*LT1*DS1DT))/DSP**2
      DE3B=(2.D0*DS1DE*(5.D0*ILT**2*DSP*(6.D0*ILT**2*DS1DT+2.D0*IT1*ILT
     1**2*DS2DT+30.D0*IT1*ILT*LT1*DS1DT)-45.D0*IT1*ILT**4*DS1DT**2))/DSP
     2**3+10.D0*ILT*DLDE*(12.D0*ILT*LT1+24.D0*IT1*LT1**2)
      DE3=DE3A+DE3B
C 
      DE3DT=120.D0*ILT**2*LT1**2.D0+(45.D0*ILT**4*DS1DT**2-5.D0*ILT**2
     1*DSP*(2.D0*ILT**2*DS2DT+30.D0*ILT*LT1*DS1DT))/DSP**2
C 
      DE3DE=DE3+DE3DT*DTE
      DE4=4.D0*ILT**3*(DLDE+5.D0*IT1*DL1DE)+12.D0*ILT**2*DLDE
     1*(ILT+5.D0*IT1*LT1)-(10.D0*IT1*ILT**4*DS2DE)/DSP-(40.D0*IT1
     2*ILT**3*DS1DT*DLDE)/DSP+(10.D0*IT1*ILT**4*DS1DT*DS1DE)
     3/DSP**2
C 
      DE4DT=20.D0*ILT**3*LT1-(10.D0*ILT**4*DS1DT)/DSP
      DE4DE=DE4+DE4DT*DTE
      DE5=4.D0*IT1*ILT**3*DLDE
      DE5DT=ILT**4
      DE5DE=DE5+DE5DT*DTE
C
      DF5DE=DE0DE*(IG511+IG521*(R-5.D0)+IG531*(R-5.D0)**2.D0+IG541*(R-5.
     1D0)**3+IG551*(R-5.D0)**4)+IG501*(DE1DE+DE2DE*(R-5.D0)+DE3DE*(R-5.D
     20)**2+DE4DE*(R-5.D0)**3+DE5DE*(R-5.D0)**4)
C
C     !%-------------------------SECOND ROOT--------------------------%!
C     Derivatives of the GR(R) shells with respect to ETA for R >=1.D0
C     Note: in subroutine GRPRLM ISP = 3(1-ETA)^2*S1, whereas in 
C     subroutine GRPRLD,DSP = S1,The "D" represents S1 being used for 
C     the derivs
C     The following derivatives are simiar to the above derivatives but 
C     in this case the second root of ti (IT2) is used throughout.
C     Note: in subroutine GRPRLM ISPC = 2(1-ETA)^2*S1, whereas in
C     subroutine GRPRLD,DISP = S1 (evaluated at IT2).
      DISP= XF1*ISPC
      DISDE=-2.D0*ETAM*IT2**3+6.D0*IT2**2*(1.D0-2.D0*ETA)+36.D0*ETA*IT2
     1-48.D0*ETA-12.D0
C     The following DITE is used frequently and is the derivative of the
C     complex root ti (IT2) with respect to density (Eta)
      DITE=-(DISDE/DISP) 
      DIS1=36.D0*ETA-12.D0*ETA*IT2+3.D0*IT2**2*(2.D0*ETA-2.D0)-12.D0
     1*IT2*(ETA-1.D0)
      DIS1DT=6.D0*IT2*(ETA-1.D0)**2+12.D0*ETA*ETAM
      DIS1DE=DIS1+DIS1DT*DITE
C     Derivatives of S2:
      DIS2=6.D0*(2.D0*ETA-2.D0)*IT2-12.D0*(2.D0*ETA-1.D0)
      DIS2DT=6.D0*(1.D0-ETA)**2
      DIS2DE=DIS2+6.D0*(ETA-1.D0)**2*DITE
C     Derivative of S3
      DIS3DE=12.D0*(ETA-1.D0)
C     Derivative of L with respect to density (ETA) and dL1/dETA
      DILDE=IT2/2.D0+DITE*LT1+2.D0
      DIL1DE=1.D0/2.D0
C     Derivatives of the function f1 ("0" = zero, not "oh")
      DIA0DT=ILTC/DISP
      DIA0=IT2/DISP*DILDE-IT2*ILTC/(DISP**2.D0)*DIS1DE
      DIF1DE=DIA0+DIA0DT*DITE
C     Derivatives of the function f2
      DIB0DE=12.D0/(DISP**2.D0)*(-ILTC-ETA*DILDE+2.D0*ETA*ILTC*DIS1DE/
     1(DISP))
      DIB1=DILDE*(1.D0-IT2*DIS1DT/DISP)+IT2*ILTC/DISP*(-DIS2DE+DIS1DE
     1*DIS1DT/DISP)+2.D0*IT2*DL1DE 
      DIB1DT=-DIS1DT*ILTC/DISP+2.D0*LT1
      DIB1DE=DIB1+DIB1DT*DITE
      DIB2=IT2*DILDE
      DIB2DE=DIB2+ILTC*DITE
      DIF2DE=DIB0DE*(IG212+IG222*(R))+IG202*(DIB1DE+DIB2DE*(R-2.D0))
C     Derivatives of the function f3
      DIA1DE=2.D0*ILTC*DILDE/DISP**2-(2.D0*ILTC**2*DIS1DE)/DISP**3
C
      DIC0DE=72.D0*ETA/(DISP**3)*(2.D0*ILTC+ETA*DILDE-3.D0*ETA*ILTC*
     1DIS1DE/DISP)
      SIINV=1.D0/DISP
      DC12=6.D0*DL1DE*(ILTC+IT2*LT1)+6.D0*LT1*(DILDE+IT2*DL1DE)-SIINV
     1*(3.D0*ILTC*DIS2DE*(ILTC+3.D0*IT2*LT1))-SIINV*(3.D0*DIS1DT*DILDE
     2*(ILTC+3.D0*IT2*LT1))-SIINV*(3.D0*ILTC*DIS1DT*(DILDE+3.D0*IT2*
     3DL1DE))-SIINV**2.D0*(IT2*ILTC**2*(DISP*DIS3DE-6.D0*DIS1DT*DIS2DE
     4+DIS2DT*DIS1DE))-SIINV**2*(2.D0*IT2*ILTC*DILDE*(DISP*DIS2DT-3.D0
     5*DIS1DT**2))+SIINV**2*(3.D0*ILTC*DIS1DT*DIS1DE*(ILTC+3.D0*IT2*LT1)
     6)+SIINV**3*(2.D0*IT2*ILTC**2.D0*DIS1DE*(DISP*DIS2DT-3.D0*DIS1DT
     7**2))
C
      DIC1DT=6.D0*LT1**2.D0-SIINV**2*(ILTC**2*(DISP*DIS2DT-3.D0*DIS1DT
     1**2))-9.D0*ILTC*LT1*DIS1DT/DISP
      DIC1DE=DC12+DIC1DT*DITE
      DIC2=DILDE*(6.D0*IT2*LT1-ILTC*((3.D0*IT2*DIS1DT)/DISP-2.D0))-ILTC
     1*(ILTC*((3.D0*IT2*DIS2DE)/DISP-(3.D0*IT2*DIS1DT*DIS1DE)/DISP**2
     2)-6.D0*IT2*DIL1DE+DILDE*((3.D0*IT2*DIS1DT)/DISP-2.D0))
      DIC2DT=ILTC*(6.D0*LT1-3.D0*ILTC*DIS1DT/DISP)
      DIC2DE=DIC2+DIC2DT*DITE
      DIC3=2.D0*IT2*ILTC*DILDE
      DIC3DT=ILTC**2
      DIC3DE=DIC3+DIC3DT*DITE
      DIF3DE=DIC0DE*(IG312+IG322*(R-3.D0)+IG332*(R-3.D0)**2)+IG302*(DIC1
     1DE+DIC2DE*(R-3.D0)+DIC3DE*(R-3.D0)**2)
C
C    Derivatives of function f4
C
      AI1=(ILTC/DISP)**2
      AI2=15.D0*DIS1DT**2-4.D0*DISP*DIS2DT
      AI3=ILTC+4.D0*IT2*LT1
      AI4=ILTC*DIS1DT/DISP
      AI5=2.D0*ILTC+3.D0*IT2*LT1
      DAI1=2.D0*ILTC*DILDE/DISP**2-2.D0*ILTC**2*DIS1DE/DISP**3
      DAI2=30.D0*DIS1DT*DIS2DE-4.D0*DISP*DIS3DE-4.D0*DIS2DT*DIS1DE
      DAI3=DILDE+4.D0*IT2*DIL1DE+4.D0*LT1*DITE
      DAI4=ILTC*DIS2DE/DISP+DIS1DT*DILDE/DISP-ILTC*DIS1DT*DIS1DE/DISP**2
      DAI5=2.D0*DILDE+3.D0*IT2*DIL1DE+3.D0*LT1*DITE
C 
      DID0DE=1152.D0*ETA**3*ILTC*DIS1DE/DISP**5-288.D0*ETA**3*DILDE/DISP
     1**4-864.D0*ETA**2*ILTC/DISP**4
      DID1=12.D0*LT1**2*(3.D0*DILDE+2.D0*IT2*DIL1DE)+AI1*AI2*DAI3+AI1
     1*AI3*DAI2+AI2*AI3*DAI1-24.D0*AI4*AI5*DIL1DE-24.D0*AI4*LT1*DAI5-
     224.D0*AI5*LT1*DAI4+24.D0*LT1*DIL1DE*(3.D0*ILTC+2.D0*IT2*LT1)+(5.D0
     3*IT2*ILTC**3*DIS1DT*(2.D0*DISP*DIS3DE-6.D0*DIS1DT*DIS2DE+2.D0*DIS2
     4DT*DIS1DE))/DISP**3+(5.D0*IT2*ILTC**3*DIS2DE*(2.D0*DISP*DIS2DT
     5-3.D0*DIS1DT**2))/DISP**3.D0+(15.D0*IT2*ILTC**2*DIS1DT*DILDE*(2.D0
     6*DISP*DIS2DT-3.D0*DIS1DT**2))/DISP**3-(15.D0*IT2*ILTC**3*DIS1DT*DI
     7S1DE*(2.D0*DISP*DIS2DT-3.D0*DIS1DT**2))/DISP**4
C
      DID1DT=24.D0*LT1**3+(5.D0*ILTC**3*DIS1DT*(2.D0*DISP*DIS2DT-3.D0*DI
     1S1DT**2))/DISP**3
      DID1DE=DID1+DID1DT*DITE
C
      DID2=DILDE*(12.D0*AI5*LT1-12.D0*AI3*AI4+IT2*AI1*AI2)+ILTC*(12.D0*
     1AI5*DIL1DE-12.D0*AI4*DAI3-12.D0*AI3*DAI4+12.D0*LT1*DAI5+IT2*AI1*
     2DAI2+IT2*AI2*DAI1)
C
      DID2DT=AI1*AI2*ILTC
      DID2DE=DID2+DID2DT*DITE
      DID3=ILTC**2*(3.D0*DAI3-6.D0*IT2*DAI4)+2.D0*ILTC*DILDE*(3.D0*
     1AI3-6.D0*IT2*AI4)
      DID3DT=-6.D0*AI4*ILTC**2
      DID3DE=DID3+DID3DT*DITE
      DID4=3.D0*IT2*ILTC**2*DILDE
      DID4DT=ILTC**3
      DID4DE=DID4+DID4DT*DITE
C
      DIF4DE=DID0DE*(IG412+IG422*(R-4.D0)+IG432*(R-4.D0)**2+IG442*(R
     1-4.D0)**3)+IG402*(DID1DE+DID2DE*(R-4.D0)+DID3DE*(R-4.D0)**2+DID4DE
     2*(R-4.D0)**3.D0)
C
C     Derivative of function f5
C
      DIE0DE=(3456.D0*ETA**3*ILTC)/DISP**5.D0+(864.D0*ETA**4*DILDE)/DISP
     1**5-(4320.D0*ETA**4*ILTC*DIS1DE)/DISP**6
C
      DIE1a=120.D0*LT1**3*(2.D0*DILDE+IT2*DIL1DE)+(10.D0*ILTC**2*(45.D0
     1*ILTC*DIS1DT**2*DIL1DE+45.D0*LT1*DIS1DT**2*DILDE+6.D0*ILTC**2*DIS1
     2DT*DIS3DE+6.D0*ILTC**2.D0*DIS2DT*DIS2DE+2.D0*IT2*ILTC*DIS2DT**2*DI
     3LDE+180.D0*IT2*LT1*DIS1DT**2*DIL1DE+180.D0*IT2*LT1**2*DIS1DT*DIS2D
     4E+2.D0*IT2*ILTC**2*DIS2DT*DIS3DE+90.D0*ILTC*LT1*DIS1DT*DIS2DE
     5+12.D0*ILTC*DIS1DT*DIS2DT*DILDE+30.D0*IT2*ILTC*LT1*DIS1DT*DIS3DE
     6+30.D0*IT2*ILTC*LT1*DIS2DT*DIS2DE+30.D0*IT2*ILTC*DIS1DT*DIS2DT*DIL
     71DE+30.D0*IT2*LT1*DIS1DT*DIS2DT*DILDE))/DISP**2
C
      DIE1b=360.D0*LT1**2*DIL1DE*(2.D0*ILTC+IT2*LT1)-(100.D0*ILTC*DIL1DE
     1*(ILTC**2*DIS2DT+6.D0*ILTC*LT1*DIS1DT+6.D0*IT2*LT1**2*DIS1DT+2.D0
     2*IT2*ILTC*LT1*DIS2DT))/DISP-(100.D0*LT1*DILDE*(ILTC**2*DIS2DT+6.D0
     3*ILTC*LT1*DIS1DT+6.D0*IT2*LT1**2*DIS1DT+2.D0*IT2*ILTC*LT1*DIS2DT))
     4/DISP-(20*ILTC**2*DIS1DE*(45.D0*ILTC*LT1*DIS1DT**2+6.D0*ILTC**2*DI
     5S1DT*DIS2DT+IT2*ILTC**2*DIS2DT**2+90.D0*IT2*LT1**2*DIS1DT**2.D0
     6+30.D0*IT2*ILTC*LT1*DIS1DT*DIS2DT))/DISP**3-(105.D0*ILTC**3*DIS1DT
     7**2*(ILTC*DIS2DE+DIS1DT*DILDE+IT2*ILTC*DIS3DE+5.D0*IT2*LT1*DIS2DE
     8+5.D0*IT2*DIS1DT*DIL1DE+IT2*DIS2DT*DILDE))/DISP**3
C
      DIE1c=-(100.D0*ILTC*LT1*(ILTC**2*DIS3DE+6.D0*ILTC*LT1*DIS2DE+6.D0
     1*ILTC*DIS1DT*DIL1DE+2.D0*ILTC*DIS2DT*DILDE+6.D0*LT1*DIS1DT*DILDE
     2+6.D0*IT2*LT1**2*DIS2DE+2.D0*IT2*ILTC*LT1*DIS3DE+2.D0*IT2*ILTC*DIS
     32DT*DIL1DE+12.D0*IT2*LT1*DIS1DT*DIL1DE+2.D0*IT2*LT1*DIS2DT*DILDE))
     4/DISP+(20.D0*ILTC*DILDE*(45.D0*ILTC*LT1*DIS1DT**2+6.D0*ILTC**2*DIS
     51DT*DIS2DT+IT2*ILTC**2*DIS2DT**2+90.D0*IT2*LT1**2*DIS1DT**2+30.D0
     6*IT2*ILTC*LT1*DIS1DT*DIS2DT))/DISP**2+(420.D0*IT2*ILTC**3*DIS1DT**
     74*DILDE)/DISP**4+(420.D0*IT2*ILTC**4*DIS1DT**3*DIS2DE)/DISP**4
C
      DIE1d=-(420.D0*IT2*ILTC**4*DIS1DT**4.D0*DIS1DE)/DISP**5-(210.D0*IL
     1TC**3*DIS1DT*DIS2DE*(ILTC*DIS1DT+IT2*ILTC*DIS2DT+5.D0*IT2*LT1*DIS1
     2DT))/DISP**3+(100.D0*ILTC*LT1*DIS1DE*(ILTC**2*DIS2DT+6.D0*ILTC*LT1
     3*DIS1DT+6.D0*IT2*LT1**2*DIS1DT+2.D0*IT2*ILTC*LT1*DIS2DT))/DISP**2
     4-(315.D0*ILTC**2*DIS1DT**2*DILDE*(ILTC*DIS1DT+IT2*ILTC*DIS2DT+5.D0
     5*IT2*LT1*DIS1DT))/DISP**3+(315.D0*ILTC**3*DIS1DT**2*DIS1DE*(ILTC
     6*DIS1DT+IT2*ILTC*DIS2DT+5.D0*IT2*LT1*DIS1DT))/DISP**4
      DIE1=DIE1a+DIE1b+DIE1c+DIE1d
C
      DIE1DT=120.D0*LT1**4+(105.D0*ILTC**4*DIS1DT**4)/DISP**4+(10.D0*ILT
     1C**2*(ILTC**2*DIS2DT**2+90.D0*LT1**2*DIS1DT**2+30.D0*ILTC*LT1*DIS1
     2DT*DIS2DT))/DISP**2-(105.D0*ILTC**3*DIS1DT**2.D0*(ILTC*DIS2DT+5.D0
     3*LT1*DIS1DT))/DISP**3-(100.D0*ILTC*LT1*(6.D0*LT1**2*DIS1DT+2.D0*IL
     4TC*LT1*DIS2DT))/DISP 
      DIE1DE=DIE1+DIE1DT*DITE
C
      DIE2a=240.D0*LT1**2*DILDE*(ILTC+IT2*LT1)-(20.D0*ILTC**2*(ILTC**2
     1*DIS3DE+15.D0*ILTC*LT1*DIS2DE+15.D0*ILTC*DIS1DT*DIL1DE+2.D0*ILTC
     2*DIS2DT*DILDE+15.D0*LT1*DIS1DT*DILDE+30.D0*IT2*LT1**2*DIS2DE+5.D0*
     3IT2*ILTC*LT1*DIS3DE+5.D0*IT2*ILTC*DIS2DT*DIL1DE+60.D0*IT2*LT1*DIS1
     4DT*DIL1DE+5.D0*IT2*LT1*DIS2DT*DILDE))/DISP+240.D0*ILTC*LT1**2*(DIL
     5DE+IT2*DIL1DE)+(30.D0*ILTC**3*DIS2DE*(3.D0*ILTC*DIS1DT+2.D0*IT2*IL
     6TC*DIS2DT+15.D0*IT2*LT1*DIS1DT))/DISP**2
C
      DIE2b=-(40.D0*ILTC*DILDE*(ILTC**2*DIS2DT+15.D0*ILTC*LT1*DIS1DT
     1+30.D0*IT2*LT1**2*DIS1DT+5.D0*IT2*ILTC*LT1*DIS2DT))/DISP+(20.D0*IL
     2TC**2*DIS1DE*(ILTC**2*DIS2DT+15.D0*ILTC*LT1*DIS1DT+30.D0*IT2*LT1**
     32*DIS1DT+5.D0*IT2*ILTC*LT1*DIS2DT))/DISP**2.D0+480.D0*ILTC*LT1*DIL
     41DE*(ILTC+IT2*LT1)+(30.D0*ILTC**3*DIS1DT*(3.D0*ILTC*DIS2DE+3.D0*DI
     5S1DT*DILDE+2.D0*IT2*ILTC*DIS3DE+15.D0*IT2*LT1*DIS2DE+15.D0*IT2*DIS
     61DT*DIL1DE+2.D0*IT2*DIS2DT*DILDE))/DISP**2-(420.D0*IT2*ILTC**3*DIS
     71DT**3*DILDE)/DISP**3-(315.D0*IT2*ILTC**4*DIS1DT**2*DIS2DE)/DISP
     8**3
C
      DIE2c=(315.D0*IT2*ILTC**4*DIS1DT**3*DIS1DE)/DISP**4+(90.D0*ILTC**2
     1*DIS1DT*DILDE*(3.D0*ILTC*DIS1DT+2.D0*IT2*ILTC*DIS2DT+15.D0*IT2*LT1
     2*DIS1DT))/DISP**2-(60.D0*ILTC**3*DIS1DT*DIS1DE*(3.D0*ILTC*DIS1DT
     3+2.D0*IT2*ILTC*DIS2DT+15.D0*IT2*LT1*DIS1DT))/DISP**3
      DIE2=DIE2a+DIE2b+DIE2c
C
      DIE2DT=240.D0*ILTC*LT1**3-(105.D0*ILTC**4*DIS1DT**3)/DISP**3-(20.D
     10*ILTC**2*(30.D0*LT1**2*DIS1DT+5.D0*ILTC*LT1*DIS2DT))/DISP+(30.D0
     2*ILTC**3*DIS1DT*(2.D0*ILTC*DIS2DT+15.D0*LT1*DIS1DT))/DISP**2
      DIE2DE=DIE2+DIE2DT*DITE
C
      DIE3a=5.D0*ILTC**2*(12.D0*ILTC*DIL1DE+12.D0*LT1*DILDE+48.D0*IT2*LT
     11*DIL1DE)-(5.D0*ILTC**2*DIS1DE*(6.D0*ILTC**2*DIS1DT+2.D0*IT2*ILTC
     2**2*DIS2DT+30.D0*IT2*ILTC*LT1*DIS1DT)+5.D0*ILTC**2*DISP*(6.D0*ILTC
     3**2*DIS2DE+12.D0*ILTC*DIS1DT*DILDE+2.D0*IT2*ILTC**2*DIS3DE+30.D0
     4*IT2*ILTC*LT1*DIS2DE+30.D0*IT2*ILTC*DIS1DT*DIL1DE+4.D0*IT2*ILTC*DI
     5S2DT*DILDE+30.D0*IT2*LT1*DIS1DT*DILDE)-90.D0*IT2*ILTC**4*DIS1DT*DI
     6S2DE-180.D0*IT2*ILTC**3.D0*DIS1DT**2.D0*DILDE+10.D0*ILTC*DISP*DILD
     7E*(6.D0*ILTC**2*DIS1DT+2.D0*IT2*ILTC**2*DIS2DT+30.D0*IT2*ILTC*LT1
     8*DIS1DT))/DISP**2
C
      DIE3b=(2.D0*DIS1DE*(5.D0*ILTC**2*DISP*(6.D0*ILTC**2*DIS1DT+2.D0*IT
     12*ILTC**2*DIS2DT+30.D0*IT2*ILTC*LT1*DIS1DT)-45.D0*IT2*ILTC**4*DIS1
     2DT**2))/DISP**3+10.D0*ILTC*DILDE*(12.D0*ILTC*LT1+24.D0*IT2*LT1**2)
      DIE3=DIE3a+DIE3b
C
      DIE3DT=120.D0*ILTC**2*LT1**2.D0+(45.D0*ILTC**4*DIS1DT**2-5.D0*ILTC
     1**2*DISP*(2.D0*ILTC**2*DIS2DT+30.D0*ILTC*LT1*DIS1DT))/DISP**2
      DIE3DE=DIE3+DIE3DT*DITE
C
      DIE4=4.D0*ILTC**3*(DILDE+5.D0*IT2*DIL1DE)+12.D0*ILTC**2*DILDE*(ILT
     1C+5.D0*IT2*LT1)-(10.D0*IT2*ILTC**4*DIS2DE)/DISP-(40.D0*IT2*ILTC**3
     2*DIS1DT*DILDE)/DISP+(10.D0*IT2*ILTC**4*DIS1DT*DIS1DE)/DISP**2
      DIE4DT=20.D0*ILTC**3*LT1-(10.D0*ILTC**4*DIS1DT)/DISP
      DIE4DE=DIE4+DIE4DT*DITE
C
      DIE5=4.D0*IT2*ILTC**3*DILDE
      DIE5DT=ILTC**4
      DIE5DE=DIE5+DIE5DT*DITE
C
      DIF5DE=DIE0DE*(IG512+IG522*(R-5.D0)+IG532*(R-5.D0)**2+IG542*(R
     1-5.D0)**3+IG552*(R-5.D0)**4)+IG502*(DIE1DE+DIE2DE*(R-5.D0)+DIE3DE
     2*(R-5.D0)**2+DIE4DE*(R-5.D0)**3+DIE5DE*(R-5.D0)**4)                                 
      RETURN                                                                   
      END 
C
C     !%--------------------------------------------------------------%!
C
      DOUBLE PRECISION FUNCTION GR(R1)                                          
C                                                                               
C     ORIGINALLY CALCULATED R TIMES RADIAL DISTRIBUTION FUNCTION SO I 
C     DIVIDED BY R TO GIVE Y(R)=G(R) FOR R>D AND HAVE ADDED Y(R) FOR 
C     R<D.  IF YOU WANT ONLY G(R) THEN IGNORE R<D.                           
C                                                                               
      IMPLICIT REAL*8(A-Z)
      COMMON/PIE/PI                                                                                                                    
      COMMON/BH/C1,C2,C3                                                                                   
      COMMON /RGR/ IG102,IG202,IG212,IG222,IG302,IG312,IG322,IG332        
     1,IG402,IG412,IG422,IG432,IG442,IG502,IG512,IG522,IG532,IG542
     2,IG552,IT2,EX12,IG101,IG201,IG211,IG221,IG301,IG311,IG321,IG331
     3,IG401,IG411,IG421,IG431,IG441,IG501,IG511,IG521,IG531,IG541
     4,IG551,IT1,EX11,ETA                                     
      COMPLEX*16 EX12,IG102,IG202,IG212,IG222,IG302,IG312,IG322,IG332        
     1,IG402,IG412,IG422,IG432,IG442,IG502,IG512,IG522,IG532,IG542
     2,IG552,IT2                                                                                    
      COMPLEX*16 CDEXP,EX2
C                                                          
      R=R1                                                                      
      IF(R.GE..99999D0.AND.R.LE.6.D0)GO TO 1010                                 
      IF(R.GT.6.D0)GO TO 1011 
C     The following calculates Y(R) for r<1.                                                 
      R3 = R*R*R                                                                
      GR = C1 + C2*R + C3*R3                                                    
      RETURN
C     The following extrapolates G(R) past r=6 empirically.
 1011 RHO = ETA*6.D0/PI
      A=3.64779052-4.5326958*RHO+2.8254675*RHO**1.5D0-1.1230973
     1*DLOG(RHO)/RHO+0.02371867*DLOG(RHO)/RHO**2
      B=7.0399098-13.686034*RHO-0.85088666*RHO**3+2.9323443
     1*DEXP(RHO)+6.6465928*DSQRT(RHO)*DLOG(RHO)
      C=5.9382949+3.7429435*RHO+1.108728*RHO**2*DLOG(RHO)-0.96469546
     1*DEXP(RHO)+0.19983335*DLOG(RHO)
      D=1.1350252-1.9259676*RHO-1.8835841*RHO**2*DLOG(RHO)+1.0381572
     1*RHO**3+0.00007112*DLOG(RHO) 
      HR=A*DEXP(-B*R)*DSIN(C*R+D)/R
      GR=HR+1.D0                                                                                                                                      
      RETURN
C     The following calculates G(R) analytically for 1<=R<=6.                                                                     
 1010 RRG = 0.D0                                                                
C                                                                               
C         FIRST SHELL CONTRIBUTION                                              
C                                                                               
      EX1=EX11*DEXP(R*IT1)                                                      
      EX2=EX12*CDEXP(R*IT2)                                                  
      RRG=RRG+IG101*EX1+2.D0*IG102*EX2                                         
      IF(R.LE.2.D0)GO TO 99                                                     
C                                                                               
C         SECOND SHELL CONTRIBUTION                                             
C                                                                               
      EX1=EX1*EX11                                                              
      EX2=EX2*EX12                                                                   
      RRG=RRG+IG201*EX1*(IG211+IG221*R)+2.D0*IG202*EX2*(IG212+IG222*R)        
      IF(R.LE.3.D0)GO TO 99                                                     
C                                                                               
C         THIRD SHELL CONTRIBUTION                                              
C                                                                               
      RX=R-3.D0                                                                 
      EX1=EX1*EX11                                                              
      EX2=EX2*EX12                                                              
      RRG=RRG+IG301*EX1*((IG331*RX+IG321)*RX+IG311)+2.D0*IG302*EX2*((IG3        
     132*RX+IG322)*RX+IG312)                                                  
      IF(R.LE.4.D0)GO TO 99                                                     
C                                                                               
C         FOURTH SHELL CONTRIBUTION                                             
C                                                                               
      RX=R-4.D0                                                                 
      EX1=EX1*EX11                                                              
      EX2=EX2*EX12                                                                          
      RRG=RRG+IG401*EX1*(((IG441*RX+IG431)*RX+IG421)*RX+IG411)+2.D0*EX2*        
     1IG402*(((IG442*RX+IG432)*RX+IG422)*RX+IG412)       
      IF(R.LE.5.D0)GO TO 99 
C                                                                               
C         FIFTH SHELL CONTRIBUTION                                             
C                                                                               
      RX=R-5.D0                                                                 
      EX1=EX1*EX11                                                              
      EX2=EX2*EX12                                                                     
      RRG=RRG+IG501*EX1*((((IG551*RX+IG541)*RX+IG531)*RX+IG521)*RX+IG511
     1)+2.D0*EX2*IG502*((((IG552*RX+IG542)*RX+IG532)*RX+IG522)*RX+IG512)                                                               
   99 GR = RRG/R                                                               
      RETURN                                                                    
      END
C                                                                     
C     !%--------------------------------------------------------------%!
C
      DOUBLE PRECISION FUNCTION DRRG(R1) 
 1003 FORMAT(2X,2F16.6)                                         
C                                                                               
C     Calculates the derivative of the Radial Distribution Function                          
C                                                                               
      IMPLICIT REAL*8(A-Z)                                                      
      COMMON/PIE/PI                                                               
      COMMON/BH/C1,C2,C3
      COMMON/DBH/DEC1,DEC2,DEC3                                                         
      COMMON /RGR/ IG102,IG202,IG212,IG222,IG302,IG312,IG322,IG332        
     1,IG402,IG412,IG422,IG432,IG442,IG502,IG512,IG522,IG532,IG542
     2,IG552,IT2,EX12,IG101,IG201,IG211,IG221,IG301,IG311,IG321,IG331
     3,IG401,IG411,IG421,IG431,IG441,IG501,IG511,IG521,IG531,IG541
     4,IG551,IT1,EX11,ETA                                                                      
C     Common block for derivatives                                   
      COMMON/DRGR1/ DF1DE,DIF1DE,DTE,DITE,DF2DE,DIF2DE,DF3DE,DIF3DE
     1,DF4DE,DIF4DE,DF5DE,DIF5DE 
      COMPLEX*16 EX12,IG102,IG202,IG212,IG222,IG302,IG312,IG322,IG332
     1,IG402,IG412,IG422,IG432,IG442,IG502,IG512,IG522,IG532,IG542,IG552
     2,IT2
      COMPLEX*16 CDEXP,FN12,FN22,FN32,FN42,EX2
      COMPLEX*16 DITE,DIF1DE,DIF2DE,DIF3DE,DIF4DE,DIF5DE,FN52                            
      R=R1                                                                      
      IF(R.GE..99999D0.AND.R.LE.6.0D0)GO TO 1010                                 
      IF(R.GT.6.0D0)GO TO 1011
C     The following calculates the density derivative of Y(R) for r<1.                                                   
      R3 = R*R*R                                                                
      DRRG = (DEC1+DEC2*R+DEC3*R3)*PI/(6.D0)                                                  
      RETURN
C     The following extrapolates the density derivative of G(R) past r=6
C     empirically.
 1011 RHO=ETA*6.D0/PI
      A=3.64779052-4.5326958*RHO+2.8254675*RHO**1.5D0-1.1230973
     1*DLOG(RHO)/RHO+0.02371867*DLOG(RHO)/RHO**2
      B=7.0399098-13.686034*RHO-0.85088666*RHO**3+2.9323443
     1*DEXP(RHO)+6.6465928*DSQRT(RHO)*DLOG(RHO)
      C=5.9382949+3.7429435*RHO+1.108728*RHO**2*DLOG(RHO)-0.96469546
     1*DEXP(RHO)+0.19983335*DLOG(RHO)
      D=1.1350252-1.9259676*RHO-1.8835841*RHO**2*DLOG(RHO)+1.0381572
     1*RHO**3+0.00007112*DLOG(RHO)
      DA=2.2578288-2.6419506/RHO-0.57462084/RHO**2-0.49861017
     1/RHO**3+0.03826629/RHO**4
      DB=-1.1150611+1.5019557*RHO-1.5365115*RHO**1.5D0+0.79801674
     1*RHO**3-1.2700068/RHO
      DC=2.5339359-5.0383713*RHO+6.0947248*RHO**1.5D0-1.3679322
     1*RHO**2.5D0+0.20489549/RHO
      DD=-18.494398+17.445882*RHO-16.716728*DSQRT(RHO)*DLOG(RHO)
     1+0.45062181/RHO**1.5D0-0.12072079/RHO**2
      DRRG=DEXP(-B*R)*(DA*DSIN(C*R+D)+DCOS(C*R+D)*A*(DC*R+DD)-R*DB*DSIN
     1(C*R+D)*A)/R                                                                                                                                                        
      RETURN
C     The following calculates the density derivative of G(R) analytically 
C     for 1<=R<=6.                                                                    
 1010 DRRG = 0.D0                                                                
C                                                                               
C         FIRST SHELL CONTRIBUTION                                              
C                                                                               
      EX1=EX11*DEXP(R*IT1)                                                      
      EX2=EX12*CDEXP(R*IT2)
      FN11=IG101
      FN12=IG102                                              
      DRRG=DRRG+EX1*(DF1DE+(R-1.D0)*FN11*DTE)+2.D0*EX2*(DIF1DE+
     1(R-1.D0)*FN12*DITE)
      IF(R.LE.2.D0)GO TO 99                                                     
C                                                                               
C         SECOND SHELL CONTRIBUTION                                             
C                                                                               
      EX1=EX1*EX11                                                              
      EX2=EX2*EX12
      FN21=IG201*(IG211+IG221*(R))
      FN22=IG202*(IG212+IG222*(R))
      DRRG=DRRG+EX1*(DF2DE+(R-2.D0)*FN21*DTE)+2.D0*EX2*(DIF2DE+(R-2.D0)
     1*FN22*DITE)                                                             
      IF(R.LE.3.D0)GO TO 99 
C                                                                               
C         THIRD SHELL CONTRIBUTION                                             
C                                                                               
      EX1=EX1*EX11                                                              
      EX2=EX2*EX12
      FN31=IG301*(IG311+IG321*(R-3.D0)+IG331*(R-3.D0)**2.D0)
      FN32=IG302*(IG312+IG322*(R-3.D0)+IG332*(R-3.D0)**2.D0)
      DRRG=DRRG+EX1*(DF3DE+(R-3.D0)*FN31*DTE)+2.D0*EX2*(DIF3DE+(R-3.D0)
     1*FN32*DITE)                                                                
      IF(R.LE.4.D0)GO TO 99 
C                                                                               
C         FOURTH SHELL CONTRIBUTION                                             
C                                                                               
      EX1=EX1*EX11                                                              
      EX2=EX2*EX12
      FN41=IG401*(IG411+IG421*(R-4.D0)+IG431*(R-4.D0)**2.D0+IG441*(R-
     14.D0)**3.D0)
      FN42=IG402*(IG412+IG422*(R-4.D0)+IG432*(R-4.D0)**2.D0+IG442*(R-
     14.D0)**3.D0)
      DRRG=DRRG+EX1*(DF4DE+(R-4.D0)*FN41*DTE)+2.D0*EX2*(DIF4DE+(R-4.D0)
     1*FN42*DITE)                                                             
      IF(R.LE.5.D0)GO TO 99
C                                                                               
C         FIFTH SHELL CONTRIBUTION                                             
C                                                                               
      EX1=EX1*EX11                                                              
      EX2=EX2*EX12
      FN51=IG501*(IG511+IG521*(R-5.D0)+IG531*(R-5.D0)**2.D0+IG541*(R-
     15.D0)**3.D0+IG551*(R-5.D0)**4.D0)
      FN52=IG502*(IG512+IG522*(R-5.D0)+IG532*(R-5.D0)**2.D0+IG542*(R-
     15.D0)**3.D0+IG552*(R-5.D0)**4.D0)
      DRRG=DRRG+EX1*(DF5DE+(R-5.D0)*FN51*DTE)+2.D0*EX2*(DIF5DE+(R-5.D0)
     1*FN52*DITE)                                                             
      IF(R.LE.6.D0)GO TO 99
   99 DRRG=PI*DRRG/(6.D0*R)
      RETURN                                                                    
      END 
