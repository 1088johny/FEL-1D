PROGRAM      Fortran
     IMPLICIT     NONE
  INTEGER i, U
  INTEGER, PARAMETER:: NSTEPS=40000, POINTS=1000
  INTEGER, PARAMETER:: vp=16 !vp=8 (precisão dupla); vp=16 (precisão quadrupla)
  REAL(KIND=vp), PARAMETER::  pi = 3.14159265358979e0_vp, c=2.99792458e10_vp,  me=9.1093897e-28_vp,   q = -4.8032e-10_vp, c2 = c*c
  REAL(KIND=vp), PARAMETER :: q2 = q*q, bm=-1.0e0_vp, rc= 1.0e0_vp, schwig=1.0e0_vp, schrad=1.0e0_vp
  REAL(KIND=vp), PARAMETER :: schplsm=1.0e0_vp, gperpon=1.0e0_vp, fermi=0.0e1_vp, reldetunningon=0.0e1_vp, detunning=0.0e1_vp
  REAL(KIND=vp), PARAMETER :: damping=0.0e0_vp, hbar = (6.6260755e-27_vp)/(2.0e0_vp*pi), lambdaccut = hbar/(me*c)
  REAL(KIND=vp):: nb, n0, as0, aw0, lambdaw, lambdas, ks, ks2, wpe, wpe2, wpen, ws2, ww, ww2, gamae, ws, theta0, l, r, r2
  REAL(KIND=vp):: sigma, alfa
  REAL(KIND=vp):: chiini, chiend, h, kks
  REAL(KIND=vp), PARAMETER :: asscale = 1.0e0_vp, awscale=1.0e0_vp, nscale1=1.0e5_vp
  REAL(KIND=vp), PARAMETER :: as2scale = 1.0e0_vp, aw2scale=1.0e0_vp, nscale2=1.0e0_vp
  REAL(KIND=vp), PARAMETER :: as2scalei = 1.0e12_vp, aw2scalei=1.0e12_vp, nscalei=1.0e0_vp
  REAL(KIND=vp), PARAMETER :: epsonchi=1.0_vp
 
  lambdaw = 1.26e-5_vp !cm (126 nm)
  gamae=5.0e0_vp
 ! gamae=1.0000000000000001e0_vp
 nb=4.1e18_vp
!  nb=5.8e22_vp
 
 !wpen= 1.125e-6_vp !wpen = wpe/c
! r = 1.0e1  
 ww = 2.0e0_vp*pi*c/lambdaw
! ww2 = ww*ww
! r = ww/(c*wpen) !r = ww/wpe
! r2 = r*r 
! nb = (me*ww2)/(r2*q2*pi*4.0e0_vp)
 n0 = 1.0e-13_vp
 !n0 = 1.0e-8_vp
 as0 = 6.0e-15_vp
 aw0 = 8.5e-5_vp

!****************INITIAL AND FINAL POSITION-block result2PW-qfeldet2nu0.dat -bxy 1:2**********************
  chiini = 0.0e1_vp
  chiend = 5.0e6_vp
!***********************************************lamp1****************

 h = (chiend - chiini)/NSTEPS
! open(unit=17,file='inicialization-laseramp-f2.dat')
! open(unit=18,file='result-laseramp-f2.dat')
! open(unit=19,file='result2-laseramp-f2.dat')
! open(unit=20,file='result2PW-laseramp-f2.dat')
! open(unit=21,file='brilouin-laseramp-f2.dat')
! open(unit=24,file='theta-laseramp-f2.dat')
! open(unit=22,file='detunning-laseramp-f2.dat')
! open(unit=23,file='rescondition-laseramp-f2.dat')
! open(unit=25,file='energy-laseramp-f2.dat')

 write(*,*) 'condition for validity = ', 32.0e0_vp*pi*lambdaccut*gamae
! open(unit=18,file='graphqfl1aw-s-new.dat')
! open(unit=19,file='graphqfel1aw2-s-new.dat')
! open(unit=20,file='graphqfel1aw2p-s-new.dat')

!WRITE(18,*) '****************************************************************************************************'

  
!150 FORMAT(A,E)
  CALL INICIALIZATION(nb, ww, gamae, as0, aw0, chiend)
 
!  CALL EXECUTION(chiini,chiend,n0,as0,aw0,h)
 
! lambdas = 2.0e0_vp*pi/ks

  CALL DATAOUTPUT()
 
!130 FORMAT(a, D17.11, a, D17.11, a, D17.11)
!************************DEFINICAO DE FUNCOES************************************
CONTAINS

  
  !EQUACOES DO MOVIMENTO
  REAL(KIND=vp) FUNCTION fnd(chi,n,as,aw,gperp,gz,theta)
    IMPLICIT NONE
    REAL(KIND=vp) chi, n,as,aw,gperp,gz, gz2, gz3, gz4, gperp2, gperp3, gperp4,theta, dl, ds, dw, deltal,  wpe
    REAL(KIND=vp) S0, S1, S2, damping_norm
    REAL(KIND=vp) A0, A1, A2, D1, D2
   
    COMMON/blmismatching/dl, deltal
    COMMON/constants/A0, A1, A2
    COMMON/denominators/D1, D2
    COMMON/blsignal/S0, S1, S2
    COMMON/beamparameters/wpe
    COMMON/damp/damping_norm
    
    

    gperp2 = gperp*gperp
    gperp3 = gperp2*gperp
    gperp4 = gperp2*gperp2
    gz2=gz*gz
    gz3=gz2*gz
    gz4=gz2*gz2

    !write(*,*) 'estoy aquí damp  =', damping_norm
    !call exit(0)

    fnd = A0*aw*as*cos(theta) - n*damping_norm
    
    !fnd = A0*aw*as*cos(theta) - n*damping_norm
    
  END FUNCTION fnd


  REAL(KIND=vp) FUNCTION fas(chi,n,as,aw,gperp,gz,theta)
    IMPLICIT NONE
    REAL(KIND=vp) chi, n,as,aw,gperp,gz, gperp2, gperp3, gperp4,theta, dl, ds, dw
    REAL(KIND=vp) S0, S1, S2
    REAL(KIND=vp) A0, A1, A2, D1, D2
    COMMON/constants/A0, A1, A2
    COMMON/blsignal/S0, S1, S2
    COMMON/denominators/D1, D2
    
    
   ! fas = -2*gz*gz*A1*aw*n*cos(theta)
    fas = -A1*aw*n*cos(theta)

  END FUNCTION fas

  REAL(KIND=vp) FUNCTION faw(chi,n,as,aw,gperp,gz,theta)
    IMPLICIT NONE
    REAL(KIND=vp) chi, n,as,aw,gperp,gz, gperp2, gperp3, gperp4,theta
    REAL(KIND=vp) S0, S1, S2
    REAL(KIND=vp) A0, A1, A2, D1, D2
    COMMON/constants/A0, A1, A2
    COMMON/denominators/D1, D2
    COMMON/blsignal/S0, S1, S2

    gperp2 = gperp*gperp
    gperp3 = gperp2*gperp
    gperp4 = gperp2*gperp2

    faw = A2*n*as*cos(theta)
    !faw = 0.5*A2*n*as*cos(theta)
    !faw = (A2/(D2*gperp*gz))*n*as*cos(theta)
     !faw = 0.0e1_vp

  END FUNCTION faw
 

 REAL(KIND=vp) FUNCTION fth(chi,n,as,aw,gperp,gz,theta)
    IMPLICIT NONE
    REAL(KIND=vp) chi, n,as,aw,gperp,gz,gperp2, gperp3, gperp4, theta, dl, deltal
    REAL(KIND=vp) S0, S1, S2, damping_norm
    REAL(KIND=vp) A0, A1, A2, D1, D2
    REAL(KIND=vp) wl, kw, kl, ws2, ww2, ks, gz2, gz3, gz4, wpe, wpe2, ve, betae9
    COMMON/constants/A0, A1, A2
    COMMON/denominators/D1, D2
    COMMON/blsignal/S0, S1, S2
    COMMON/blmismatching/dl, deltal
    COMMON/fields/wl, kw, kl, ws2, ww2, ks
    COMMON/beamparameters/gz2, wpe, wpe2, ve
   REAL(KIND=vp)  aux1, aux2, aux3 
  
    ww = sqrt(ww2)
    ws = sqrt(ws2)
    gz = sqrt(gz2)
    gz2=gz*gz
    gz3=gz2*gz
    gz4=gz2*gz2
    !dl = F(ws, ww, gz, ve, wpe, nb, gperp)

aux1 = A0*(aw*as/n)
!aux2 = A1*(aw*n/as)*2*gz*gz
!aux3 = 0.5*A2*(as*n/aw)

aux2 = A1*(aw*n/as)
aux3 = A2*(as*n/aw)

!write(*,*) 'estoy aquí delta l =', deltal
!call exit(0)
!fth = - 2*gz*gz*deltal + (-aux1 + aux2 + aux3 )*sin(theta) 
!fth = (detunning -(aux1 + aux2 + aux3 )*sin(theta) )
 fth =  deltal - (aux1 - aux2 + aux3 )*sin(theta) 

 END FUNCTION fth

  !SUBROTINAS

  !SUBROTINA DE INICIALIZACAO


REAL(KIND=vp) FUNCTION F(ws, ww, gz, ve, wpe, nb, gperp)
   IMPLICIT NONE
   INTEGER*8 i, j
   REAL(KIND=vp) gz, gz2, gz3, ve, ks, ks2, kw, ws, ww, wl, wl2, Omegal2, Omegab, Omegal, ww2, wpe2, kl
   REAL(KIND=vp) kl2, sigma, sigma2, ws2, wpe, gperp,gperp2, gperp3, gperp4, betae, vf, vf2, kf, nb, kks
   COMMON/fields/wl, kw, kl, ws2, ww2, ks
   COMMON/beammode/Omegab
   COMMON/counterbl/i, j
   
	
   gz2 = gz*gz
   gz3 = gz2*gz
   gperp2 = gperp*gperp
   gperp3 = gperp2*gperp
   gperp4 = gperp2*gperp2
   ws2 = ws*ws
   ww2 = ww*ww
   wpe2 = wpe*wpe
   betae = ve/c
   kw = (ww/c)*sqrt(1.0e0_vp - wpe2/(ww2*gz))
   ks = (ws/c)*sqrt(1.0e0_vp - wpe2/(ws2*gz))
   kks = ks
   kl = ks + kw
   kl2 = kl*kl
   wl = ws - ww - wpe*detunning
   sigma = (hbar*kl/(me*gz3))
   sigma2 = sigma*sigma
   !vf2 = ((3.0e0_vp*pi*pi*nb)**(1.0e0_vp/3.0e0_vp))*(hbar*c)/(3.0e0_vp*me*gperp3*gz3) !relativistic fermi
   vf = (hbar/me)*(3.0e0_vp*pi*pi*nb)**(1.0e0_vp/3.0e0_vp)
   vf2 = vf*vf
   
  
   Omegab = bm*sqrt(wpe2/gz3  + sigma2*kl2/4.0e0_vp)
    
  ! F = -ww + ws - kl*ve - Omegab !LASER AMPLIFICATION
    F = wl - kl*ve - Omegab  !QUANTUM FEL
  
      
  END FUNCTION F
   
  SUBROUTINE INICIALIZATION(nb, ww, gamae, as0, aw0, chiend)

    IMPLICIT NONE
    INTEGER*8 I, J, AUX, L
    REAL(KIND=vp) A0, A1, A2, A0S, A1S, A2S, A0T, A1T, A2T, D1, D2, S0, S1, S2
    REAL(KIND=vp) lambdas, lambdaw, aw0, as0, theta0, chiend, gperp0, gperp02, gperp03, gperp04
    REAL(KIND=vp) ws, ww, wl, wl2, kl2, ks, ks2, kw2, kw, kl, ws2, ww2, k, w, ke, vph, vph2, betaph, betaph2, betaphi
    REAL(KIND=vp) Omegas, Omegaw, Omegal, Omegaq2, qs, qw, ql, b, b2, wpeold, wpeini, wpeend, dwpe, pe, wsend, wsini, dws, wsold
    REAL(KIND=vp) deltal, lambdacutc, damping_norm, omegaq
    REAL(KIND=vp) nb, wpe, wpe2, Omegape, Omegape2, ve, betae, betae2, betas, betas2, alfa, alfa2, alfaw, alfaw2, betaw, betaw2
    REAL(KIND=vp) betaq, betar, betar2, betaf
    REAL(KIND=vp) gamae, gz, gz2, gz3, gz4, dl, dwsddelta, deltaws
    REAL(KIND=vp) sigma, sigma2, sigmanew, sigmanew2, FAUXOLD, FAUX, kwold, kwini, kwend, kf, vf, vf2
    REAL(KIND=vp) Deltaq, Deltas, Deltaw, OMEGALOLD, OMEGAL2OLD, OMEGAL2, Omegab, PM1, PM2, PM3, PM4, I0S, I0W
    
   COMMON/blmismatching/dl, deltal
    COMMON/blsignal/S0, S1, S2
    COMMON/fields/wl, kw, kl, ws2, ww2, ks
    COMMON/fieldswp/Omegas, Omegaw, Omegal, qs, qw, ql, betas, betaw
    COMMON/beamparameters/gz2, wpe, wpe2, ve, betae
    COMMON/constants/A0, A1, A2
    COMMON/denominators/D1, D2
    COMMON/beammode/Omegab
    COMMON/counterbl/ i, j
    COMMON/fauxbl/FAUX
    COMMON/bltheta0/theta0
    
      
  gz = gamae
    gz2 = gz*gz
    gz3 = gz2*gz
    gz4 = gz2*gz2
 
   gperp0= 1.0e0_vp 
   gperp02= gperp0*gperp0
   gperp03= gperp02*gperp0

    betae2 = 1.0e0_vp - (1.0e0_vp/gz2)
    betae  = sqrt(betae2)
    ve = c*betae
    pe = gz*me*ve
    ke = pe/hbar
    wpe2 = (q2*pi*4.0e0_vp*nb)/me
    wpe = sqrt(wpe2)

    
     wsini = ww
     wsend = 4.0e1_vp*ww*gz2

    aux= 0
    
    DO J=1,100

    dws = (wsend - wsini)/(4.0e5_vp)
    ws = wsini
    FAUX = F(ws, ww, gz, ve, wpe, nb, gperp0)

       
    DO I=0, 400000

      	wsold = ws
	FAUXOLD = FAUX
    	
    	ws = ws + dws
        FAUX = F(ws, ww, gz, ve, wpe, nb, gperp0)
         
     
 	IF (FAUXOLD/FAUX < 0.0e0_vp) THEN
!    		WRITE(*,*) '[',J,I,']'
!		WRITE(*,*) 'EXISTE RAIZ ENTRE '
!    		WRITE(*,*) wsold, '<WS<', ws
!		WRITE(*,*) '|ABS(FAUX)-ABS(FAUXOLD)| =', ABS(ABS(FAUX)-ABS(FAUXOLD))
!		WRITE(*,*) ' FAUX =', FAUX
!		WRITE(*,*)
    		wsini = wsold
    		wsend = ws
                aux = 1
		EXIT
    	END IF
      
    	     
    END DO

	IF (AUX .EQ. 0) THEN
	WRITE(*,*) 'PARAMETROS INICIAIS NÃO PODEM SATISFAZER  A RELAÇÃO DE DISPERSÃO OMEGAL2 = WPE2/GAMAE3 + SIGMA2 KL2'
        STOP
	END IF

    IF ( ABS( (ABS(FAUX)-ABS(FAUXOLD)) ) .LE. 1.0e-30_vp) THEN
    !write(*,*) ABS( (ABS(FAUX)-ABS(FAUXOLD)))
    deltal = wpe*detunning
    !WRITE(*,*) deltal
    !WRITE(*,*) 'Detunning =', deltal
    EXIT
    END IF


 
    END DO

    ws = (wsend + wsini)/2.0e0_vp
    ws2 = ws*ws
  
  

    WRITE(*,*) 'initial relativistic detunning =           ',  F(ws, ww, gz, ve, wpe, nb, gperp0)
   
    WRITE(*,*) 'normalized inicial relativistic detunning =', F(ws, ww, gz, ve, wpe, nb, gperp0)/wl
   ! theta0 = 0.0e1_vp
     theta0 = pi/2.0e0_vp

    kl2 = kl*kl
    sigma = hbar*kl/(gz3*me)
    sigma2 = sigma*sigma
    betaw = c*kw/ww
    betaw2 = betaw*betaw
    betas = c*ks/ws
    betas2 = betas*betas
    Omegal = wl - kl*ve
    betaq = (c*lambdaccut*lambdaccut*kl2*kl)/(2.0e0*gz3*gz3*Omegab)
    lambdas = 2.0e0_vp*pi/ks
    lambdaw = 2.0e0_vp*pi/kw
    vph = wl/kl
    vph2 = vph*vph
    betaph = vph/c
    betaph2 = betaph*betaph
    vf = (hbar/me)*(3.0e0_vp*pi*pi*nb)**(1.0e0_vp/3.0e0_vp)
    vf2 = vf*vf
    
    
    dwsddelta = 1.0e0_vp/(1.0e0_vp - (ws/(c2*ks))*(ve + bm*sigma2*kl/Omegab) ) 
    deltaws = dwsddelta*detunning*wpe2

   
    WRITE(*,*) '#', ' ws/ww =', ws/ww
    WRITE(*,*) '#', ' ws/wl =', ws/Omegal
    WRITE(*,*) '#', ' ks/kw =', ks/kw
    WRITE(*,*) '#', ' wq2/wpe2 =', (sigma2*kl2*gz3/4.0e0_vp)/wpe2

    write(*,*) 'lambdas   =', lambdas
    write(*,*) 'dettuning     =', detunning
    write(*,*) 'diflbds =',( 2.0e0_vp*pi*c/ws2)*dwsddelta
    write(*,*) 'ratio q/p =', (hbar*hbar*(ks+kw)*(ks+kw)*(ks+kw)*(ks+kw)/(4.0e0_vp*me*me*gz3*gz3))/(wpe2/gz3)
      
    !D1 = 1.0e0_vp - betas
    !D2 = 1.0e0_vp + betaw
 !    D1 = ( betas - betae)  
 
    betaphi= 1.0e0_vp

    PM1 = (hbar*kl)/(me*gz*c)
    PM2 = (sigma2*kl2)/(4.0e0_vp)
    PM3 = ((3.0e0_vp/5.0e0_vp))*(vf2*kl2)
    PM4 = wpe2/(gz3*gperp03)

   
     I0S = (2.0e0_vp*1.0e12_vp)/((8.6e0_vp)*(8.6e0_vp)*(lambdas)*(lambdas))
     I0W = (2.0e0_vp*1.0e12_vp)/((8.6e0_vp)*(8.6e0_vp)*(lambdaw)*(lambdaw))

  write(*,*) 'IOS   =', I0S
    write(*,*) 'IOW     =', I0W
!****************TEMPORAL LIMIT*******************************************
!    A0T = c2*kl2/(4.0e0_vp*gz4*wpe*Omegab*betaq)
     A0T = c2*kl2/(4.0e0_vp*gz4*wpe*Omegab)
    A1T = wpe/(2*gz*ws)
    A2T = wpe/(2*ww*gz)

  !  A0S = ((kl2*c2)/(4.0e0_vp*wpe*Omegab))
  !  A1S = -(wpe/(2.0e0_vp*ws))
  !  A2S =  wpe/(2.0e0_vp*ww)


    A0 = A0T	
    A1 = A1T
    A2 = A2T

    write(*,*) 'CS   =', A1T
    write(*,*) 'CD   =', A0T
    write(*,*) 'bs   =', betas
    write(*,*) 'be   =', betae
    write(*,*) 'bq   =', betaq
 

!*************************************************************************



   DO U=17, 25
   
   WRITE(U,*)   '#', '  CONDICOES INICIAIS '
   WRITE(U,110) '#',  ' n(0)=', n0, ' as(0)=', as0, ' aw(0)=', aw0
   WRITE(U,*)   '#'
   WRITE(U,110) '#', '  lbds= ', lambdas, '  lbdw= ', lambdaw, '  nb= ', nb
   WRITE(U,110) '#', '  ks= ', ks, '  kw= ', kw, '  ks+kw= ', ks+kw
   WRITE(U,112) '#', '  ww/wpe= ', ww/wpe, '  ws/wp= ', ws/wpe
   WRITE(U,110) '#', '  wpe2= ', wpe2, '  wpe= ', wpe, '  gz=', gz
   WRITE(U,112) '#', '  wpe2/(gz*gperp0)= ', wpe2/(gz*gperp0), '  sqrt(wpe2/(gz*gperp0))= ', sqrt(wpe2/(gz*gperp0))
   WRITE(U,111) '#', '  ww    =    ', ww
   WRITE(U,111) '#', '  ws    =    ', ws
   WRITE(U,111) '#', '  ws-ww =    ', ws-ww
   WRITE(U,111) '#', '  wl    =    ', wl
   WRITE(U,111) '#', '  klve  =    ', kl*ve
   WRITE(U,111) '#', '  wl-klve=   ', wl - kl*ve
   WRITE(U,111) '#', '  wl-klve - Omegab =', wl - kl*ve - Omegab
   WRITE(U,111) '#', '  Omegal - Omegab    ', Omegab
   WRITE(U,111) '#', ' Omegab+klve=', Omegab + kl*ve
   WRITE(U,111) '#', '  Norm. Detunning delta/wpe =', detunning
   WRITE(U,111) '#', '  Non Norm. Detunning delta           =', detunning*wpe
   WRITE(U,111) '#', '  Damping                             =', damping
   WRITE(U,111) '#', '  ws(Detunning = delta) =', ws 
   WRITE(U,111) '#', '  ws(Detunning = 0)     ~', ws - dwsddelta*detunning*wpe
   WRITE(U,111) '#', '  betae= ', betae
   WRITE(U,111) '#', '  betaq= ', betaq
   WRITE(U,111) '#', '  betas= ', betas
   WRITE(U,111) '#', '  betaw= ', betaw
   WRITE(U,111) '#', '  vf = ', vf
   WRITE(U,111) '#', '  PM1 = hbar*kl/ge*me*c =', PM1
   WRITE(U,111) '#', '  PM2 = sigma2*kl2/4 =   ', PM2
   WRITE(U,111) '#', '  PM3 = 3/5*vf2*kl2 =    ', PM3
   WRITE(U,111) '#', '  PM4 = wpe2/ge3=        ', PM4
   WRITE(*,*)   '#', '  RATIO PM2/PM4 =        ', PM2/PM4
   WRITE(U,111) '#', '  RATIO PM2/PM4 =        ', PM2/PM4
   WRITE(U,111) '#', '  RATIO PM3/PM4 =        ', PM3/PM4
   WRITE(U,111) '#', '  RATIO PM2/PM3 =        ', PM2/PM3
   WRITE(U,111) '#', '  4 kwlambdacgamae =     ', 1.6e1_vp*lambdaccut*kw*gz
   WRITE(U,111) '#', '  Modo LAMBDACCUT    =        ', -2.0e0_vp*pi*(lambdaccut/gz)                           
 !  WRITE(U,110) '#', '  lambdae =', 2.0e0_vp*pi*hbar/pe, '  lambdac =',  2.0e0_vp*pi*hbar/(me*c), '  RATIO lambdae/lambdac =', (2.0e0_vp*pi*hbar/pe)/(2.0e0_vp*pi*hbar/(me*c))
   WRITE(U,110) '#', '  ke =', ke, '  kl =',  kl, '  RATIO kl/ke =', kl/ke
   WRITE(U,110) '#', '  C0= ', A0, '  C1= ', A1, '  C2= ', A2
   WRITE(U,111) '#', '  phase velocity/c (wl/ckl): ',  wl/(c*kl)
   WRITE(U,111) '#', '  phase velocity/c (Omegal/ckl)[1]: ',  Omegal/(c*kl)
!   WRITE(U,111) '#', '  phase velocity/c (Omegal/ckl)[2]: ',  sqrt((wpe2/(c2*kl2) + (lambdaccut*lambdaccut*kl2)/(4.0e0_vp*gz3) )/gz3)
   WRITE(U,111) '#', '  wl2/c2 =', wl*wl/c2
   WRITE(U,111) '#', '  kl2 =', kl*kl
   WRITE(U,111) '#', '  wl2/c2k2 =', (wl*wl)/(c2*kl2)
   WRITE(U,110) '#', '  Gama : ', gz*gperp0, '  Gama phase: ', sqrt(1.0e0_vp/(1.0e0_vp - betaph2)), '  Gama/Gama phase: ', gz*gperp0/(sqrt(1.0e0_vp/(1.0e0_vp - betaph2)))
   !WRITE(U,111) '#', '  Initial detunning =           ',  F(ws, ww, gz, ve, wpe, nb, gperp0)
   !WRITE(U,111) '#', '  normalized inicial detunning =', F(ws, ww, gz, ve, wpe, nb, gperp0)/wl
   WRITE(U,110) '#', '  Power Is Scale=', as2scalei, '  Power Iw Scale=', aw2scalei, 'plasmon scale', nscalei
   WRITE(U,110) '#', '  as2 Scale=', as2scale, '  aw2 Scale=', aw2scale, 'plasmon scale', nscale2
   WRITE(U,110) '#', '  as Scale=', asscale, '  aw Scale=', awscale, 'plasmon scale', nscale1
   WRITE(U,110) '#', '  horizontal scale= ', epsonchi
   WRITE(U,111) '#', '  Distance (z) of interaction =  ', epsonchi*c*chiend/wpe
   WRITE(U,111) '#', '  Time     (t) of interaction =  ', epsonchi*chiend/wpe
   WRITE(U,111) '#', '  lambdas(vac) =', 2.0e0_vp*pi*c/ws
   WRITE(U,111) '#', '  lambdas      =', 2.0e0_vp*pi/ks
   WRITE(U,111) '#', '  lambdaw(vac) =', 2.0e0_vp*pi*c/ww
   WRITE(U,111) '#', '  lambdaw      =', 2.0e0_vp*pi/kw
   WRITE(U,111) '#', '  (wp~0, quantum term~0) lbdw=lbs*4gz2 =', (2.0e0_vp*pi*c/ws)*4.0e0_vp*gz2*gperp02
   WRITE(U,113) '#', '  schrad =', schrad, '  schwig= ', schwig, '  schplsm =', schplsm, '  mode= ', bm
   WRITE(U,113) '#', '  rc= ', rc, '  gperpon= ', gperpon, '  fermi= ', fermi
   

   END DO

   WRITE(*,*) 'max as2 calc =', (ww/ws)*aw0*aw0
   WRITE(*,*) 'max n2  calc =', ((c2*kl*kl*ww)/(2.0e0_vp*gz3*abs(Omegab)*wpe2))*aw0*aw0


   110 FORMAT(a, a, D17.11, a, D17.11, a, D17.11)
   111 FORMAT(a, a, D36.30)
   112 FORMAT(a, a, D24.18, a, D24.18)
   113 FORMAT(a, a, D6.1, a, D6.1, a, D6.1, a, D6.1)
   120 FORMAT(I7, 4X, D17.11, 4X, D17.11, 4X, D17.11)
  ! 130 FORMAT(a, D17.11, a, D17.11, a, D17.11) 
 !  140 FORMAT(I, 4x, D17.11, 4X, D17.11)
 !  150 FORMAT(A,E)

END SUBROUTINE INICIALIZATION


!****************************************************************************


SUBROUTINE EXECUTION(chi0,chiend,n0,as0,aw0,h)
  IMPLICIT NONE
  INTEGER i
  REAL(KIND=vp), DIMENSION(1:4):: RKFND, RKFAS, RKFAW, RKFTH
  REAL(KIND=vp) chi,n,as,aw,theta, h, wl, ks, kw, kl, ws2, ww2, ww, ws,lambdas, lambdaw
  REAL(KIND=vp) n0rn, as0rn, aw0rn, theta0
  REAL(KIND=vp) chi0, n0, as0, aw0, chiend
  REAL(KIND=vp) A0, A1, A2, gperp, gz3, gz4, D1, D2, I0W, I0S
  REAL(KIND=vp) dnd, das, daw, dth, gperp0, gperp02, gperp04
  REAL(KIND=vp) gz2, gz, wpe, wpe2, ve, betae, damping_norm
  
  COMMON/dynvariablesrn0/n0rn, aw0rn, as0rn
  COMMON/constants/A0, A1, A2
  COMMON/fields/wl, kw, kl, ws, ks, ww
  COMMON/beamparameters/gz2, wpe, wpe2, ve, betae
  COMMON/relativisticfactor/gperp
  COMMON/denominators/D1, D2
  COMMON/damp/damping_norm
    
    
!COMMON BLOCKS WITH OUTPUT SUBROUTINE
 COMMON/dynvariables/n, aw, as, theta
 COMMON/space/chi

!COMMON BLOCKS WITH INICIALIZATION SUBROUTINE
 COMMON/bltheta0/theta0


 ! ww = sqrt(ww2)
 ! ws = sqrt(ws2)
  gz = sqrt(gz2)
  gz4 = gz2*gz2
  
  
  chi = chi0
  n = n0
  as = as0
  aw = aw0
  theta = theta0

  

  write(*,*) 'detunning = ', detunning
  write(*,*) 'Constant couplings'
  write(*,*) 'C0= ', A0
  write(*,*) 'C1= ', A1
  write(*,*) 'C2= ', A2
  write(*,*) '******************************************'
  write(*,*) 'initial fields'
  write(*,*) 'n= ', n
  write(*,*) 'as= ', as
  write(*,*) 'aw= ', aw
  write(*,*) '******************************************'
  write(*,*) 'A0*aw*as= ', A0*aw*as
  write(*,*) 'A1*aw*n= ', A1*aw*n
  write(*,*) 'A2*as*n= ', A2*as*n
  write(*,*) '******************************************'
  write(*,*) 'A0*aw*as/n= ', A0*aw*as/n
  write(*,*) 'A1*aw*n/as= ', A1*aw*n/as
  write(*,*) 'A2*as*n/aw= ', A2*as*n/aw

  

  WRITE(19,111) '#'
  WRITE(19,111) '#', 'chi',   'as^2',         'aw^2',   'n',    'theta'

!  WRITE(18,110) '#', '  nrn(0)" =', n0rn, '  asrn(0)" =', as0rn, '  awrn(0)" =', aw0rn
  WRITE(18,111) '#'
  WRITE(18,111) '#', 'chi',  'as^2',         'aw^2',   'n',    'theta'

  WRITE(20,111) '#'
  WRITE(20,111) '#', 'chi',  'as^2',         'aw^2',   'n',    'theta'


  110 FORMAT(a, a, D17.11, a, D17.11, a, D17.11)
  111 FORMAT(a, a, 16X, a, 16X, a, 16X, a, 16X, a, 16X)

  gperp = sqrt(1.0e0 + as*as + aw*aw)*gperpon + 1.0e0_vp - gperpon   
  CALL OUTPUT (0)

  DO i=1, NSTEPS

	 !write(*,*) 'detunning =  ', F(ws, ww, gamae, ve, wpe, gperp)       
    
         RKFND(1)  = h*fnd(chi,n,as,aw,gperp,gz,theta)
         RKFAS(1)  = h*fas(chi,n,as,aw,gperp,gz,theta)
         RKFAW(1)  = h*faw(chi,n,as,aw,gperp,gz,theta)
         RKFTH(1)  = h*fth(chi,n,as,aw,gperp,gz,theta)
	
         RKFND(2)  = h*fnd(chi+h/2.0e0_vp,n+RKFND(1)/2.0e0_vp,as,aw,gperp,gz,theta)
         RKFAS(2)  = h*fas(chi+h/2.0e0_vp,n,as+RKFAS(1)/2.0e0_vp,aw,gperp,gz,theta)
         RKFAW(2)  = h*faw(chi+h/2.0e0_vp,n,as,aw+RKFAW(1)/2.0e0_vp,gperp,gz,theta)
         RKFTH(2)  = h*fth(chi+h/2.0e0_vp,n,as,aw,gperp,gz,theta+RKFTH(1)/2.0e0_vp)

         RKFND(3)  = h*fnd(chi+h/2.0e0_vp,n+RKFND(2)/2.0e0_vp,as,aw,gperp,gz,theta)
         RKFAS(3)  = h*fas(chi+h/2.0e0_vp,n,as+RKFAS(2)/2.0e0_vp,aw,gperp,gz,theta)
         RKFAW(3)  = h*faw(chi+h/2.0e0_vp,n,as,aw+RKFAW(2)/2.0e0_vp,gperp,gz,theta)
         RKFTH(3)  = h*fth(chi+h/2.0e0_vp,n,as,aw,gperp,gz,theta+RKFTH(2)/2.0e0_vp)
         
         RKFND(4)  = h*fnd(chi+h,n+RKFND(3),as,aw,gperp,gz,theta)
         RKFAS(4)  = h*fas(chi+h,n,as+RKFAS(3),aw,gperp,gz,theta)
         RKFAW(4)  = h*faw(chi+h,n,as,aw+RKFAW(3),gperp,gz,theta)
         RKFTH(4)  = h*fth(chi+h,n,as,aw,gperp,gz,theta+RKFTH(3))
         
         dnd  = (1.0e0_vp/6.0e0_vp)*(RKFND(1)+2.0e0_vp*RKFND(2)+2.0e0_vp*RKFND(3)+RKFND(4))
         das  = (1.0e0_vp/6.0e0_vp)*(RKFAS(1)+2.0e0_vp*RKFAS(2)+2.0e0_vp*RKFAS(3)+RKFAS(4))
         daw  = (1.0e0_vp/6.0e0_vp)*(RKFAW(1)+2.0e0_vp*RKFAW(2)+2.0e0_vp*RKFAW(3)+RKFAW(4))
         dth  = (1.0e0_vp/6.0e0_vp)*(RKFTH(1)+2.0e0_vp*RKFTH(2)+2.0e0_vp*RKFTH(3)+RKFTH(4))

         n  =  n + dnd
         as = as + das
         aw = aw + daw
        
!         if(i>10) then
!         stop
!         end if
!         WRITE(*,*) ' n=', dnd, ' as=', das,' aw=',daw, ' theta=',dth
!         write(*,*) 'theta, dth 1 =', theta, dth
         theta = theta + dth
!         write(*,*) 'theta, dth 2 =', theta, dth
!         exit

         chi = chi + h !IT IS THE SPACE POINT RELATED TO THE STEP ABOVE
         gperp = sqrt(1.0e0 + as*as + aw*aw)*gperpon + 1.0e0_vp - gperpon

         CALL OUTPUT(i)

  END DO

    
END SUBROUTINE EXECUTION


SUBROUTINE OUTPUT(i)
  IMPLICIT NONE
  INTEGER i, l, s, aux
  REAL(KIND=vp) NSTEPSTEMP, POINTSTEMP, l2
  REAL(KIND=vp) chi, n, as, aw, theta, ratio, auxratio, ratio2, auxratio2, venergy, wenergy,qenergy, gperp,gz
  REAL(KIND=vp) wl, kl, ks, kw, ws, ww, ws2, ww2, FAUX, sigma, sigma2, gperp2, gperp3, gperp4, gz2, wpe, wpe2
  REAL(KIND=vp) ve, betae, wl2, kl2, nb, vf2
  REAL(KIND=vp), DIMENSION(0:POINTS):: chiout, nout, awout, asout, thetaout, gperpout, dlout, ratioout, ratio2out
  REAL(KIND=vp), DIMENSION(0:POINTS):: venergyout, wenergyout, qenergyout
  REAL(KIND=vp), DIMENSION(0:100):: fauxout, wsout
  COMMON/relativisticfactor/gperp
 
!COMMON BLOCKS WITH EXECUTION SUBROUTINE ***********************
 COMMON/dynvariables/n, aw, as, theta
 COMMON/space/chi
 COMMON/beamparameters/gz2, wpe, wpe2, ve, betae

!COMMON BLOCKS WITH DATAOUTPUT SUBROUTINE **********************
 COMMON/output1/nout, asout, awout, thetaout, ratioout, ratio2out, venergyout, wenergyout, qenergyout
 COMMON/output2/chiout
 COMMON/output3/gperpout
 COMMON/output4/dlout

!COMMON BLOCKS WITH INITIALIZATION SUBROUTINE
 COMMON/fields/wl, kw, kl, ks, ww, ws
 COMMON/fauxbl/FAUX

  ww2 = ww*ww
  ws2 = ws*ws
 ! kl2 = kl*kl
!  wl2 = wl*wl
  ws = sqrt(ws2)
  ww = sqrt(ww2)
  kl2 = kl*kl

  IF(ABS(n) > 4.0e-1_vp) THEN
  WRITE (*,*) 'Plasma density evolution is out of the theory validity.  |n| > 1.0.  Execution interrupted.'
  STOP
  END IF   

  venergy = me*c2*aw*as/(2.0e0_vp*gz)
  wenergy = -me*wpe2*n/kl2
  qenergy = -((hbar*hbar*kl2)/(4.0e0_vp*me*gz*gz*gz))*n
  ratio = venergy/wenergy
  ratio2 = qenergy/wenergy
  

  POINTSTEMP= POINTS
  NSTEPSTEMP   = NSTEPS
  
  l2 = i*POINTSTEMP/NSTEPSTEMP

 
!       IF (venergy > wenergy) THEN
!       auxi = 1
!       ELSE
!       auxi = 0
!       END IF

!       auxold = auxi      
!       IF (auxi == 1) THEN
!       WRITE(*,*) 'Compton regime'
!       ELSE
!       WRITE(*,*) 'Raman regime'
!       END IF
      
       
!  END IF
    
  IF (l2==int(l2)) THEN

     
     l=int(l2)    
     nout(l) = n
     asout(l) = as
     awout(l) = aw
     thetaout(l) = theta
     chiout(l) = chi
     gperpout(l) = gperp
     !dlout(l) = F(ws, ww, gz, ve, wpe, nb,
    ! dlout(l) = 
     venergyout(l) = venergy
     wenergyout(l) = wenergy
     qenergyout(l) = qenergy
     ratioout(l) = ratio
     ratio2out(l) = ratio2
 
  END IF

 !   WRITE(*,*) 'initial detunning (out) =           ',  F(ws, ww, gamae, ve, wpe, gperp)
 !   WRITE(*,*) 'normalized inicial detunning (out) =', F(ws, ww, gamae, ve, wpe, gperp)/wl
  
END SUBROUTINE OUTPUT


!***********************


SUBROUTINE DATAOUTPUT()
  IMPLICIT NONE
  INTEGER i
  REAL(KIND=vp), DIMENSION(0:POINTS):: chiout, nout, awout, asout, thetaout, gperpout, dlout, ratioout, ratio2out
  REAL(KIND=vp) venergyout, wenergyout, qenergyout
  REAL(KIND=vp) A0, A1, A2, kks
  REAL(KIND=vp) I0S, I0W, ws, ww, wl, ks, kw, kl, kl2, ws2, ww2, lambdas, lambdaw, gz2, gz, wpe, wpe2, ve, betae
  REAL(KIND=vp) asout2old, awout2old,nout2old
!COMMON BLOCKS WITH OUTPUT SUBROUTINE ********************************
!COMMON BLOCKS WITH OUTPUT SUBROUTINE ********************************
 COMMON/output1/nout, asout, awout, thetaout, ratioout, ratio2out, venergyout, wenergyout, qenergyout
 COMMON/beamparameters/gz2, wpe, wpe2, ve, betae
 COMMON/output2/chiout
 COMMON/constants/A0, A1, A2
 COMMON/fields/wl, kw, kl, ws2, ww2, ks
 COMMON/output3/gperpout
 COMMON/output4/dlout
 

 lambdas = 2.0e0_vp*pi/ks
  lambdaw = 2.0e0_vp*pi/kw

  gz = sqrt(gz2)
  kl2 = kl*kl

  I0S = (2.0e12_vp)/((8.6)*(8.6)*(lambdas)*(lambdas))
  I0W = (2.0e12_vp)/((8.6)*(8.6)*(lambdaw)*(lambdaw))
 
    write(*,*) 'I0S= ', I0S
  write(*,*) 'lambdas= ', lambdas

  asout2old = asout(0)*asout(0)
  awout2old = awout(0)*awout(0)
  nout2old  = nout(0)*nout(0)

  DO i=0, POINTS

  IF ( asout(i)*asout(i) > asout2old ) THEN
  asout2old = asout(i)*asout(i)
  END IF
  IF ( awout(i)*awout(i) < awout2old ) THEN
  awout2old = awout(i)*awout(i)
  END IF
  IF ( nout(i)*nout(i) > nout2old ) THEN
  nout2old = nout(i)*nout(i)
  END IF



   !write(*,*) 'ios = ', I0S
   !write(*,*) 'iow = ', I0W
  !chiout(i) = chiout(i)*c*epsonchi/wpe
   
  WRITE(18,100) chiout(i), abs(asout(i))/asscale, abs(awout(i))/awscale, nout(i)/nscale1, thetaout(i)
  WRITE(19,100) chiout(i), asout(i)*asout(i), awout(i)*awout(i), abs(nout(i)), thetaout(i)
  WRITE(20,100) chiout(i), I0S*asout(i)*asout(i)/as2scalei, I0W*awout(i)*awout(i)/aw2scalei, abs(nout(i)/nscalei), thetaout(i)
  WRITE(24,101) chiout(i), (awout(0)*awout(0) - awout(i)*awout(i))/(asout(i)*asout(i)) 
!  WRITE(25,102) chiout(i), abs(venergyout(i)/(me*c2))*1.0e12_vp, abs(wenergyout(i)/(me*c2))*1.0e9_vp, abs(qenergyout(i)/(me*c2))*1.0e9_vp, abs(venergyout(i)/wenergyout(i)), log(abs(venergyout(i)/wenergyout(i)))
 
  WRITE(22,100) chiout(i), gperpout(i)-gperpout(0), dlout(i)

 
  END DO
  WRITE(*,*) 'Final position   =', chiout(POINTS)
  WRITE(*,*) 'Final as         =', asout(POINTS)
  WRITE(*,*) 'Final aw         =', awout(POINTS)
  WRITE(*,*) 'Final n          =', nout(POINTS)
  WRITE(*,*) 'Final theta      =', thetaout(POINTS)
  WRITE(*,*) 'Final norm detunn=', dlout(POINTS)/wl 
  WRITE(*,*) wpe2, kl*kl
  WRITE(*,*) 'max as2 plot =', asout2old
  WRITE(*,*) 'max aw2 plot =', awout2old
  WRITE(*,*) 'max n2  plot =', nout2old
  
100 FORMAT(E12.7,1X,E17.11,1X,E30.24,1X,E17.11,1X,E18.12)
101 FORMAT(E12.7,4X,E34.28)
102 FORMAT(E12.7,2X,E18.12,2X,E18.12,2X,E18.12,2X,E18.12,2X,E18.12)
END SUBROUTINE DATAOUTPUT



END PROGRAM  Fortran



