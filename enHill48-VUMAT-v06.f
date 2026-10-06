c====================================================================
c          Program for 3D anisotropic hardening 
c          Anisotropic Hill48 yield function
c          Hill7 - anisotropic hardening & r-value evolution
c          by Junhe Lian
c          junhe.lian@iehk.rwth-aachen.de
c          February 2016 
c          Abaqus 6.11-1
c
c          !DO NOT DISTRIBUTE WITHOUT AUTHOR'S PERMISSION!
c====================================================================
c
       subroutine vumat(
C Read only -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     3  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     5  stressNew, stateNew, enerInternNew, enerInelasNew)
C
      include 'vaba_param.inc'
C
C All arrays dimensioned by (*) are not used in this algorithm
      dimension props(nprops), density(nblock),
     1  coordMp(nblock,*),
     2  charLength(*), strainInc(nblock,ndir+nshr),
     3  relSpinInc(*), tempOld(*),
     4  stretchOld(*), defgradOld(*),
     5  fieldOld(*), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(*),
     8  stretchNew(*), defgradNew(*), fieldNew(*),
     9  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     1  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C
C Defining numerical constants
      parameter(zero=0.,one=1.,two=2.,three=3.,four=4.,six=6.,
     1  half=0.5,threeHalves = 1.5,third = one/three,Tol=1.D-8,
     1  twoThirds=two/three,const1=sqrt(threeHalves),NTENS=6,NDI=3)
      REAL A0,e0,xn0,A45,e45,xn45,A90,e90,xn90,AEB,eEB,xnEB,r0,r45,r90,
     1  rEB,F,G,H,L,M,N,F1,G1,H1,L1,M1,N1	 
C     User defined tensors
      DIMENSION dfds(6),Ydirection(6),dgds(6),DSTRESS(6),
     1 ST(6),EELAS(6),EPLAS(6),EELT(6),EPLT(6),STRESST(6),
     2 EELAST(6),EPLAST(6),DDSDDE(6,6),DSTRAN(6) 
C Define material properties constants
c ELASTIC properties
       EMOD = PROPS(1)
       ENU  = PROPS(2)
c ELASTIC STIFFNESS
       EBULK3=EMOD/(one-two*ENU)
       EG2=EMOD/(one+ENU)
       EG=EG2/two
       EG3=three*EG
       ELAM=(EBULK3-EG2)/three
	   PI= four*atan(one)
c Define dimension of flow curve and r-value input table	   
       nvalue = (nprops-8)/8	   
C----------------------------------------------------------------------C
C Fill in the elastic stiffness tensor
      DO 20 K1=1,NTENS
        DO 10 K2=1,NTENS
          DDSDDE(K2,K1)=zero
 10     CONTINUE
 20   CONTINUE
C
      DO 40 K1=1,NDI
        DO 30 K2=1,NDI
          DDSDDE(K2,K1)=ELAM
 30     CONTINUE
        DDSDDE(K1,K1)=EG2+ELAM
 40   CONTINUE
      DO 50 K1=NDI+1,NTENS
        DDSDDE(K1,K1)=EG
 50   CONTINUE
C The stiffness is material-constant; build it once per VUMAT call.
C Computation per material point starts here
      do 9000 i = 1,nblock
C  Elastic predictor and plastic corrector scheme in conjuction with the
C  Tangent Cutting Plane (TCP) algorithm 
C  Elastic predictor
      eqplastrue = stateOld(i,1)
      eqplas = eqplastrue
C Transite DSTRAN from strainInc, note diff. between UMAT and VUMAT.  
      DSTRAN(1)=strainInc(i,1)
      DSTRAN(2)=strainInc(i,2)
      DSTRAN(3)=strainInc(i,3)
      DSTRAN(4)=strainInc(i,4)*two
      DSTRAN(5)=strainInc(i,5)*two
      DSTRAN(6)=strainInc(i,6)*two
C Calculate Trial stress
      DO 70 K1=1,NTENS
          STRESST(K1)=stressOld(i,K1)+DDSDDE(K1,1)*DSTRAN(1)
     1   +DDSDDE(K1,2)*DSTRAN(2)+DDSDDE(K1,3)*DSTRAN(3)
     1   +DDSDDE(K1,4)*DSTRAN(4)+DDSDDE(K1,5)*DSTRAN(5)
     1   +DDSDDE(K1,6)*DSTRAN(6)
          ST(K1)=STRESST(K1)
 70   CONTINUE
       ST1=STRESST(1)
       ST2=STRESST(2)
       ST3=STRESST(3)
       ST4=STRESST(4)
       ST5=STRESST(5)
       ST6=STRESST(6)  
C Calculate anisotropic parameters in Hill48 based on flow stresses	          
C Fetching yield stress from flow curve
       call ahard(sigmayield0,hard0,sigmayield45,sigmayield90,sigmayieldEB,
     1	  rvalue0,rvalue45,rvalue90,eqplas,PROPS(9),nvalue)	   
       sigmayield=sigmayield0
	   xhard=hard0
	   s0=sigmayield
       s45=sigmayield45
       s90=sigmayield90
       sEB=sigmayieldEB
       RS90=(s0/s90)**2
       RSEB=(s0/sEB)**2
       RS45=(s0/s45)**2
       F=(RS90-one+RSEB)/two
       G=(one-RS90+RSEB)/two
       H=(one+RS90-RSEB)/two
       N=(four*RS45-RSEB)/two
	   L=threeHalves
	   M=threeHalves
C Calculate subvalue of the trial Hill48 equivalent stress		   
       subsigtrial= F*(ST2-ST3)**2+G*(ST3-ST1)**2
     1             +H*(ST1-ST2)**2+two*N*ST4**2
     2             +two*L*ST5**2+two*M*ST6**2
c To prevent subsigtrial =0
       if(subsigtrial.lt.Tol) subsigtrial = Tol
C Calculate the trial Hill48 equivalent stress 
       sigtrial = sqrt(subsigtrial)
c To prevent sigtrial =0	   
       if(sigtrial.lt.Tol .and. sigtrial.ge.zero) sigtrial = Tol
       if(sigtrial.gt.(zero-Tol) .and. sigtrial.le.zero) sigtrial = zero-Tol	   
C Check yielding
       radius = sigmayield
	   Xhardnew= xhard
       xf = sigtrial - radius
       if(xf.le.(Tol)) then            		 
          goto 1000
       endif
       dgama=zero
       iflag=zero
C Plastic corrector and Begin iterations
5000   continue
       iflag=iflag+one   
c Calculate anisotropic parameters in Hill48 based on flow stresses
C The current curve values were already evaluated either by the elastic
C predictor (first iteration) or at the end of the previous iteration.
       sigmayield=sigmayield0
	   xhard=hard0
	   s0=sigmayield
       s45=sigmayield45
       s90=sigmayield90
       sEB=sigmayieldEB
       r0=rvalue0
       r45=rvalue45
       r90=rvalue90
       RS90=(s0/s90)**2
       RSEB=(s0/sEB)**2
       RS45=(s0/s45)**2
       F=(RS90-one+RSEB)/two
       G=(one-RS90+RSEB)/two
       H=(one+RS90-RSEB)/two
       N=(four*RS45-RSEB)/two
	   L=threeHalves
	   M=threeHalves
C Calculate subvalue of the trial Hill48 equivalent stress		   
       subsigtrial=F*(ST2-ST3)**2+G*(ST3-ST1)**2
     1             +H*(ST1-ST2)**2+two*N*ST4**2
     2             +two*L*ST5**2+two*M*ST6**2
c To prevent subsigtrial =0	 
       if(subsigtrial.lt.Tol) subsigtrial = Tol
C Calculate the trial Hill48 equivalent stress	   
       sigtrial = sqrt(subsigtrial)
c To prevent sigtrial =0	   
       if(sigtrial.lt.Tol .and. sigtrial.ge.zero) sigtrial = Tol
       if(sigtrial.gt.(zero-Tol) .and. sigtrial.le.zero) sigtrial = zero-Tol
c Calculate anisotropic parameters in Hill48 based on r-values	 
       F1=r0/(r90)/(one+r0)
       G1=one/(one+r0)
       H1=r0/(one+r0)
       L1=threeHalves
       M1=threeHalves
       N1=(r0+r90)*(one+two*r45)/two/r90/(one+r0)
C Calculate subvalue of the trial Hill48 flow potential
	   subsigmaHill=F1*(ST2-ST3)**2+G1*(ST3-ST1)**2
     1          +H1*(ST1-ST2)**2+two*N1*ST4**2
     2          +two*L1*ST5**2+two*M1*ST6**2
c To prevent subsigmaHill =0	 
       if(subsigmaHill.lt.Tol) subsigmaHill = Tol
C Calculate the trial Hill48 flow potential	   
       sigmaHill=sqrt(subsigmaHill)
c To prevent sigmaHill =0	   
       if(sigmaHill.lt.Tol .and. sigmaHill.ge.zero) sigmaHill = Tol
       if(sigmaHill.gt.(zero-Tol) .and. sigmaHill.le.zero) sigmaHill = zero-Tol		   
C Update parameters with Non-AFR and calculate sigmayield and sigmaHill at every iteration 
C Calculate the flow direction vector based on Flow potential
       dgds(1)=(H1*(ST1-ST2)-G1*(ST3-ST1))/sigmaHill
       dgds(2)=(F1*(ST2-ST3)-H1*(ST1-ST2))/sigmaHill
       dgds(3)=(G1*(ST3-ST1)-F1*(ST2-ST3))/sigmaHill    
       dgds(4)=(two*N1*ST4)/sigmaHill
       dgds(5)=(two*L1*ST5)/sigmaHill
       dgds(6)=(two*M1*ST6)/sigmaHill
C Calculate the flow direction vector based on Yield function 
       dfds(1)=(H*(ST1-ST2)-G*(ST3-ST1))/sigtrial
       dfds(2)=(F*(ST2-ST3)-H*(ST1-ST2))/sigtrial
       dfds(3)=(G*(ST3-ST1)-F*(ST2-ST3))/sigtrial 
       dfds(4)=(two*N*ST4)/sigtrial
       dfds(5)=(two*L*ST5)/sigtrial
       dfds(6)=(two*M*ST6)/sigtrial
C Calculate the flow direction in the TCP algorithm
       DGSUM=dgds(1)+dgds(2)+dgds(3)
       Ydirection(1)=EG2*dgds(1)+ELAM*DGSUM
       Ydirection(2)=EG2*dgds(2)+ELAM*DGSUM
       Ydirection(3)=EG2*dgds(3)+ELAM*DGSUM
       Ydirection(4)=EG*dgds(4)
       Ydirection(5)=EG*dgds(5)
       Ydirection(6)=EG*dgds(6)
C Calculate the first term in the TCP algorithm
c dsdgama caculation=CNN, dot product, prepare for iteration
C Calculate MCN
       dsdga=dfds(1)*Ydirection(1)+dfds(2)*Ydirection(2)
     1      +dfds(3)*Ydirection(3)+dfds(4)*Ydirection(4)
     1      +dfds(5)*Ydirection(5)+dfds(6)*Ydirection(6)
C Calculate the second term in the TCP algorithm
C Correspongding flow direction, Nf
	   radius = sigmayield
	   Xhardnew= xhard
	   Xdeltagamma = dsdga+Xhardnew	   
c To prevent Xdeltagamma =0	   
       if(Xdeltagamma.lt.Tol .and. Xdeltagamma.ge.zero) Xdeltagamma = Tol
       if(Xdeltagamma.gt.(zero-Tol) .and. Xdeltagamma.le.zero) Xdeltagamma=zero-Tol	   
	   xddgamma = xf/Xdeltagamma
       dgama = dgama+xddgamma
c Exact algebraic form: sigma = elastic trial - gamma*C:dgds.
      DO 160 K1=1,NTENS
         STRESST(K1)=ST(K1)-dgama*Ydirection(K1)
160   CONTINUE
c Update stress tensors
      ST1=STRESST(1)
      ST2=STRESST(2)
      ST3=STRESST(3)
      ST4=STRESST(4)
      ST5=STRESST(5)
      ST6=STRESST(6)		  
C Update the equivalent plastic strain	  
      eqplastrue = stateOld(i,1)
      eqplas = eqplastrue+dgama*(sigmaHill/sigtrial)
	  deqplas = dgama*(sigmaHill/sigtrial)	  
C Update anisotropic parameters in Hill48 based on flow stresses
C Fetching yield stress from flow curve
       call ahard(sigmayield0,hard0,sigmayield45,sigmayield90,sigmayieldEB,
     1	  rvalue0,rvalue45,rvalue90,eqplas,PROPS(9),nvalue)	   
       sigmayield=sigmayield0
	   xhard=hard0
	   radius=sigmayield	   
	   s0=sigmayield
       s45=sigmayield45
       s90=sigmayield90
       sEB=sigmayieldEB
       RS90=(s0/s90)**2
       RSEB=(s0/sEB)**2
       RS45=(s0/s45)**2
       F=(RS90-one+RSEB)/two
       G=(one-RS90+RSEB)/two
       H=(one+RS90-RSEB)/two
       N=(four*RS45-RSEB)/two
	   L=threeHalves
	   M=threeHalves	   
C Update subvalue of the trial Hill48 equivalent stress	   
       subsigtrial=F*(ST2-ST3)**2+G*(ST3-ST1)**2
     1            +H*(ST1-ST2)**2+two*N*ST4**2
     2            +two*L*ST5**2+two*M*ST6**2
c To prevent subsigtrial =0	 
       if(subsigtrial.lt.Tol) subsigtrial = Tol
C Update the trial Hill48 equivalent stress 	   
       sigtrial = sqrt(subsigtrial)
c To prevent sigtrial =0	   
       if(sigtrial.lt.Tol .and. sigtrial.ge.zero) sigtrial = Tol
       if(sigtrial.gt.(zero-Tol) .and. sigtrial.le.zero) sigtrial = zero-Tol
C Update the radius	   
	   radius = sigmayield
	   Xhardnew= xhard	   
       xf = sigtrial - radius
C Check if the stress state is located on the updated yield loci
      if(iflag.gt.(10.+Tol)) then
         goto 1000
      endif
      if(xf.le.Tol .and. xf.ge.(zero-Tol)) then
         goto 1000
      else 
         goto 5000
      endif		   

1000  continue
C     Update PEEQ		 
      stateNew(i,1) = eqplas
c      write(*,*) "stateNew(i,1)", stateNew(i,1)	 
C     Update equivalent stress
	  stateNew(i,2) = sigtrial
C     Update stress
      DO 180 K1=1,NTENS          
         stressNew(i,K1)=STRESST(K1)
180   CONTINUE
C   
9000  continue
      return
      end
c
      subroutine ahard(sigmayield0,hard0,sigmayield45,sigmayield90,sigmayieldEB,
     1	  rvalue0,rvalue45,rvalue90,eqplas,table,nvalue)
C
      include 'vaba_param.inc'
      dimension table(8,nvalue)
C
C     Set yield stress to second value of table, hardening to zero
      sigmayield0=table(2,nvalue)
      hard0=zero
	  sigmayield45=table(3,nvalue)
	  hard45=zero
	  sigmayield90=table(4,nvalue)
	  hard90=zero
	  sigmayieldEB=table(5,nvalue)
	  hardEB=zero
	  rvalue0=table(6,nvalue)
	  hardr0=zero
	  rvalue45=table(7,nvalue)
	  hardr45=zero
	  rvalue90=table(8,nvalue)
	  hardr90=zero
C
C     If more than one entry, search table
      if(nvalue.gt.1) then
        do 10 k1=1,nvalue-1
          eqpl1=table(1,k1+1)
          if(eqplas.lt.eqpl1) then
            eqpl0=table(1,k1)
            if(eqpl1.le.eqpl0) then
              write(6,7)
 7            format(//,30X,'***ERROR - PLASTIC STRAIN MUST BE ',
     1               'ENTERED IN ASCENDING ORDER,')
C
C             Subroutine XIT terminates execution and closes all files
              call xplb_exit
            endif
            deqpl=eqpl1-eqpl0
            sigmayield0a=table(2,k1)
            sigmayield0b=table(2,k1+1)
            dsigmayield0=sigmayield0b-sigmayield0a
            hard0=dsigmayield0/deqpl
            sigmayield0=sigmayield0a+(eqplas-eqpl0)*hard0
            sigmayield45a=table(3,k1)
            sigmayield45b=table(3,k1+1)
            dsigmayield45=sigmayield45b-sigmayield45a
            hard45=dsigmayield45/deqpl
            sigmayield45=sigmayield45a+(eqplas-eqpl0)*hard45
            sigmayield90a=table(4,k1)
            sigmayield90b=table(4,k1+1)
            dsigmayield90=sigmayield90b-sigmayield90a
            hard90=dsigmayield90/deqpl
            sigmayield90=sigmayield90a+(eqplas-eqpl0)*hard90
            sigmayieldEBa=table(5,k1)
            sigmayieldEBb=table(5,k1+1)
            dsigmayieldEB=sigmayieldEBb-sigmayieldEBa
            hardEB=dsigmayieldEB/deqpl
            sigmayieldEB=sigmayieldEBa+(eqplas-eqpl0)*hardEB
            rvalue0a=table(6,k1)
            rvalue0b=table(6,k1+1)
            drvalue0=rvalue0b-rvalue0a
            hardr0=drvalue0/deqpl
            rvalue0=rvalue0a+(eqplas-eqpl0)*hardr0
            rvalue45a=table(7,k1)
            rvalue45b=table(7,k1+1)
            drvalue45=rvalue45b-rvalue45a
            hardr45=drvalue45/deqpl
            rvalue45=rvalue45a+(eqplas-eqpl0)*hardr45
            rvalue90a=table(8,k1)
            rvalue90b=table(8,k1+1)
            drvalue90=rvalue90b-rvalue90a
            hardr90=drvalue90/deqpl
            rvalue90=rvalue90a+(eqplas-eqpl0)*hardr90			
            goto 20
          endif
 10     continue
 20     continue
        if(eqplas.gt.table(1,nvalue)) then
          hard0=(table(2,nvalue)-table(2,nvalue-1))
     1        /(table(1,nvalue)-table(1,nvalue-1))
          sigmayield0=table(2,nvalue)+(eqplas-table(1,nvalue))*hard0
          hard45=(table(3,nvalue)-table(3,nvalue-1))
     1        /(table(1,nvalue)-table(1,nvalue-1))
          sigmayield45=table(3,nvalue)+(eqplas-table(1,nvalue))*hard45
          hard90=(table(4,nvalue)-table(4,nvalue-1))
     1        /(table(1,nvalue)-table(1,nvalue-1))
          sigmayield90=table(4,nvalue)+(eqplas-table(1,nvalue))*hard90
          hardEB=(table(5,nvalue)-table(5,nvalue-1))
     1        /(table(1,nvalue)-table(1,nvalue-1))
          sigmayieldEB=table(5,nvalue)+(eqplas-table(1,nvalue))*hardEB
          hardr0=(table(6,nvalue)-table(6,nvalue-1))
     1        /(table(1,nvalue)-table(1,nvalue-1))
          rvalue0=table(6,nvalue)+(eqplas-table(1,nvalue))*hardr0
          hardr45=(table(7,nvalue)-table(7,nvalue-1))
     1        /(table(1,nvalue)-table(1,nvalue-1))
          rvalue45=table(7,nvalue)+(eqplas-table(1,nvalue))*hardr45
          hardr90=(table(8,nvalue)-table(8,nvalue-1))
     1        /(table(1,nvalue)-table(1,nvalue-1))
          rvalue90=table(8,nvalue)+(eqplas-table(1,nvalue))*hardr90	  
        endif
      endif
      return
C
C Iteration ends here
      end

	  
