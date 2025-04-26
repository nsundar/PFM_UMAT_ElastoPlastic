!************************************************************************
! This UMAT is adapted and extended from the work of: 
! 
! Yousef Navidtehrani, Covadonga Betegón, and Emilio Martínez-Pañeda.
! "A unified Abaqus implementation of the phase field fracture method 
! using only a user material subroutine."
! Materials, 14(8), 2021. 
! ISSN: 1996-1944.
! DOI: 10.3390/ma14081913
! URL: https://www.mdpi.com/1996-1944/14/8/1913
!************************************************************************

          subroutine umat(stress,statev,ddsdde,sse,spd,scd,rpl,ddsddt, &
      drplde,drpldt,stran,dstran,time,dtime,temp,dtemp,predef,dpred, &
      cmname,ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt, &
      celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,jstep,kinc)

      include 'aba_param.inc'

      character*8 cmname
      dimension stress(ntens),statev(nstatv),ddsdde(ntens,ntens), &
      ddsddt(ntens),drplde(ntens),stran(ntens),dstran(ntens),time(2), &
      predef(1),dpred(1),props(nprops),coords(3),drot(3,3),dfgrd0(3,3), &
      dfgrd1(3,3),jstep(4), eplas(ntens), eelas(ntens),PS(3)
               
      real*8, parameter :: pi = 3.1415926535897932384626433832d0
      real*8, parameter :: newton = 100
      real*8 :: stressi(ntens), flow(ntens)

      ddsdde=0.d0
      Hmin=0.d0
      E=props(1) ! Young's modulus
      xnu=props(2) ! Poisson's ratio
      xl=props(3) ! Phase field length scale parameter
      Gc=props(4) !Fracture Toughness
      xk=1.d-07 !well conditioning parameter
     
      phi=temp+dtemp     !phase field parameter 
      
!     Degradation function
      call KDegFun(phi,xk,g,dg,ddg)     

	eg=E/(1.d0+xnu)/2.d0
      elam=(E/(1.d0-2.d0*xnu)-2.d0*eg)/3.d0
      bk=E/(1.d0-2.d0*xnu)/3.d0
      g=(1.d0-phi)**2+xk
      eg3=eg*3.d0
      
        do i=1,3
         do j=1,3
          ddsdde(j,i)=elam
         end do
         ddsdde(i,i)=eg*2.d0+elam
        end do
        do i=4,ntens
         ddsdde(i,i)=eg
        end do
	call ROTSIG(statev(1),drot,eelas,2,ndi,nshr)
	call ROTSIG(statev(ntens+1),drot,eplas,2,ndi,nshr)
	
	!Recover undegraded stress   
	do k1=1,ntens
          stressi(k1)= statev(k1+2*ntens)        
        enddo
        
      eqplas=statev(1+3*ntens)  ! equivalent plastic  strain
      spi = statev(3+3*ntens)   ! specific plastic dissipation
      psit = statev(2+3*ntens)  ! elastic strain energy
 !calculate predictor stress and elastic strain
      do k1=1,ntens
         eelas(k1)=eelas(k1)+dstran(k1)           
       do k2=1,ntens
            stressi(k2) = stressi(k2) + ddsdde(k2,k1) * dstran(k1)
        enddo                
      enddo
      
       !     calculate mises stress       
        smises=(stressi(1)-stressi(2))**2 &
            +(stressi(2)-stressi(3))**2 &
            +(stressi(3)-stressi(1))**2    
!                
        do k1=ndi+1,ntens
           smises=smises+(6.d0 *stressi(k1)**2)
        enddo 
	   smises=sqrt(smises/2.d0)
!
!     hardening curve to get yield stress

       call hardsub (syield0,hard,eqplas,props)
       
!  Determine if actively yielding
	
	if (smises.gt. (1.d0+xk)*syield0) then        
          shydro=(stressi(1)+ stressi(2)+ stressi(3))/(3.d0) !hydrostatic stress   
          do k1=1,ndi
              flow(k1)=(stressi(k1)-shydro)/smises   !flow rule
          enddo
          do k1=ndi+1,ntens
              flow(k1)= stressi(k1)/smises
          enddo

       !     solve for equivalent stress, newton iteration
          syield=syield0
          deqpl=0.d0
          do kewton=1,newton
             rhs=smises-eg3*deqpl-syield
             deqpl=deqpl+rhs/(eg3+hard)
      
             call hardsub(syield,hard,eqplas+deqpl,props)
             if(abs(rhs).lt.xk*syield0) goto 110
          enddo
        write(6,109) int(newton)
 109      format(//, 'plasticity algorithm did not converge after ',i3,' iterations')
 110      continue
    
!     calculate stress and update strains
          do k1=1,ndi
            stressi(k1)= flow(k1)*syield+shydro
            eplas(k1)= eplas(k1)+3.d0* flow(k1)*deqpl/2.d0
            eelas(k1)= eelas(k1)-3.d0* flow(k1)*deqpl/2.d0
          enddo
          do k1=ndi+1,ntens
              stressi(k1)= flow(k1)*syield
              eplas(k1)= eplas(k1)+3.d0* flow(k1)*deqpl
              eelas(k1)= eelas(k1)-3.d0* flow(k1)*deqpl
          enddo

          eqplas=eqplas+deqpl          
          spd=spi+deqpl*(syield0+syield)/2.d0! update specific plastic dissipation         
       
!     material jacobian
	  ebulk3=bk*3.d0
          effg=eg*syield/smises
          effg2=2.d0*effg
          effg3=(3.d0*effg2)/2.d0
          efflam=(ebulk3-effg2)/3.d0
          effhrd=eg3*hard/(eg3+hard)-effg3
!          
          do k1=1,ndi
             do k2=1,ndi
                ddsdde(k2,k1)=efflam
             enddo
             ddsdde(k1,k1)=effg2+efflam
          enddo 
!          
          do k1=ndi+1,ntens
             ddsdde(k1,k1)=effg
          enddo   

          do k1=1,ntens
             do  k2=1,ntens
                ddsdde(k2,k1)=ddsdde(k2,k1)+ flow(k2) &
                           * flow(k1)*effhrd
              enddo
          enddo
        endif
    
!      
      stress= stressi*g  ! degraded stress
      ddsdde= ddsdde*g    ! degraded material stiffness matrix
!
   call SPRINC(stran,PS,2,ndi,nshr)
        trp1=(PS(1)+PS(2)+PS(3)+abs(PS(1)+PS(2)+PS(3)))/2.d0
        trn1=(PS(1)+PS(2)+PS(3)-abs(PS(1)+PS(2)+PS(3)))/2.d0
        trp2=0.d0
        trn2=0.d0
        do i=1,3
         trp2=trp2+(PS(i)+abs(PS(i)))**2.d0/4.d0
         trn2=trn2+(PS(i)-abs(PS(i)))**2.d0/4.d0
        end do
        psip=xnu*eg/(1d0-2d0*xnu)*trp1**2d0+eg*trp2
        psin=xnu*eg/(1d0-2d0*xnu)*trn1**2d0+eg*trn2

         He=max(psip,psit,Hmin)      ! Elastic history field

      evol = stran(1)+stran(2)+stran(3)
            if (evol .ge. 0d0) then
      Hp=spd    !Plastic history field
      else
      Hp=spi
      end if
           
        do k1=1,ntens
        statev(k1)= eelas(k1)
        statev(k1+ntens)= eplas(k1)
        statev(k1+2*ntens)= stressi(k1)
      enddo 
        statev(1+3*ntens)=eqplas
        statev(2+3*ntens)=He
        statev(3+3*ntens)=Hp   
        statev(4+3*ntens)=phi
   !     phase field
       w=phi**2
       dw=2.d0*phi
       ddw=2.d0
       cw=0.5d0             
  !     heat transfer analogy
         rpl=-(dg*(He+Hp)*2.d0*cw/(xl*Gc)+dw/(2.d0*xl**2))
      drpldt=-(ddg*(He+Hp)*2.d0*cw/(xl*Gc)+ddw/(2.d0*xl**2))
             
      RETURN 
      end
!***********************************************************************
!     Degradation function
      subroutine KDegFun(phi,xk,g,dg,ddg)
      
      include 'aba_param.inc'
      
      !AT2 & AT1 model
       g=(1.d0-phi)**2+xk
       dg=-2.d0*(1.d0-phi)
       ddg=2.d0                

      end
!***********************************************************************
	subroutine hardsub(syield,hard,eqplas,props)
	
	include 'aba_param.inc'
	
	dimension props(nprops)
	!     Current yield stress and hardening    
	 
	 eqpl0=props(5)
	 eqpl1=props(6)
	 deqpl=eqpl1-eqpl0
         syield0=props(7)
         syield1=props(8)
         dsyield=syield1-syield0
         hard=dsyield/deqpl
         syield=syield0+(eqplas-eqpl0)*hard        
         
         end     
	
	
	
