PROGRAM GGsimul 


implicit none


! Project: Gender Gaps 
! Author: Clement Joubert
! Date: December 14 2015

! This code simulates salary distributions at the world bank for the purpose of decomposing pay gaps as des-
! cribed in :
! 	"gender-pay-gaps.pdf"

! Tasks:
! 	(1) Set program specifications
!	(2) Import salary distributions
!	(3) Loop over years and types simulating salary distributions
!	(4) Output simulations

! Specification variables 
character*200,parameter:: path=""
character*200:: pathinput
character*200:: pathoutput
integer,parameter:: tenure      = 0   				! toggles the Career Gap decomposition
integer,parameter:: yp          = 0				! treats yp as a separate category
integer,parameter:: year1       = 1987
integer,parameter:: yearT       = 2027 
integer,parameter:: lastyeardata= 2015 
integer,parameter:: nyrs        = yearT-year1+1 
integer,parameter:: nyrsdata    = lastyeardata-year1+1 
integer,parameter:: ntypes      = 24				! number of employee groups
integer,parameter:: nbins       = 500				! number of salary bins
real(8),parameter:: binsize     = 1.0
integer,parameter:: nclones     = 10
real(8),parameter:: seed        = 145677.0d-0
integer,parameter:: regdum      = 1

! Import_salary variables
integer:: nend(nyrs,ntypes), nentry(nyrs,ntypes), nexit(nyrs,ntypes),ninter(nyrs,ntypes),tempint,gap,gaptype(ntypes)
real(8):: muend(nyrs,ntypes),     sigend(nyrs,ntypes),&
        & muentry(nyrs,ntypes), sigentry(nyrs,ntypes),&
        & muexit(nyrs,ntypes),   sigexit(nyrs,ntypes), &
        & muint(nyrs,ntypes), sigint(nyrs,ntypes),tempr
real(8):: alphagrowth(nyrs,ntypes), &
         & betagrowth(nyrs,ntypes), &
         &  siggrowth(nyrs,ntypes), &
         & noincrease(nyrs,ntypes), &
         & entryrate(nyrs,ntypes), &
         & exitrate(nyrs,ntypes), &
         & alphaexit(nyrs,ntypes), betaexit(nyrs,ntypes), beta2exit(nyrs,ntypes)

! Output statistics
integer:: nend_sim(nyrs,ntypes),nint_sim(nyrs,ntypes),nintdif(nyrs,ntypes)
real(8):: muend_sim(nyrs,ntypes),muenddif(nyrs,ntypes),&
        & sigend_sim(nyrs,ntypes), sigenddif(nyrs,ntypes), &
        & muint_sim(nyrs,ntypes),muintdif(nyrs,ntypes),&
        & sigint_sim(nyrs,ntypes), sigintdif(nyrs,ntypes),nenddif(nyrs,ntypes)

! internal variables
integer:: ss,typecor(ntypes),yy,tt,bin,idum,newbin,i,yr,tp,nolegacy,noentrygap,noexitgap,nogrowthgap, nocompgap 

integer::   sal_sim(nyrs,ntypes,nbins),&
          &               sal_sim_int(nyrs,ntypes,nbins), sal_data(nyrs,ntypes,nbins), &
          &                 entry_sim(nyrs,ntypes,nbins),  entry_data(nyrs,ntypes,nbins),&
          &                  exit_sim(nyrs,ntypes,nbins), exit_data(nyrs,ntypes,nbins)
real(8)::median_sal(nbins),sal,raise,newsal,sig,mu,draw(1), temp(nbins),bub,blb
real(8),allocatable:: raisedraw(:),saldraw(:),noincreasedraw(:)
real(8)::probexit, res,test, test2,test3,test4,test5
character*6::specification

! misc.

    median_sal=(/(i,i=1,nbins)/)*binsize+binsize/2
    if (tenure.eq.1) then
        pathinput=adjustl(trim(path))//"./Estimation/"//"CareerGap/"
        pathoutput=adjustl(trim(path))//"Simulations/"//"CareerGap/"
    elseif (tenure.eq.0) then
        pathinput=adjustl(trim(path))//"./Estimation/"//"AggregateGap/"
        pathoutput=adjustl(trim(path))//"Simulations/"//"AggregateGap/"
    endif
    

! TASK 1 - Set program specifications
    do gap=1,2

    ! type correspondence for counterfactual (i.e. simulate type x with parameters from type y)
    
    if (gap.eq.1) then ! gender gap
            typecor=(/(0,i=1,(ntypes/2)),(i,i=1,(ntypes/2))/)
    else if (gap.eq.2) then ! country part gap
            typecor=(/(0,i=1,(ntypes/4)),(i,i=1,(ntypes/4)),(0,i=1,(ntypes/4)),(i+ntypes/2,i=1,(ntypes/4))/)
    endif
    
    do noexitgap=0,1
    do nogrowthgap=0,1
    do noentrygap=0,1
    do nocompgap=0,1
    do nolegacy=0,1
    
    idum = dint(-1*abs(seed))
    
    if (noentrygap+noexitgap+nogrowthgap+nolegacy+nocompgap.gt.0) then 
        cycle !skip counterfactuals
    endif

    if (noexitgap.lt.nogrowthgap.or.nogrowthgap.lt.noentrygap.or.noentrygap.lt.nocompgap.or.nocompgap.lt.nolegacy) then 
        cycle !do only preferred order
    endif

    print*, noexitgap,nogrowthgap, noentrygap, nocompgap, nolegacy
    write(fmt='(6I1)',unit=specification) gap,nolegacy, nocompgap, noentrygap, nogrowthgap, noexitgap


! TASK 2 - Import data salary distributions
!   note: eventually/alternatively may import whole data set of salaries
!
!   nend, muend, sigend:                   distrib. of salaries at end of period
!   nentry, muentry, sigentry:             distrib. of salaries among entrants
!   nexit, muexit, sigexit:                distrib. of salaries among attriters
!   alphagrowth, betagrowth, siggrowth:    distrib. of logsalary growth
!   alphaexit, betaexit, beta2exit:        coefficients in probability of exit model
   
    call import_salary  (nbins,year1,nyrs,nyrsdata, ntypes, &
                    &    ninter, muint, sigint, &         !   end = end of period
                    &    nend, muend, sigend, &         !   end = end of period
                    &    nentry, muentry, sigentry, &
                    &    nexit, muexit, sigexit, &
                    &    alphagrowth, betagrowth, siggrowth, &
                    &    alphaexit, betaexit, beta2exit, &
                    &    noincrease,pathinput,entry_data,exit_data,sal_data)
    close(1)
    
  
    ! Simulate "nclones" copies of each employee
    sal_data=sal_data*nclones
    entry_data=entry_data*nclones
    exit_data=exit_data*nclones
    nend=nend*nclones
    ninter=ninter*nclones
    nentry=nentry*nclones
    nexit=nexit*nclones
    ! nentry=sum(entry_data,3)
    ! nexit=sum(exit_data,3)
    
    ! Change units to 000s
!    muentry=muentry/1000
!    sigentry=sigentry/1000
!    muexit=muexit/1000
!    sigexit=sigexit/1000
!    muend=muend/1000
!    sigend=sigend/1000
!    muint=muint/1000
!    sigint=sigint/1000
    
!    do yy=1,nyrs
!    do tt=1,ntypes
!        muentry(yy,tt)=dot_product(entry_data(yy,tt,:),median_sal(:))
!        muentry(yy,tt)=muentry(yy,tt)/nentry(yy,tt)
!        muexit(yy,tt)=dot_product(exit_data(yy,tt,:),median_sal(:))
!        muexit(yy,tt)=muexit(yy,tt)/nexit(yy,tt)
!    enddo
!    enddo

    ! Compute entry and exit rates
    entryrate=nentry*1.0d-0/ninter
    exitrate=nexit*1.0d-0/ninter


! TASK 3 - Simulate salary distributions one-year-ahead for each type of employee
!
! Simulations are stored in sal_sim(ntypes,nyrs,nbins). 
! Entries correspond to the number of individuals of a given type, 
! with salaries within the range of a given bin in a given year


    ! 0. Initialize year 1 and entries using data
    sal_sim=0
    sal_sim_int=0
    entry_sim(1,:,:)=entry_data(1,:,:)

    ! 0bis. Use counterfactual distributions if applicable
    	
    do tt=1,ntypes
    if (typecor(tt).gt.0) then
        if (nocompgap==1) then
            nentry       (:,tt)=nentry      (:,typecor(tt))            
        endif
        if (noentrygap==1) then 
            muentry      (:,tt)=muentry     (:,typecor(tt))
            sigentry     (:,tt)=sigentry    (:,typecor(tt))
        endif
        if (noexitgap==1) then 
            alphaexit    (:,tt)=alphaexit   (:,typecor(tt))
            betaexit     (:,tt)=betaexit    (:,typecor(tt))
            beta2exit    (:,tt)=beta2exit   (:,typecor(tt))
        endif
        if (nogrowthgap==1) then 
            alphagrowth (:,tt)=alphagrowth  (:,typecor(tt))
            betagrowth  (:,tt)=betagrowth   (:,typecor(tt))
            siggrowth   (:,tt)=siggrowth    (:,typecor(tt))
            noincrease  (:,tt)=noincrease   (:,typecor(tt))
        endif
        if (nolegacy==1) then 
            entry_sim (1,tt,:)=entry_sim  (1,typecor(tt),:)
        endif
    endif
    enddo

    ! 1. Simulate salaries one year ahead
    ! 
    ! Use data for number of entries and exits as well as distrib. of salaries for
    ! entrants, exiters and distrib. of salary growth

    do yy=1,nyrs       
    yr=yy+year1-1       !yr denotes calendar year, while yy is a counter
    
    ! If yy is larger than last year of data, use parameters from last year of data

    if (yr.gt.lastyeardata) then
	nentry(yy,:)=nentry(yy-1,:)
	muentry(yy,:)=muentry(yy-1,:)
	sigentry(yy,:)=sigentry(yy-1,:)
	alphaexit(yy,:)=alphaexit(yy-1,:)
	betaexit(yy,:)=betaexit(yy-1,:)
	beta2exit(yy,:)=beta2exit(yy-1,:)
	alphagrowth(yy,:)=alphagrowth(yy-1,:)
	betagrowth(yy,:)=betagrowth(yy-1,:)
	siggrowth(yy,:)=siggrowth(yy-1,:)
	noincrease(yy,:)=noincrease(yy-1,:)
    endif
    

    do tt=1,ntypes

    ! a. Salary growth
    
        if (yy>1)    call salgrowth(nyrs,ntypes,binsize,nbins,sal_sim,sal_sim_int,idum,noincrease,alphagrowth,betagrowth,siggrowth,tt,yy)
       
    ! b. Entries and exits
 
        res=0.0d-0            
        do bin=1,nbins

            ! Add new entrants
            if (yy>1.or.tenure==1) then   !For 1987, use exact salary bins, for later years, simulate based on estimated average and sd 
                    blb= 0.0d-0 + (bin-1)*binsize
                    if (bin.eq.1) blb=-999999999999.0d-0
                    bub= bin*binsize
                if (sigentry(yy,tt).gt.0.0d-0)  then 
                    if (bin.eq.nbins) bub=999999999999.0d-0
                    call hist(entry_sim(yy,tt,bin),nentry(yy,tt),muentry(yy,tt),sigentry(yy,tt),blb,bub,res)
                elseif (sigentry(yy,tt).eq.0.0d-0.and.muentry(yy,tt).gt.blb.and.muentry(yy,tt).le.bub) then
                    entry_sim(yy,tt,bin)=nentry(yy,tt) !in case the distribution of salary entrants is degenerate
                else
                    entry_sim(yy,tt,bin)=0
                endif
            endif
            
            sal_sim(yy,tt,bin)=sal_sim_int(yy,tt,bin)+entry_sim(yy,tt,bin)
            
            ! For each salary bin, simulate how many employees leave
            tempr=probexit(median_sal(bin),alphaexit(yy,tt),betaexit(yy,tt),beta2exit(yy,tt),regdum)
            tempr=tempr*sal_sim(yy,tt,bin)
            tempint=idnint(tempr)
            exit_sim(yy,tt,bin)=tempint
            
            
            ! Combine salaries of entrants exiters and employee stock 
            sal_sim(yy,tt,bin)=sal_sim(yy,tt,bin)-exit_sim(yy,tt,bin)
        
        enddo
    enddo
    enddo
    
! TASK 4 - Compare actual and simulated end-of-year distributions
    
    call fit(pathoutput,sal_sim, sal_data, sal_sim_int,entry_sim,entry_data,exit_sim,exit_data,nend,ninter,nexit,nentry,binsize,nbins,nyrs, ntypes, muend,sigend,muint,sigint,muentry,muexit,specification,tenure)
    
enddo
enddo
enddo
enddo
enddo
enddo


print*,"the end"

END PROGRAM GGsimul
	
subroutine hist(freq,n,mu,sigma,blb,bub,res)

external alnorm

!*****************************************************************************80
!  Modified:
!
!    17 December 2015
!
!  Author:
!
!    Clement Joubert
!
!  Parameters:
!   Inputs
    integer:: n !is the size of the population
    real(8):: mu, sigma !are the parameters of the normal distribution
    real(8):: res !residual frequency from rounding the frequency of the previous bin
!
!    Output
    integer:: freq !frequency of the bin
!****************************************************************************80
    real(8):: binproba,blb,bub,sblb, sbub, alnorm, lcdf, ucdf

    
    sblb=(blb-mu)/sigma   !standardized bin lower bound
    sbub=(bub-mu)/sigma   !standardized bin upper bound
    
    lcdf=alnorm(sblb,.false.)   !cdf at lower bound
    ucdf=alnorm(sbub,.false.)   !cdf at upper bound
    
    binproba= ucdf-lcdf 
    
    if (blb.gt.18) then
        freq = nint(n*binproba + res)  
    else
        freq = 0
    endif
    res = n*binproba + res - freq
    
endsubroutine


function probexit(sal,alpha,beta,beta2,regdum)

!*****************************************************************************80
!  Modified:
!
!
!  Author:
!
!    Clement Joubert
!
!  Parameters:
!   Inputs
    real(8):: temp, alpha, beta, beta2, sal
    integer:: regdum !if using a linear probability model
!    Output

    real(8):: probexit 
!****************************************************************************80   

if (regdum.eq.0) then
    temp=exp(-alpha- beta*sal - beta2*sal*sal)  
    probexit = 1/(1+temp)  
elseif (regdum.eq.1) then
    temp=beta*(sal-alpha)+beta2
    probexit= dmin1(1.0d-0,dmax1(0.0d-0,temp))
endif

    
endfunction


