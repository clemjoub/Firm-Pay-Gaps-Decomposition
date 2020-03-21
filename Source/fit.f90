subroutine fit(pathoutput,sal_sim,sal_data, sal_sim_int,entry_sim,entry_data,exit_sim,exit_data,nend,ninter,nexit,nentry,binsize,nbins,nyrs, ntypes, muend,sigend,muint,sigint,muentry,muexit,specification,tenure)
implicit none

character*200:: pathoutput 
integer:: nyrs 
integer:: ntypes 
integer:: nbins
real(8):: binsize
character*6:: specification
  
   ! Import_salary variables
integer:: nend(nyrs,ntypes), nentry(nyrs,ntypes), nexit(nyrs,ntypes),ninter(nyrs,ntypes), tenure
real(8):: muend(nyrs,ntypes),     sigend(nyrs,ntypes),&
        & muentry(nyrs,ntypes), sigentry(nyrs,ntypes),&
        & muexit(nyrs,ntypes),   sigexit(nyrs,ntypes), &
        & muint(nyrs,ntypes), sigint(nyrs,ntypes)
real(8):: alphagrowth(nyrs,ntypes), &
         & betagrowth(nyrs,ntypes), &
         &  siggrowth(nyrs,ntypes), &
         & noincrease(nyrs,ntypes)

! Output statistics
integer:: nend_sim(nyrs,ntypes),nenddif(nyrs,ntypes),nint_sim(nyrs,ntypes),nintdif(nyrs,ntypes),nentry_sim(nyrs,ntypes),nexit_sim(nyrs,ntypes),nentrydif(nyrs,ntypes),nexitdif(nyrs,ntypes)
real(8):: muend_sim(nyrs,ntypes),muenddif(nyrs,ntypes),&
        & sigend_sim(nyrs,ntypes), sigenddif(nyrs,ntypes), &
        & muint_sim(nyrs,ntypes),muintdif(nyrs,ntypes), &
        & sigint_sim(nyrs,ntypes), sigintdif(nyrs,ntypes), &
        & muexit_sim(nyrs,ntypes),muentry_sim(nyrs,ntypes), muentrydif(nyrs,ntypes), muexitdif(nyrs,ntypes)

! internal variables
integer:: ss,yy,tt,bin,hist,idum,newbin,i,yr,tp
integer::   sal_sim(nyrs,ntypes,nbins), sal_data(nyrs,ntypes,nbins),&
          &               sal_sim_int(nyrs,ntypes,nbins),&
          &                 entry_sim(nyrs,ntypes,nbins),entry_data(nyrs,ntypes,nbins),&
          &                  exit_sim(nyrs,ntypes,nbins),exit_data(nyrs,ntypes,nbins)
real(8)::median_sal(nbins),seed,sal,raise,newsal,sig,mu,draw(1), temp(nbins)
real(8),allocatable:: raisedraw(:),saldraw(:),noincreasedraw(:)

   
    ! Compute number of simulated individuals per type for each year
    nend_sim(:,:)=sum(sal_sim,3)
    nenddif=nend_sim-nend
    
    nint_sim(:,:)=sum(sal_sim_int,3)
    nintdif=nint_sim-ninter
    
    nentry_sim(:,:)=sum(entry_sim,3)
    nentrydif=nentry-nentry_sim
    
    nexit_sim(:,:) =sum(exit_sim,3)
    nexitdif=nexit-nexit_sim

    median_sal=(/(i,i=1,nbins)/)*binsize-binsize/2

    ! Compute means of the distribution of salaries before and after entries and exits   
    do yy=1,nyrs
    do tt=1,ntypes
        muend_sim(yy,tt)=dot_product(sal_sim(yy,tt,:),median_sal(:))
        muend_sim(yy,tt)=muend_sim(yy,tt)/nend_sim(yy,tt)
        muint_sim(yy,tt)=dot_product(sal_sim_int(yy,tt,:),median_sal(:))
        muint_sim(yy,tt)=muint_sim(yy,tt)/nint_sim(yy,tt)
        muexit_sim(yy,tt)=dot_product(exit_sim(yy,tt,:),median_sal(:))
        muexit_sim(yy,tt)=muexit_sim(yy,tt)/nexit_sim(yy,tt)
        muentry_sim(yy,tt)=dot_product(entry_sim(yy,tt,:),median_sal(:))
        muentry_sim(yy,tt)=muentry_sim(yy,tt)/nentry_sim(yy,tt)

    enddo
    enddo
    muenddif=muend_sim-muend
    muintdif=muint_sim-muint
    muentrydif=muentry_sim-muentry
    muexitdif=muexit_sim-muexit
    
    ! Compute standard deviations of the distribution of salaries before and after entries and exits
    do yy=1,nyrs
    do tt=1,ntypes
        temp(:)=(median_sal(:)*1.0d-0-muend_sim(yy,tt))
        temp(:)=temp(:)**2
        temp(:)=sal_sim(yy,tt,:)*temp(:)
        sigend_sim(yy,tt)=sum(temp(:))/nend_sim(yy,tt)
        sigend_sim(yy,tt)=sqrt(sigend_sim(yy,tt))
        
        temp(:)=(median_sal(:)*1.0d-0-muint_sim(yy,tt))
        temp(:)=temp(:)**2
        temp(:)=sal_sim_int(yy,tt,:)*temp(:)
        sigint_sim(yy,tt)=sum(temp(:))/nint_sim(yy,tt)
        sigint_sim(yy,tt)=sqrt(sigint_sim(yy,tt))
    enddo
    enddo
    sigenddif=sigend_sim-sigend
    sigintdif=sigint_sim-sigint
    
    
    open(unit=1,file=trim(adjustl(pathoutput))//"simulationfit"//specification//".txt")
    write(1,'(A14,A14,19A18)')      "year","  type  ",&
                                  & "meanend_d","meanend_s","meanendfit",  "sdendfit", &
                                  & "nend_d","nend_sim",&
                                  & "muentry_d", "muentry_sim","muentrydif", &
                                  & "nentry_d","nentry_s",&
                                  & "muexit_d", "muexit_s","muexitdif",&
                                  & "nexit_d", "nexit_s",& 
                                  & "meanint_d","meanint_s  ", "meanintfit "
    do tt=1,ntypes
    do yy=1,nyrs
        write(1,fmt='(I4,I4,4f10.1,2I8,3f10.1,2I8,3f10.1,2I8,3f10.1)')1986*(1-tenure)+yy,tt,&
                    & muend(yy,tt),muend_sim(yy,tt),muenddif(yy,tt), sigenddif(yy,tt)    ,&  
                    & nend(yy,tt),nend_sim(yy,tt),&
                    & muentry(yy,tt),muentry_sim(yy,tt),muentrydif(yy,tt),&
                    & nentry(yy,tt),nentry_sim(yy,tt), &
                    & muexit(yy,tt),muexit_sim(yy,tt),muexitdif(yy,tt),&
                    & nexit(yy,tt),nexit_sim(yy,tt),&
                    & muint(yy,tt),muint_sim(yy,tt),muintdif(yy,tt)
                  
    enddo
    enddo
    close(1)
    
end subroutine
