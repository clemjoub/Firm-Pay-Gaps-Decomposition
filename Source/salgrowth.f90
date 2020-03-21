subroutine salgrowth(nyrs,ntypes,binsize,nbins,sal_sim,sal_sim_int,idum,noincrease,alphagrowth,betagrowth,siggrowth,tt,yy)
implicit none

integer:: nyrs 
integer:: ntypes 
integer:: nbins
real(8):: binsize

! Import_salary variables
integer:: nend(nyrs,ntypes), nentry(nyrs,ntypes), nexit(nyrs,ntypes),ninter(nyrs,ntypes)
real(8):: muend(nyrs,ntypes),     sigend(nyrs,ntypes),&
        & muentry(nyrs,ntypes), sigentry(nyrs,ntypes),&
        & muexit(nyrs,ntypes),   sigexit(nyrs,ntypes), &
        & muint(nyrs,ntypes), sigint(nyrs,ntypes)
real(8):: alphagrowth(nyrs,ntypes), &
         & betagrowth(nyrs,ntypes), &
         &  siggrowth(nyrs,ntypes), &
         & noincrease(nyrs,ntypes)

! Output statistics
integer:: nend_sim(nyrs,ntypes),nenddif(nyrs,ntypes),nint_sim(nyrs,ntypes),nintdif(nyrs,ntypes)
real(8):: muend_sim(nyrs,ntypes),muenddif(nyrs,ntypes),&
        & sigend_sim(nyrs,ntypes), sigenddif(nyrs,ntypes), &
        & muint_sim(nyrs,ntypes),muintdif(nyrs,ntypes),&
        & sigint_sim(nyrs,ntypes), sigintdif(nyrs,ntypes)

! internal variables
integer:: ss,yy,tt,bin,hist,idum,newbin,i,yr,tp,binfreq
integer::   sal_sim(nyrs,ntypes,nbins),&
          &               sal_sim_int(nyrs,ntypes,nbins),&
          &                 entry_sim(nyrs,ntypes,nbins),&
          &                  exit_sim(nyrs,ntypes,nbins)
real(8)::median_sal(nbins),sal,raise,newsal,sig,mu,draw(1), temp(nbins)
real(8),allocatable:: raisedraw(:),saldraw(:),noincreasedraw(:)

 
        GROWTHLOOP: do bin=1,nbins
	    binfreq=sal_sim(yy-1,tt,bin)
            if (binfreq.gt.0) then
            
                allocate(raisedraw(binfreq))
                allocate(saldraw(binfreq))
                allocate(noincreasedraw(binfreq))
            
                call gasdev(idum,raisedraw,binfreq)
                call ran1(idum,saldraw,binfreq)
                call ran1(idum,noincreasedraw,binfreq)
            
                ! Consider each observation in each bin of last year's distribution
                BINSHIFTS: do ss=1,binfreq
                    ! Draw a salary from the bin
                    sal=(bin-1+saldraw(ss))*binsize
                    ! Draw salary growth and compute new salary
                    mu=alphagrowth(yy,tt)+sal*betagrowth(yy,tt)
                    sig=siggrowth(yy,tt)
                    raise=1+exp(raisedraw(ss)*sig+mu)
                    if (noincreasedraw(ss).lt.noincrease(yy,tt)) raise= 1       
                    newsal=sal*raise
                    ! Identify new salary bin for that simulated employee
                    newbin=min(nbins,1+int(newsal/binsize))
                    sal_sim_int(yy,tt,newbin)=sal_sim_int(yy,tt,newbin)+1
                    
                enddo BINSHIFTS
                
                deallocate(raisedraw)
                deallocate(saldraw)  
                deallocate(noincreasedraw) 
                 
            endif
           
         enddo GROWTHLOOP
         

endsubroutine salgrowth
