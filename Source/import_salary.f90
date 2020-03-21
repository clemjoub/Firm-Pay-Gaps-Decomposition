subroutine Import_salary (nbins, year1, nyrs, nyrsdata, ntypes, &
                    &    ninter,muinter,siginter, &
                    &    nend, muend, sigend, &         !   end = end of period
                    &    nentry, muentry, sigentry, &
                    &    nexit, muexit, sigexit, &
                    &    alphagrowth, betagrowth, siggrowth, &
                    &    alphaexit, betaexit, beta2exit, &
                    &    noincrease,pathinput,entry_data,exit_data,sal_data)
implicit none


! arguments
integer:: nyrs,nyrsdata, ntypes, nstats,nbins,bin, year1
integer:: nend(nyrs,ntypes), nentry(nyrs,ntypes), nexit(nyrs,ntypes),ninter(nyrs,ntypes)
real(8):: muend(nyrs,ntypes), sigend(nyrs,ntypes),&
        & muentry(nyrs,ntypes), sigentry(nyrs,ntypes),&
        & muexit(nyrs,ntypes), sigexit(nyrs,ntypes), &
        & muinter(nyrs,ntypes), siginter(nyrs,ntypes), &
        & alphagrowth(nyrs,ntypes), betagrowth(nyrs,ntypes), siggrowth(nyrs,ntypes),sebetagrowth(nyrs,ntypes), &
        & noincrease(nyrs,ntypes), &
        & alphaexit(nyrs,ntypes), betaexit(nyrs,ntypes), beta2exit(nyrs,ntypes)

integer::   sal_data(nyrs,ntypes,nbins),&
          &               sal_data_int(nyrs,ntypes,nbins),&
          &                 entry_data(nyrs,ntypes,nbins),&
          &                  exit_data(nyrs,ntypes,nbins)
character*200:: pathinput,filename
character*2::ttchar


! internal
integer:: yy,tt,stat
character:: blah
    
    entry_data=0
    exit_data=0
    sal_data=0
    
    
    open(file=adjustl(trim(pathinput))//"simstats.txt",unit=1,position="rewind")
    read(1,*)
    do tt=1,ntypes
    !do tt=1,24
        do yy=1,nyrsdata
            write(ttchar,'(I2)') tt
            read(1,*) blah, nentry(yy,tt), muentry(yy,tt), sigentry(yy,tt), noincrease(yy,tt),  betagrowth(yy,tt),  sebetagrowth(yy,tt), alphagrowth(yy,tt), siggrowth(yy,tt),  betaexit(yy,tt), beta2exit(yy,tt), alphaexit(yy,tt)
        enddo
    enddo
       
    open(file=adjustl(trim(pathinput))//"salhist.txt",unit=2)
    bin=0
    stat=0      
    do while(bin.le.nbins.and.stat.ge.0)
    read(2,iostat=stat,fmt=*) bin, tt, yy, entry_data(yy-year1+1,tt,bin), exit_data(yy-year1+1,tt,bin), sal_data(yy-year1+1,tt,bin)
    !write(*,*) bin
    enddo
    close(2)
    
endsubroutine Import_salary
