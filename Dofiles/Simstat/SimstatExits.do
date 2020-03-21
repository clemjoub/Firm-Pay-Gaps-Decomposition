			
			quietly su upi					//number of employees in the cell
			mat A=r(N)			
			quietly su upi if endspell==1 	//number of exiters
			mat B=r(N) 
			
			mat M=0			
			if A[1,1]>0 { 					// need non-empty cell
				quietly su salarycurr000 
				mat M=r(mean)				//store mean salary in the cell
			}
			mat temp=(0,0,M)
			
			
			global Y = "endspell"			// dummy for exiting
			global X = "salarycurr000centered"	// centered salary variable
			
			if A[1,1]>20 & B[1,1]>5 { 	//require cellsize of at least 20 and at least 5 exits to run regression
				quietly reg $Y $X 
				local rc=_rc			
				if `rc'==0 {			//if no problem with regression
					mat temp=(e(b),M)
				}		
			}
			
			else if A[1,1]>0{ 			//if cellsize too small, just use fraction of exiters
				quietly tabstat $Y, s(mean) c(s) save
				quietly tabstatmat temp	
				mat temp=(0,temp',M)	
			}
			
			if $j==1 {					//append parameters for this cell to the rest
				mat exit=temp
			}
				else{
				mat exit=exit\temp
			}

			
