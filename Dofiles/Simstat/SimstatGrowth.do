		
			
			
			
			
			mat temp=(0.0,0.0,0.0,0.0)
			quietly su upi if logincr_tr!=. & lagsal!=. & newspell!=1 //number incumbents with non-zero salary increases
			mat F=r(N)
			
			global Y= "logincr_tr"
			global X= "lagsal"
			
			if F[1,1]>1 {	/*need at least 1 observation*/
				
				quietly reg $Y $X if newspell!=1
				mat D=e(b)
				mat E=e(V)
				
				if  sqrt(E[1,1])>0 { // non-degenerate case
					mat temp=(D[1,1],sqrt(E[1,1]),D[1,2],e(rmse))			
				}
				else{               // degenerate case
					quietly reg $Y if newspell!=1
					mat D=e(b)
					mat E=e(V)
					mat temp=(0.0,0.0,D[1,1],e(rmse)) 
								
				}
			}
			
	
			if $j==1{ // append estimates for this group to "growth" matrix
					mat growth=temp
				}
				else{
					mat growth=growth\temp
			}
