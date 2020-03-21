			
			
			quietly su upi if $yr!=$ct1  & newspell!=1 // get the number of incumbents for that cell
			mat C=r(N)
			mat temp=(0.0)
			
			if $yr!=$ct1 & C[1,1]>1 {
				quietly tabstat noincrease if newspell!=1 , s(mean) c(s) save //fraction of incumbents whose salaries do not increase
				quietly tabstatmat tempp
				mat temp=tempp'		 
			}			
			if $j==1 {
				mat noincrease=temp
			}
			else{
				mat noincrease=noincrease\temp
			}
