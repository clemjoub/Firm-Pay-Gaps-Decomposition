			mat temp=(0.0,0.0,0.0)
			
			global Y = "salarycurr000"
			
			quietly su $Y if entries ==1 // average salary among entrants
			if r(N)==1 {
				mat temp=(r(N),r(mean),0.0)
			}
			else if r(N)>1 {
				mat temp=(r(N),r(mean),r(sd))
			}
			
			if $j==1 { //append parameters for this cell to the others
				mat entry=temp
			}
			else{
				mat entry=entry\temp
			}
