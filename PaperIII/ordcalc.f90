! This program calculates the experimental order of convergence given two grid resolutions and the respective error
program ordcalc
        real :: h1, h2
        real :: err1, err2
        real :: ord

        write (*,*) "First h?:"
        read (*,*) h1
        write (*,*) "First error?"
        read (*,*) err1
        write (*,*) "Second h?:"
        read(*,*) h2
        write (*,*) "Second error:"
        read(*,*) err2


        ord = log(err1 / err2) / log(h1 / h2)
        write (*,*) "Experimental order of accuracy is", ord


end program ordcalc
