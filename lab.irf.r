lab.irf <- function(ip, x = NULL, D = NULL)
 {	# this is a modified irf function in irtoys including D  
	if (is.null(x))
        x = seq(-4, 4, length = 101)
	if (is.null(dim(ip))) 
        dim(ip) = c(1, 3)
	if (is.null(D))	
		D = 1.0
	ni <- dim(ip)[1]
	a <- ip[,1]; 
	b <- ip[,2]; 
	g <- ip[,3]
 	f = sapply(1:ni, function(i) g[i] + (1-g[i])* 1/(1 + exp(-1*D*a[i]*(x-b[i]))))
	r = list(x = x, f = f, d=D)
	class(r)="irf"
	return(r)
}	