
#
# Script to generate automatically the fillQ routines
#

GenerateCodonsList <- function() {

	b <- c("T", "C", "A", "G")
	
	codon <- character(3)
	codons <- list()
	n <- 0
	for(i in b) {
	
		for(j in b) {

			for(k in b) {
			
				if(i == "T" & j == "A" & k == "A") next
				if(i == "T" & j == "A" & k == "G") next
				if(i == "T" & j == "G" & k == "A") next
				
				codon[1] <- i
				codon[2] <- j
				codon[3] <- k
				
				n <- n+1
				codons[[n]] <- codon
			}
		}
	}
	
	return(codons)
}

CountChanges <- function(c1, c2) {

	n <- 0
	if(c1[1] != c2[1]) n <- n+1
	if(c1[2] != c2[2]) n <- n+1
	if(c1[3] != c2[3]) n <- n+1
	
	return(n)
}

IsTransition <- function(c1, c2) {

	if(c1[1] != c2[1]) n <- 1
	if(c1[2] != c2[2]) n <- 2
	if(c1[3] != c2[3]) n <- 3
	
	if(c1[n] == "A" & c2[n] == "G") return(TRUE)
	if(c1[n] == "G" & c2[n] == "A") return(TRUE)
	if(c1[n] == "C" & c2[n] == "T") return(TRUE)
	if(c1[n] == "T" & c2[n] == "C") return(TRUE)
	
	return(FALSE)
}

IsSynonimous <- function(i1, i2) {

#    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
#  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG

	AAs <- "FFLLSSSSYYCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

	a1 <- substr(AAs, i1, i1)
	a2 <- substr(AAs, i2, i2)
	return(a1 == a2)
}


GenerateQMatrix <- function(omega.eq.1=FALSE, k="k", omega="omega", p="codonFreq", idx.base=0) {

	codons <- GenerateCodonsList()
	
	Q <- matrix("0", 61, 61)
 
	for(from in 1:60) {
	
		for(to in (from+1):61) {
		
			#cf <- paste(codons[[from]], collapse="")
			#ct <- paste(codons[[to]], collapse="")

			if(CountChanges(codons[[from]], codons[[to]]) > 1) next
			
			ist <- IsTransition(codons[[from]], codons[[to]])
			
			if(IsSynonimous(from, to)) {
			
				if(ist) {
					#cat(sprintf("%2d %2d %s %s k\n", from, to, cf, ct))
					Q[from, to] <- k
				} else {
					#cat(sprintf("%2d %2d %s %s 1\n", from, to, cf, ct))
					Q[from, to] <- "1"
				}
			
			} else {
			
				if(ist) {
					#cat(sprintf("%2d %2d %s %s omega*k\n", from, to, cf, ct))
					Q[from, to] <- ifelse(omega.eq.1, k, paste(k,"*",omega, sep=""))
				} else {
					#cat(sprintf("%2d %2d %s %s omega\n", from, to, cf, ct))
					Q[from, to] <- ifelse(omega.eq.1, "1", omega)
				}
			}
		}
	}
	
   # Complete the matrix
    for(from in 1:60) {

        for(to in (from+1):61) {

            Q[to, from] <- Q[from, to]

            if(Q[from, to] == "0") next
            if(Q[from, to] == "1") {

                Q[from, to] <- paste(p, "[",   to-1+idx.base, "]", sep="")
                Q[to, from] <- paste(p, "[", from-1+idx.base, "]", sep="")

            } else {

                Q[from, to] <- paste(Q[from, to], "*", p, "[",   to-1+idx.base, "]", sep="")
                Q[to, from] <- paste(Q[to, from], "*", p, "[", from-1+idx.base, "]", sep="")

            }
        }
    }

	for(i in 1:61) {
	
		Q[i, i] <- gsub("-0", "", paste(c("", Q[i, -i]), collapse="-"))
	}
	
	cn <- character(61)
	for(i in 1:61) cn[i] <- paste(codons[[i]], collapse="")
	colnames(Q) <- cn
	rownames(Q) <- cn
	
	return(Q)
}


GenerateSMatrix <- function(omega.eq.1=FALSE, k="k", omega="omega", p="codonFreq") {

	codons <- GenerateCodonsList()
	
	S <- matrix("0", 61, 61)
 
	for(from in 1:60) {
	
		for(to in (from+1):61) {

			if(CountChanges(codons[[from]], codons[[to]]) > 1) next
			
			ist <- IsTransition(codons[[from]], codons[[to]])
			
			if(IsSynonimous(from, to)) {
			
				if(ist) {
					S[from, to] <- k
				} else {
					S[from, to] <- "1"
				}
			
			} else {
			
				if(ist) {
					S[from, to] <- ifelse(omega.eq.1, k, paste(k,"*",omega, sep=""))
				} else {
					S[from, to] <- ifelse(omega.eq.1, "1", omega)
				}
			}

            S[to, from] <- S[from, to]
		}
	}

	# The S diagonal elements are premultiplied by the corresponding codon frequency value
	for(i in 1:61) {
	
		S[i, i] <- ""
		for(j in 1:61) {
		
			if(i == j) next
			if(S[i, j] == "0") next
			if(S[i, j] == "1") {
				S[i, i] <- paste(S[i, i], "-", p, "[", j-1, "]", sep="")
			} else {
				S[i, i] <- paste(S[i, i], "-", S[i, j], "*", p, "[", j-1, "]", sep="")
			}
		}
	}

	cn <- character(61)
	for(i in 1:61) cn[i] <- paste(codons[[i]], collapse="")
	colnames(S) <- cn
	rownames(S) <- cn
	
	return(S)
}

WriteFillQ <- function(c.source.file, mode="w", omega.eq.1=FALSE) {

	fd <- file(c.source.file, mode)

	if(omega.eq.1) {
	    cat("double TransitionMatrix::fillMatrix(double aK)\n{\n", sep="", file=fd)
	} else {
	    cat("double TransitionMatrix::fillMatrix(double aOmega, double aK)\n{\n", sep="", file=fd)
	}

	Q <- GenerateQMatrix(omega.eq.1, k="aK", omega="aOmega", p="mCodonFreq", idx.base=0)

	for(i in 1:61) {

		for(j in 1:61) {
		
			if(Q[i, j] == "0") next
			
			cat(sprintf("\tmQ[%4d] = %s;\n", j+i*61-62, Q[i, j]), sep="", file=fd)
		}
	}

	S <- GenerateSMatrix(omega.eq.1, k="aK", omega="aOmega")

	cat("\n\t// Compute the scale factor\n", sep="", file=fd)
	cat("\tdouble scale_q =     ", sep="", file=fd)

	is.first <- TRUE
	for(i in 1:60) {
		for(j in (i+1):61) {

			if(S[i, j] == "0") next
			if(!is.first) cat("\n\t             + ", sep="", file=fd)
			is.first <- FALSE
			if(S[i, j] == "1") {
				cat("  mCodonFreq[",i-1,"]*mCodonFreq[", j-1,"]", sep="", file=fd)
			} else {
				cat(S[i, j],"*mCodonFreq[",i-1,"]*mCodonFreq[", j-1,"]", sep="", file=fd)
			}
		}
	}

	cat(";\n\n\treturn scale_q*2.;\n", sep="", file=fd)

    cat("}\n\n", sep="", file=fd)

	close(fd)
}


WriteFillS <- function(c.source.file, mode="w", omega.eq.1=FALSE) {

	fd <- file(c.source.file, mode)

	if(omega.eq.1) {
	    cat("double TransitionMatrix::fillMatrix(double aK)\n{\n", sep="", file=fd)
	} else {
	    cat("double TransitionMatrix::fillMatrix(double aOmega, double aK)\n{\n", sep="", file=fd)
	}

	S <- GenerateSMatrix(omega.eq.1, k="aK", omega="aOmega", p="mCodonFreq")

	for(i in 1:61) {

		#for(j in 1:61) {
		for(j in 1:i) {
		
			if(S[i, j] == "0") next
			
			cat(sprintf("\tmS[%4d] = %s;\n", j+i*61-62, S[i, j]), sep="", file=fd)
		}
		cat("\n", sep="", file=fd)
	}

	cat("\n\t// Compute the scale factor\n", sep="", file=fd)
	cat("\tdouble scale_q =     ", sep="", file=fd)

	is.first <- TRUE
	for(i in 1:60) {
		for(j in (i+1):61) {

			if(S[i, j] == "0") next
			if(!is.first) cat("\n\t             + ", sep="", file=fd)
			is.first <- FALSE
			if(S[i, j] == "1") {
				cat("  mCodonFreq[",i-1,"]*mCodonFreq[", j-1,"]", sep="", file=fd)
			} else {
				cat(S[i, j],"*mCodonFreq[",i-1,"]*mCodonFreq[", j-1,"]", sep="", file=fd)
			}
		}
	}

	cat(";\n\n\treturn scale_q*2.;\n", sep="", file=fd)

    cat("}\n\n", sep="", file=fd)

	close(fd)
}

WriteFillQR <- function(R.source.file, mode="w", omega.eq.1=FALSE) {

	fd <- file(R.source.file, mode)

	if(omega.eq.1) {
	    cat("fillQ <- function(k) {\n", sep="", file=fd)
	} else {
	    cat("fillQ <- function(omega, k) {\n", sep="", file=fd)
	}

	cat("\tp <- rep(1/61, 61)\n\n", sep="", file=fd)

	Q <- GenerateQMatrix(omega.eq.1, k="k", omega="omega", p="p", idx.base=1)

	cat("\tQ <- matrix(0, 61, 61)\n\n", sep="", file=fd)
	for(i in 1:61) {

		cat("\tQ[", i, ", ", i, "] <- ", Q[i,i], "\n", sep="", file=fd)

		if(i == 61) break
		for(j in (i+1):61) {
		
			if(Q[i, j] == "0") next
			
			cat("\tQ[", i, ", ", j, "] <- Q[", j, ", ", i, "] <- ", Q[i, j], "\n", sep="", file=fd)
		}
	}

	S <- GenerateSMatrix(omega.eq.1, k="k", omega="omega")

	cat("\tscale_q <-     ", sep="", file=fd)

	is.first <- TRUE
	for(i in 1:61) {
		for(j in 1:61) {

			if(i == j) next
			if(S[i, j] == "0") next
			if(!is.first) cat("\n\t             + ", sep="", file=fd)
			is.first <- FALSE

			if(S[i, j] == "1") {
				cat("p[",i,"]*p[", j,"]", sep="", file=fd)
			} else {
				cat(S[i, j],"*p[",i,"]*p[", j,"]", sep="", file=fd)
			}
		}
	}			
	
	
	cat("\n\n\treturn(list(scale=scale_q, Q=Q))\n", sep="", file=fd)

    cat("}\n\n", sep="", file=fd)

	close(fd)
}

VisualizeStructure <- function(Q) {

	V <- ifelse(Q != "0", 1, 0)
	W <- V[1:61,61:1]
	image(1:61, 1:61, W, col=c("white", "blue"), xaxt='n', yaxt='n', xlab="", ylab="")
	axis(2, at=c(0,62), label=rep("",2))
	axis(3, at=c(0,62), label=rep("",2))
}

#
# Transform the Q matrix into a numeric one substituting the values given
#
GenerateNumericMatrix <- function(Q, k, omega, codonFreq) {

	res <- matrix(0, 61, 61)
	for(i in 1:61) {
		for(j in 1:61) {
			if(Q[i, j] == "0") next
			cat(i, j, Q[i, j], "\n")
			e <- parse(text=Q[i, j])
			res[i, j] <- eval(e)
		}
	}
	return(res)
}

#
# Function to be called to generate fillQ.cpp
#
generateCpp.old <- function() {
	if(.Platform$OS.type == "windows") {
		file.src <- "C:/mv/Projects/Selectome/codeml/slim/trunk/Codeml_Variants/Fastcodeml/FillMatrix.cpp"
	} else {
		file.src <- "/users/mvalle/codeml/slim/trunk/Codeml_Variants/Fastcodeml/FillMatrix.cpp"
	}
	cat("// Automatically generated file using genQ.r script, don't modify!\n\n", sep="", file=file.src)
	cat("#include \"TransitionMatrix.h\"\n\n", sep="", file=file.src, append=TRUE)
	WriteFillQ(file.src, mode="a", omega.eq.1=FALSE)
	WriteFillQ(file.src, mode="a", omega.eq.1=TRUE)

	# WriteFillQR("C:/temp/ex.r")
	# source("C:/temp/ex.r")
}

#
# Function to be called to generate fillQ.cpp (new experimental version)
#
generateCpp <- function(file.src="./FillMatrix.cpp") {

	cat("// File automatically generated by the genQ.r script, don't modify!\n\n", sep="", file=file.src)
	cat("#include \"TransitionMatrix.h\"\n\n", sep="", file=file.src, append=TRUE)
	cat("#ifndef USE_LAPACK\n\n", sep="", file=file.src, append=TRUE)
	WriteFillQ(file.src, mode="a", omega.eq.1=FALSE)
	WriteFillQ(file.src, mode="a", omega.eq.1=TRUE)
	cat("#else\n\n", sep="", file=file.src, append=TRUE)
	WriteFillS(file.src, mode="a", omega.eq.1=FALSE)
	WriteFillS(file.src, mode="a", omega.eq.1=TRUE)
	cat("#endif\n", sep="", file=file.src, append=TRUE)
}

generateCpp()
