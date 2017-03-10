#' Class "optsol_dGSMEF"  
#'
#' Structure of class \code{optsol_dGSMEF}. Objects of that class are returned by the function \code{dGSMEF}. 
#'
#' @param concmat Object of class "matrix" containing concentrations of extracellular metabolite and recombinant protein.
#' @param exRxn Object of class "matrix" containing names of exchange reactions for the EC metabolites.
#' @param tmVec Object of class "numeric" Vector of time points
#' @param bmVec Object of class "numeric" Vector of biomass concentrations
#' @param flux Object of class "fluxDistribution" containing the solutons flux distributions.
#' @param redCost Object of class "matrix" containing the solutons flux distributions.
#' @param Rxn Object of class "character" Vector of names of exchange reactions for the EC metabolites.
#'
#' @name optsol_dGSMEF-class
#' @rdname optsol_dGSMEF-class
#' @exportClass optsol_dGSMEF
#'
#' @seealso \code{dGSMEF}
#'
#' @examples
#' showClass("optsol_dGSMEF")

setClass("optsol_dGSMEF",
      slots = list(
		concentrationMatrix="numeric", # Matrix of extracellular metabolite concentrations
		excRxnNames="character",       # Names of exchange reactions for the EC metabolites
		timeVec="numeric",             # Vector of time points
		biomassVec="numeric",           # Vector of biomass values
		fluxMatrix = "numeric", 		# the flux distribution
		redCostMatrix = "numeric", 		# the reduced cost distribution
		RxnNames="character"			# Names of all reactions
      )
)

#showClass("optsol_dGSMEF")

#' Constructor method for optsol_dGSMEF class
#'
#' @name optsol_dGSMEF
#' @rdname optsol_dGSMEF-class
optsol_dGSMEF <- function(concmat,exRxn,tmVec,bmVec,flux,redCost,Rxn) {
    if (missing(concmat) || 
        missing(exRxn) ||
        missing(tmVec)  ||
        missing(bmVec)  ||
        missing(flux)  ||
        missing(redCost) ||
        missing(Rxn)
       ) {
        stop("Not enough arguments for creating an object of class optsol_dGSMEF!")
    }

    new("optsol_dGSMEF", 
		fluxMatrix           = flux,
		redCostMatrix        = redCost,
		concentrationMatrix  = concmat,
		excRxnNames          = exRxn,
		timeVec              = tmVec,  
		biomassVec           = bmVec,
		RxnNames             = Rxn
		)
}

#' @name optsol_dGSMEF-class
#' @rdname optsol_dGSMEF-class
#'
#' @param x An object of class \code{optsol_dGSMEF}.  
#' @param y not used but kept for compitability with generic plot.
#' @param plotRxns List of reaction id's to be ploted (metabolite concentrations).
#' @param plotFlux List of reaction id's to be ploted (reaction fluxes).
#'
#' @exportMethod plot

setMethod("plot", signature("optsol_dGSMEF","missing"),
          function(x,y,plotRxns=NULL,plotFlux=NULL) {
                if(missing(plotRxns) && missing(plotFlux)){
                   plot(spline(x@timeVec,x@biomassVec, n = 201, method = "natural"), col = 1
			                    ,main='Cell density',xlab='Time(hrs)',ylab="X(g/L)",type="l",lwd=2);
					points(x@timeVec,x@biomassVec, col = "red",lwd=2);
                   }
                else if (!missing(plotRxns)) {
					def.par <- par(no.readonly = TRUE);
					layout(matrix(c(1,2,1,2), 2, 2, byrow = TRUE))
					#layout.show(2);
					# first plot biomass
					plot(spline(x@timeVec,x@biomassVec, n = 201, method = "natural"), col = 1
			                    ,main='Cell density',xlab='Time(hrs)',ylab="X(g/L)",type="l",lwd=2);
					points(x@timeVec,x@biomassVec, col = "red",lwd=2);
					for ( i in 1:length(plotRxns) ){
						plotInd=(x@excRxnNames %in% plotRxns[i]);
						if(plotRxns[i] == 'V'){
								par(new=T);
								ymin <- min(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns[i]], function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
								ymax <- max(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns[i]], function(x) max(x, na.rm = TRUE)), na.rm = TRUE) 
								plot(spline(x@timeVec, x@concentrationMatrix[plotInd], n = 201, method = "natural"),
                                     type="l",lwd=2,col =i,axes=F, xlab = "", ylab = "",ylim=c(ymin-0.3,ymax));
								axis(4);
								mtext("Volume(L)",side=4,cex=0.8);
								#points(x@timeVec,x@concentrationMatrix[plotInd], col = i,pch=2,lwd=2);
								legend(1.8,ymax, c("Biomass","Volume"), col=c(1,i), lty=1, pch=c(NA,NA),text.font = 1.2, cex = 0.6);
								legend(1.8,ymax, c("",""), col=c("red",i), lty=0,pch=c(1,NA),text.font = 1.2, cex = 0.6,bty='n');
								plotRxns <- plotRxns[-i];
						}
					}

					# plot concentrations
					##plot(x@timeVec,2*x@biomassVec,main='Biomass',xlab='Time',ylab=ylab);
					## define min/max ()plot(x@timeVec,
					ymin <- min(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns], function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
					ymax <- max(sapply(x@concentrationMatrix[x@excRxnNames %in% plotRxns], function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
					for ( i in 1:length(plotRxns) ){
						plotInd=(x@excRxnNames %in% plotRxns[i]);
						#print( x@concentrationMatrix[plotInd]);
						if(i==1){
								   plot(spline(x@timeVec, x@concentrationMatrix[plotInd], n = 201, method = "natural"),
                                   type="l", col =i,main="Concentrations",ylab = "mmol/L",xlab='Time(hrs)'
                                   ,ylim=c(ymin,ymax));
                        }
						else{
							if(plotRxns[i] == 'pH'){
								par(new=T);
								ymin <- 6;
								ymax <- 8;
								plot(spline(x@timeVec, x@concentrationMatrix[plotInd], n = 201, method = "natural"),
                                   type="l", col =i,axes=F, xlab = "", ylab = "",ylim=c(ymin,ymax));
								axis(4);
								mtext("pH",side=4,cex=0.8);
							}
							else{ lines(spline(x@timeVec, x@concentrationMatrix[plotInd], n = 201, method = "natural"), col =i);}
						}
					}
					legend("topright", plotRxns, col=1:length(plotRxns), lty=1, text.font = 1.2, cex = 0.6);
                }
				else{
					if(length(x@fluxMatrix) <= length(x@RxnNames)){ stop("No flux distribution argument for creating an object of class optsol_dGSMEF!"); }
					# plot flux distributions
					##plot(x@timeVec,2*x@biomassVec,main='Biomass',xlab='Time',ylab=ylab);
					## define min/max
					ymin <- min(sapply(x@fluxMatrix[x@RxnNames %in% plotFlux], function(x) min(x, na.rm = TRUE)), na.rm = TRUE)
					ymax <- max(sapply(x@fluxMatrix[x@RxnNames %in% plotFlux], function(x) max(x, na.rm = TRUE)), na.rm = TRUE)  
					for ( i in 1:length(plotFlux) ){
						plotInd=(x@RxnNames %in% plotFlux[i]);
						#print( x@concentrationMatrix[plotInd]);
						if(i==1){
								   plot(spline(x@timeVec, x@fluxMatrix[plotInd], n = 201, method = "natural"),
                                   type="l", col =i,main="Flux distributions",ylab = "mmol/gDCW*hrs",xlab='Time(hrs)'
                                   ,ylim=c(ymin,ymax));
                        }
						else{
							lines(spline(x@timeVec, x@fluxMatrix[plotInd], n = 201, method = "natural"), col =i);
						}
					}
					legend("topright", plotFlux, col=1:length(plotFlux), lty=1, text.font = 1.2, cex = 0.6);
					
				}
                  
		#if (!missing(plotRxns)){		   }
          }
)

