#' A dynamic genome-scale modelling of E. Coli fermentation (dGSMEF)
#'
#' \code{dGSMEF} performs a dynamic flux balance analysis under batch, fed-batch and induction cultures.
#' 
#' The function is implemented with static optimization approach (SOA)
#' First, the simulation time is divided into small periods which are 
#' assumed to be in quasi-steady state; Second, for every time step, 
#' an FBA problem is solved and the fluxes are integrated over the 
#' time period and extracellular concentrations are calculated, accordingly.
#' Furthermore, the flux optimizaiton is achieved by parsimonius FBA (pFBA),
#' which minimize all fluxes within the solution space of the growth optimum,
#' to minimize the requirement for enzyme expression.
#'
#' In addition, for the bi-objective optimization of both maximum growth and
#' recombinant protein synthesis rates, a novel proteome resource allocation
#' theory is applied. According to this theory, the growth reduction induced
#' by heterologous protein expression is a simple consequence of proteome
#' allocations. And more heterologous proteins, less growth rate.
#'
#' @param model Sybil model structure (class \code{\link[sybil]{modelorg}}).
#' @param substrateRxns List of exchange reaction names for substrates
#'                       initially in the media that may change (e.g. not
#'                       h2o or co2).
#' @param initConcentrations Initial concentrations of substrates (in the same
#'                       structure as substrateRxns).
#' @param initBiomass Initial biomass (must be non zero).
#' @param timeStep Time step size.
#' @param nSteps Maximum number of time steps.
#' @param u_fix	Fixed specific growth rate during fed-batch phase before induction.
#' @param x_ind	Biomass centration at induction.
#' @param feedSubstrateRxns	List of exchange reaction names for substrates in the feeding solution. Default: "EX_glc__D(e)".
#' @param feedConcentrations Substrate concentrations in the feeding solution (same structure as feedSubstrateRxns). Default: 2775.3 mmol/L, 500 g/L.
#' @param yield_rate Biomass yield rate on glucose, unit: g DW/g glucose.
#' @param ngam	Non-growth associated maintenance, to characterize the metabolic burden for the recombinant cell.
#' @param initRatio	Intial ratio of resources allocated to the recombinant protein synthesis in the beginning time of induction.
#' @param feedRate_ind	Feeding rate during induction, unit: ml/L.
#' @param initVolume Initial total volume of bioreactor, default: 1.
#' @param initpH Initial pH in the bioreactor, default: 7.3.
#' @param exclUptakeRxns List of uptake reactions whose substrate
#'                       concentrations do not change (Default =
#'                       {'EX_co2(e)','EX_h2o(e)','EX_h(e)'}).
#' @param rcd Boolean. Save the reduced costs. Default: FALSE.
#' @param fld Boolean. Save the resulting flux distribution. Default: FALSE.
#' @param verboseMode An integer value indicating the amount of output to stdout: 
#'						0: nothing, 
#'						1: status messages, 
#'						2: like 1 plus a progress indicator, 
#'						3: a table containing the reaction id's and the corresponding min max values.
#'						Default: 2.
#'
#' @return  returns class \code{optsol_dGSMEF}  
#' \itemize{
#' \item concentrationMatrix:    Matrix of extracellular metabolite concentrations.
#' \item excRxnNames:            Names of exchange reactions for the EC metabolites.
#' \item timeVec:                Vector of time points.
#' \item biomassVec:             Vector of biomass values.
#' \item fluxMatrix:             Matrix of flux distributions.
#' \item redCostMatrix:          Matrix of reduced costs for each reaction.
#' \item RxnNames:               Names of all reactions for the E.Coli metabolites.
#' }
#'
#' @seealso \code{\link[sybil]{modelorg}}, \code{optsol_dGSMEF}.
#'
#' @examples
#' library(sybil)
#' library(sybilSBML)
#' library(glpkAPI)
#' library(deSolve)
#' library(dGSMEF)
#' model = readSBMLmod("../data/HNC47_without_fbc.xml");
#' lowbnd(model)[react_id(model) == 'EX_glc__D(e)'] = -9.58; ## for HNC47 with Fc
#' lowbnd(model)[react_id(model) == 'EX_o2(e)']=-11.5; # for HNC47 with Fc
#' 
#' ## Add Fc synthesis reaction into the model
#' model <- addReact(model,id="RecombinantProtein",
#'				met=c('ala__L[c]','cys__L[c]','asp__L[c]','glu__L[c]','phe__L[c]',
#'                    'gly[c]','his__L[c]','ile__L[c]','lys__L[c]','leu__L[c]',
#'                    'met__L[c]','asn__L[c]','pro__L[c]','gln__L[c]','arg__L[c]',
#'                    'ser__L[c]','thr__L[c]','val__L[c]','trp__L[c]','tyr__L[c]',
#'                     'atp[c]','adp[c]','pi[c]'),
#'				Scoef=c(-c(6,6,9,17,9,11,5,4,17,17,3,11,19,11,6,24,14,23,4,8),
#'                      -4.306*224,964.544,965.544),
#'				lb=0,ub=1000,obj=0);
#'
#' ## Define FDM recipe
#' substrateRxns = c('EX_glc__D(e)','EX_nh4(e)','EX_pi(e)','EX_so4(e)','EX_mg2(e)',
#'                   'EX_k(e)','EX_ca2(e)','EX_cl(e)','EX_fe2(e)','EX_mn2(e)',
#'                   'EX_cu2(e)','EX_zn2(e)','EX_mobd(e)');
#' initConcentrations = c(55.56,2000,67.6,63.08,2,67.6,1000,1000,2000,0.09171,1000,0.05912,0.01078);
#' 
#' ## EX_ni2(e) and EX_cobalt2(e) both are required to support the growth of cell
#' lowbnd(model)[react_id(model) == 'EX_cobalt2(e)'] = -1000;
#' lowbnd(model)[react_id(model) == 'EX_fe3(e)'] = 0;
#' lowbnd(model)[react_id(model) == 'EX_tungs(e)'] = 0;
#' lowbnd(model)[react_id(model) == 'EX_ni2(e)'] = -1000;
#' lowbnd(model)[react_id(model) == 'EX_sel(e)'] = 0;
#' lowbnd(model)[react_id(model) == 'EX_slnt(e)'] = 0;
#'
#' ## 10% YE + 25% glucose for feeding during induction
#' AA <- c('EX_ala__L(e)','EX_arg__L(e)','EX_asn__L(e)','EX_asp__L(e)','EX_cys__L(e)','EX_gln__L(e)',
#'         'EX_glu__L(e)','EX_gly(e)','EX_his__L(e)','EX_ile__L(e)','EX_leu__L(e)','EX_lys__L(e)',
#'         'EX_met__L(e)','EX_phe__L(e)','EX_pro__L(e)','EX_ser__L(e)','EX_thr__L(e)',
#'          'EX_trp__L(e)','EX_tyr__L(e)','EX_val__L(e)')
#'
#' conc_aa_feeding <- 2.5*c(15.018,8.795,0.000,6.655,0.000,0.000,16.459,7.423,0.000,6.955,
#'                          11.957,5.917,2.207,5.588,3.082,4.451,6.036,0.000,1.175,9.289)	
#' feedSubstrateRxns = c('EX_glc__D(e)','EX_nh4(e)','EX_so4(e)','EX_mg2(e)','EX_fe2(e)',
#'                       'EX_mn2(e)','EX_cu2(e)','EX_zn2(e)','EX_mobd(e)',AA);
#' feedConcentrations = c(2775.3,0.00924,40.23174,40,0.0719,0.09171,0.00881,0.05912,0.01078);
#' feedConcentrations = feedConcentrations/2;
#' feedConcentrations = c(feedConcentrations,conc_aa_feeding);
#'
#' Ec_df2 <- dGSMEF(model,substrateRxns = substrateRxns,initConcentrations = initConcentrations,
#'                  initBiomass=0.0142,u_fix=0.15,x_ind= 31.63,
#'					feedSubstrateRxns = feedSubstrateRxns,
#'                  feedConcentrations = feedConcentrations,
#'					yield_rate = 0.6586, # 25% glucose + 10% YE
#'					ngam = 6.34,
#'					initRatio = 0.70,
#'					feedRate_ind = 2.8,
#'					timeStep=0.25,nSteps=180,fld=TRUE);
#'
#' ## compare predictions with experimental results during induction phase
#' locs <- Ec_df2@biomassVec >= 31.63
#' ind_time = Ec_df2@timeVec[locs][1]
#' plot(spline(Ec_df2@timeVec[locs]-ind_time,Ec_df2@biomassVec[locs], n = 201, method = "natural"), 
#'      col = 2,main='Concentrations',xlab='Time(hrs)',ylab="Biomass(g/L)",ylim = c(31,42),
#'      type="l",lwd=2);
#'
#' ind_time = 0
#' points(c(ind_time,ind_time+4,ind_time+8,ind_time+12,ind_time+16,ind_time+20),
#'        c(31.63,37.84164381,40.34110603,40.00462097,41.14031209,41.68691719),pch=2,col=2)
#'
#' par(new=T);	
#'
#' ind_time = Ec_df2@timeVec[locs][1]	
#' plot(spline(Ec_df2@timeVec[locs]-ind_time, 
#'             Ec_df2@concentrationMatrix[Ec_df2@excRxnNames %in% "Protein"][locs]*25163.41/1000, 
#'             n = 201, method = "natural"),axes=F,col = 4,xlab="",ylab="",type="l",lwd=2);
#'
#' ind_time = 0
#' points(c(ind_time,ind_time+4,ind_time+8,ind_time+12,ind_time+16,ind_time+20),
#'        c(0,2.198356194,5.63889397,8.805379035,10.99968791,12.85308281),pch=3,col=4)
#' axis(4);
#' mtext("Recombinant protein(g/L)",side=4,cex=0.8);
#' legend("bottomright", c("Biomass","Recombinant protein"), col=c(2,4), lty=1, 
#'        pch=c(2,3),text.font = 1.2, cex = 0.6);

dGSMEF <- function (model,substrateRxns,initConcentrations,initBiomass,
						feedSubstrateRxns,feedConcentrations,initVolume = 0.1, initpH = 7.25,
						u_fix = 0.15,
						yield_rate = 0.6586, ## 25% glucose + 10%YE
						x_ind = 51.6,
						ngam = 11.36, ## Fc,maintenance energy for the recombinant cell,
						initRatio = 0.8,
						feedRate_ind = 2.9,   ## feeding rate during induction, ml/L
						timeStep,nSteps,exclUptakeRxns,
						rcd = FALSE,fld = FALSE,verboseMode = 2){


##########################################################################################################

# If no initial concentration is given for a substrate that has an open
# uptake in the model (i.e. model.lb < 0) the concentration is assumed to
# be high enough to not be limiting. If the uptake rate for a nutrient is
# calculated to exceed the maximum uptake rate for that nutrient specified
# in the model and the max uptake rate specified is > 0, the maximum uptake 
# rate specified in the model is used instead of the calculated uptake
# rate.

##--------------------------------------------------------------------------##
 # check prerequisites 
    if (!is(model, "modelorg")) {
      stop("needs an object of class modelorg!")
    }
##--------------------------------------------------------------------------##

# Uptake reactions whose substrate concentrations do not change
if (missing(exclUptakeRxns)){
    exclUptakeRxns = c('EX_co2(e)','EX_h2o(e)','EX_h(e)');
    if (verboseMode > 2){
       print('Default extra cellular uptake reactions will be used: ')
       print(exclUptakeRxns);
    }
}

# Find exchange reactions
excReact = findExchReact(model);
excReactInd=(react_id(model) %in% react_id(excReact));#excReact$exchange
#represent extra cellular reaction with boolean vector.
exclUptakeRxnsInd=is.element(react_id(model) ,exclUptakeRxns);
#Exclude reactions with concentrations that will not be changed 
excReactInd = excReactInd & !exclUptakeRxnsInd;   #excInd & ~ismember(model.rxns,exclUptakeRxns);
#get reaction names
RxnNames = react_id(model);
excRxnNames = RxnNames[excReactInd];                #excRxnNames = model.rxns(excInd);

substrateRxnsInd=(react_id(model) %in% substrateRxns)
# Figure out if substrate reactions are correct: all substrate reactions should be exchange reactions.
missingSub = substrateRxnsInd & !excReactInd;
if (sum(missingSub)!=0){
    print(sum(missingSub));
    print(react_id(model)[missingSub]);
    print('Invalid substrate uptake reaction!');
}
## 	***********************************************************     ##

### define bioreatcor model
### biomass: biomass concentration
### glc_feed_rates: glucose volumetric feed rates
### C_glc_feed: glucose concentration in the feeding solution: 500g/L = 2775.3mmol/L
### C_glc: glucose concentration in the bioreactor
### C_o2: dissolved oxygen centration in the bioreactor
### V: volumn of bioreactor
### D: dilution rate
### E: volumetric evaporation rate of water: 1.6L/h
### k: mass transfer coefficient for oxygen
### u: specific growth rate of the cellular
### 

GlucoseInd = which(feedSubstrateRxns == "EX_glc__D(e)");
C_glc_feed = feedConcentrations[GlucoseInd]; ### mMol/L
k = 7.5;   ### /h
E = 1.6;   ### L/h

V = initVolume;   ### L
pH = initpH;

# Initialize concentrations

protein = 0; ## initial recombinant protein cencentration

#substrateMatchInd = intersect(excRxnNames,substrateRxns);
concentrations=rep(0,length(react_id(model)))#table(excRxnNames);##vector(length=length(excRxnNames),mode="numeric");
#concentrations[1:length(concentrations)]=0;
#concentrations[substrateRxnsInd] = initConcentrations;
initconcentrationInd = c()
for (react in substrateRxns){
	if(any(react_id(model) == react)){ 
		initconcentrationInd = c(initconcentrationInd,which(react_id(model) == react))
	}
}
concentrations[initconcentrationInd] = initConcentrations;

# Deal with reactions for which there are no initial concentrations
originalBound = -lowbnd(model);# take all to be able to directly update
noInitConcentration = (concentrations==0)&(lowbnd(model)<0)#(concentrations == 0 & originalBound > 0);
concentrations[noInitConcentration] = 1000;

#biomass = initBiomass*V;
biomass = initBiomass; ## biomass concentration

# Initialize bounds
uptakeBound =  concentrations/(biomass*timeStep);

uptakeBound[uptakeBound > 1000] = 1000;
# Make sure bounds are not higher than what are specified in the model

aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
uptakeBound[aboveOriginal] = originalBound[aboveOriginal];
lowbnd(model)[excReactInd]  = -uptakeBound[excReactInd];

concentrationMatrix = c(concentrations[excReactInd],V,protein,pH);
biomassVec = initBiomass;
timeVec = 0;
fluxMatrix = rep(0,length(react_id(model)));
redCostMatrix = rep(0,length(react_id(model)));

##-----------------------------------------------------------------------------##
 if (verboseMode > 2) print('Step number    Biomass\n');
# Inititialize progress bar ...');
#if (verboseMode == 2)  progr <- .progressBar();

##-----------------------------------------------------------------------------##
# define ODE for each process variable under batch phase 

batch_bioreactor <- function(t,state,parameters){
	with(as.list(c(state, parameters)), {
		dbiomass <- u*biomass;
		res <- c(dbiomass);
		for(i in 1:length(excRxnNames)){
			if(excRxnNames[i] == "EX_o2(e)"){ dS <- k*(0.25-state[i+1]) + uptakeRate[i]*biomass; }
			else if(excRxnNames[i] == "EX_ac(e)"){ 
				dS <- uptakeRate[i]*biomass;		
				dpH = -10^(pH-10)*(5*dS + 25*exp(pH-7.25))/log(10,exp(1));
			}
			#else if(excRxnNames[i] == 'EX_na1(e)'){
			#	dS <- -84.7*dbiomass*state[i+1]^0.2/(1+16.94*biomass*state[i+1]^(-0.8));	
			#}
			else{ dS <- uptakeRate[i]*biomass; } 
			res <- c(res,dS);
		};
		
		list(c(res,dpH));
	})
}

# define ODE for each process variable under batch phase
fedbatch_bioreactor <- function(t,state,parameters){
	with(as.list(c(state, parameters)), {
		### define glucose feeding rates according to biomass concentration in the bioreactor
		glc_feed_rates <- 1/(C_glc_feed)*(-uptakeRate[excRxnNames == feedSubstrateRxns[GlucoseInd]])*biomass*V;
		#D <- (glc_feed_rates - E)/V
		D <- glc_feed_rates/V;
		dbiomass <- (u-D)*biomass;
		dV = D*V;
		dprotein <- q_p*biomass -D*protein;
		res <- c(dbiomass,dV,dprotein);
		feedSubstrateRxnsInd = which(excRxnNames %in% feedSubstrateRxns);
		for(i in 1:length(excRxnNames)){
			if(excRxnNames[i] == "EX_o2(e)"){ dS <- k*(0.25-state[i+3]) + uptakeRate[i]*biomass + D*(100-state[i+3]); }
			else if(excRxnNames[i] == "EX_ac(e)"){ 
				dS <- uptakeRate[i]*biomass - D*state[i+3];		
				dpH = -10^(pH-10)*(5*dS + 25*exp(pH-7.25))/log(10,exp(1));
			}
			#else if(excRxnNames[i] == 'EX_na1(e)'){
			#	dS <- -D*state[i+3] - 84.7*dbiomass*state[i+3]^0.2/(1+16.94*biomass*state[i+3]^(-0.8));	
			#}
			else if(any(feedSubstrateRxnsInd == i)){
				j = which(feedSubstrateRxns == excRxnNames[i]);
				dS <- uptakeRate[i]*biomass + D*(feedConcentrations[j] - state[i+3]);
			}
			else{ dS <- uptakeRate[i]*biomass - D*state[i+3]; }
			res <- c(res,dS);
		};

		res<-c(res,dpH);
		list(res);
	})
}

##-----------------------------------------------------------------------------##
# define AA, RNA and DNA compositions in biomass and required energy

growth_rate <- c()	
	
	## 0.027017 dgtp_c + 0.056843 trp__L_c + 2e-06 btn_c + 0.26316 glu__L_c + 0.423162 val__L_c + 0.09158 cys__L_c + 0.000223 ribflv_c + 
	## 0.153686 met__L_c + 48.601527 h2o_c + 0.000223 10fthf_c + 0.005205 cl_c + 0.000691 mn2_c + 0.144104 utp_c + 0.000223 thf_c + 
	## 0.094738 his__L_c + 0.295792 arg__L_c + 5.5e-05 udcpdp_c + 0.26316 gln__L_c + 2.5e-05 cobalt2_c + 0.450531 leu__L_c + 0.000223 pheme_c + 
	## 0.026166 datp_c + 0.133508 ctp_c + 0.241055 asn__L_c + 0.513689 ala__L_c + 0.185265 phe__L_c + 0.017868 pe160_c + 0.000323 ni2_c + 
	## 0.026166 dttp_c + 0.004338 so4_c + 0.013894 murein5px4p_p + 0.000122 bmocogdp_c + 2.6e-05 2fe2s_c + 7e-06 mobd_c + 0.221055 pro__L_c + 
	## 0.054154 pe161_c + 0.612638 gly_c + 0.000576 coa_c + 0.000223 pydx5p_c + 0.000223 fad_c + 0.008675 mg2_c + 0.241055 asp__L_c + 
	## 0.000223 sheme_c + 0.000223 amet_c + 0.253687 thr__L_c + 0.007808 fe3_c + 54.124831 atp_c + 0.000447 nadp_c + 0.013013 nh4_c + 
	## 0.000223 mlthf_c + 0.00026 4fe4s_c + 0.000709 cu2_c + 0.005205 ca2_c + 0.137896 tyr__L_c + 0.000223 thmpp_c + 0.001831 nad_c + 
	## 0.019456 kdo2lipid4_e + 0.045946 pe160_p + 0.000341 zn2_c + 0.215792 ser__L_c + 0.000223 2ohph_c + 0.290529 ile__L_c + 0.195193 k_c + 
	## 0.02106 pe161_p + 0.006715 fe2_c + 0.215096 gtp_c + 0.343161 lys__L_c + 0.027017 dctp_c --> 53.95 h_c + 53.945662 pi_c + 53.95 adp_c + 0.773903 ppi_c

core_biomass_coefficients <- c(-0.000223, -0.000026, -0.000223, -0.00026, 53.95, -0.513689, -0.000223, -0.295792, -0.241055, -0.241055, -54.124831, 
									-0.000122, -0.000002, -0.005205, -0.005205, -0.000576, -0.000025, -0.133508, -0.000709, -0.09158, -0.026166, -0.027017, 
									-0.027017, -0.026166, -0.000223, -0.006715, -0.007808, -0.26316, -0.26316, -0.612638, -0.215096, 53.95, -48.601527, 
									-0.094738, -0.290529, -0.195193, -0.450531, -0.343161, -0.153686, -0.008675, -0.000223, -0.000691, -0.000007, -0.001831, 
									-0.000447, -0.013013, -0.000323, -0.017868, -0.054154, -0.185265, -0.000223, 53.945662, 0.773903, -0.221055, -0.000223, 
									-0.000223, -0.215792, -0.000223, -0.004338, -0.000223, -0.000223, -0.253687, -0.056843, -0.137896, -0.000055, -0.144104, 
									-0.423162, -0.000341, -0.019456, -0.013894, -0.045946, -0.02106
									);
core_biomass_mets <- c('10fthf[c]', '2fe2s[c]', '2ohph[c]', '4fe4s[c]', 'adp[c]', 'ala__L[c]', 'amet[c]', 'arg__L[c]', 'asn__L[c]', 'asp__L[c]', 'atp[c]', 
							'bmocogdp[c]', 'btn[c]', 'ca2[c]', 'cl[c]', 'coa[c]', 'cobalt2[c]', 'ctp[c]', 'cu2[c]', 'cys__L[c]', 'datp[c]', 'dctp[c]', 'dgtp[c]', 
							'dttp[c]', 'fad[c]', 'fe2[c]', 'fe3[c]', 'gln__L[c]', 'glu__L[c]', 'gly[c]', 'gtp[c]', 'h[c]', 'h2o[c]', 'his__L[c]', 'ile__L[c]', 
							'k[c]', 'leu__L[c]', 'lys__L[c]', 'met__L[c]', 'mg2[c]', 'mlthf[c]', 'mn2[c]', 'mobd[c]', 'nad[c]', 'nadp[c]', 'nh4[c]', 'ni2[c]', 
							'pe160[c]', 'pe161[c]', 'phe__L[c]', 'pheme[c]', 'pi[c]', 'ppi[c]', 'pro__L[c]', 'pydx5p[c]', 'ribflv[c]', 'ser__L[c]', 'sheme[c]', 
							'so4[c]', 'thf[c]', 'thmpp[c]', 'thr__L[c]', 'trp__L[c]', 'tyr__L[c]', 'udcpdp[c]', 'utp[c]', 'val__L[c]', 'zn2[c]', 'kdo2lipid4[e]', 
							'murein5px4p[p]', 'pe160[p]', 'pe161[p]');

AA <- c('ala__L[c]','cys__L[c]','asp__L[c]','glu__L[c]','phe__L[c]','gly[c]','his__L[c]','ile__L[c]','lys__L[c]','leu__L[c]',
			'met__L[c]','asn__L[c]','pro__L[c]','gln__L[c]','arg__L[c]','ser__L[c]','thr__L[c]','val__L[c]','trp__L[c]','tyr__L[c]');	
RNA <- c('utp[c]','ctp[c]', 'gtp[c]', 'atp[c]');
DNA <- c('dgtp[c]', 'datp[c]', 'dttp[c]', 'dctp[c]');
	
aaInds <- core_biomass_mets %in% AA;
rnaInds <- core_biomass_mets %in% RNA;
dnaInds <- core_biomass_mets %in% DNA;

aa_composition = core_biomass_coefficients[aaInds];
## compute the RNA atp composition
core_biomass_coefficients[core_biomass_mets == 'atp[c]'] = core_biomass_coefficients[core_biomass_mets == 'atp[c]'] + 53.95;
rna_composition = core_biomass_coefficients[rnaInds];
core_biomass_coefficients[core_biomass_mets == 'atp[c]'] = core_biomass_coefficients[core_biomass_mets == 'atp[c]'] - 53.95;
	
dna_composition = core_biomass_coefficients[dnaInds];
	
total_aa = -sum(aa_composition);
total_rna = -sum(rna_composition);
total_dna = -sum(dna_composition);
	
macro_synthesis_atp = total_aa*4.306 +  total_rna*0.4 + total_dna*1.372;
other_growth_atp = 53.95 - macro_synthesis_atp;
	
non_growth_atp = lowbnd(model)[react_id(model) == 'ATPM'];
	
##-----------------------------------------------------------------------------##
#
# extract the recombinant protein reaction coefficients and reaction susbtractes

react_tmp <-  printReaction(model,react='RecombinantProtein');
equation <- unlist(strsplit(react_tmp,"\t"))[2];
equation_left <- unlist(strsplit(equation, "-->"));
metabolites <- unlist(strsplit(equation_left[1]," \\+ "));
metabolites <- unlist(strsplit(metabolites," "));
met_coeff <- metabolites[seq(1,length(metabolites),2)];
met_coeff <- gsub("[(]","",met_coeff);
met_coeff <- gsub("[)]","",met_coeff);
met_coeff <- as.numeric(met_coeff);
mets <- metabolites[seq(2,length(metabolites),2)];

# amino acids and coefficients
mets_aaInd <- mets %in% AA;
mets_aa_coeff <- met_coeff[mets_aaInd];

#mass_aa <- c(71.0788,156.1875,114.1038,115.0886,103.1388,129.1155,128.1307,57.0519,137.1411,113.1594,113.1594,128.1741,131.1926,147.1766,97.1167,87.0782,101.1051,186.2132,163.176,99.1326);
mass_aa <- c(89.09,174.2,132.12,133.1,121.16,146.15,147.13,75.07,155.16,131.17,131.17,146.19,149.21,165.19,115.13,105.09,119.12,204.23,181.19,117.15);

##-----------------------------------------------------------------------------##
#
# define a flag, indicating whether glucose is consumed out in the batch phase
C_glc_flag = 0;
i = 0; ### indicate feed_batch phase is coming after the transition from batch phase
for (stepNo in 1:nSteps){
	
	##-----------------------------------------------------------------------------##
	## biomass composition changes along with growth rate
	## biomass content of a given molecular species  = slope * mu +  intercept

	## compute the average smoothed growth rate over the previous 2 h
	mean_mu = 0;
	if(length(growth_rate) == 0 ) { 
		sol = sybil::optimizeProb(model,algorithm="fba",retOptSol=TRUE);
		if ( length(checkSolStat(lp_stat(sol),solver(sol)))!=0 ){## checkSolStat
			print('No feasible solution - nutrients exhausted\n');
			break;
		}
		objective = lp_obj(sol);  ##objvalue sol.f
		#print(paste("initial growth rate: ",objective));
		mean_mu = objective; 
	}
	else if (length(growth_rate) >= 2/timeStep) { mean_mu = mean(growth_rate[ (length(growth_rate)-2/timeStep + 1):length(growth_rate) ]  );  }
	else{ mean_mu = mean(growth_rate); }
	
	## for protein
	protein_quantity <- -1.76*mean_mu + 7.09;
	
	## for RNA
	rna_quantity <- 0.277*mean_mu + 0.312;

	## for DNA
	dna_quantity <- -0.0853*mean_mu + 0.186;
	
	## glycogen: only in BIOMASS_Ec_iJO1366_WT_53p95M, and assumed constant in all conditions 
	## cell wall: assumed constant
	
	core_biomass_coefficients[aaInds] <- protein_quantity/total_aa*aa_composition;	
	core_biomass_coefficients[rnaInds] <- rna_quantity/total_rna*rna_composition;
	core_biomass_coefficients[dnaInds] <- dna_quantity/total_dna*dna_composition;
	
	## compute the coefficient of atp[c]: engergy ATP needed and RNA nucleotide A
	core_biomass_coefficients[core_biomass_mets == 'atp[c]'] = core_biomass_coefficients[core_biomass_mets == 'atp[c]']-(protein_quantity*4.306 + rna_quantity*0.4 + dna_quantity*1.372 + other_growth_atp);
	#print(core_biomass_coefficients[core_biomass_mets == 'atp[c]']);
	
	## compute the coefficient of atp[c]'s products: only engergy ATP
	core_biomass_coefficients[core_biomass_mets == 'adp[c]'] = protein_quantity*4.306 + rna_quantity*0.4 + dna_quantity*1.372 + other_growth_atp;
	core_biomass_coefficients[core_biomass_mets == 'h[c]'] = protein_quantity*4.306 + rna_quantity*0.4 + dna_quantity*1.372 + other_growth_atp;
	core_biomass_coefficients[core_biomass_mets == 'pi[c]'] = protein_quantity*4.306 + rna_quantity*0.4 + dna_quantity*1.372 + other_growth_atp;
		
	model <- addReact(model,id="BIOMASS_Ec_iJO1366_core_53p95M",
				met = core_biomass_mets,
				Scoef = core_biomass_coefficients,
				lb=0,ub=1000,obj=1);

	# run FBA
	
	## inducing protein expression when biomass grows to the predefined cut-off
	if(biomass >= x_ind) {
		
		## non-growth rate dependent energy requirement is changed from the host to the recombinant cell
		lowbnd(model)[react_id(model) == 'ATPM'] = ngam; ## 32.075 for MIC-1
		#lowbnd(model)[react_id(model) == 'ATPM'] = 11.36; ## for Fc
		
		## proteome resource allocation model
		# Run FBA first to get growth rate without recombinant proteins
		model = changeObjFunc(model,c('BIOMASS_Ec_iJO1366_core_53p95M','RecombinantProtein'),c(1,0));
		sol = sybil::optimizeProb(model,algorithm="fba",retOptSol=TRUE);
		u_max = lp_obj(sol);
		
		# then u/u_max = 1 - mass fraction of recombinant protein/0.48
		recomb_pro_mass = protein*(sum(mets_aa_coeff*mass_aa) - (sum(mets_aa_coeff) - 1)*18.01524)/1000;
		biomass_pro_mass = biomass*(abs(sum(core_biomass_coefficients[aaInds]*mass_aa)))/1000;
		proteome_mass = biomass_pro_mass + recomb_pro_mass;
		fixed_pro = biomass_pro_mass*0.48/proteome_mass;
		fixed_pro = 0.48;
		u = u_max*(1 - recomb_pro_mass/proteome_mass/fixed_pro);
		if(u > initRatio*u_max){ u = initRatio*u_max; }
		
		if(u <= 0){
			return('Cells stop to growth -- recombinant protein proportion exceeds the limit\n');
		}
		
		# set the growth rate to u
		lowbnd(model)[react_id(model) == 'BIOMASS_Ec_iJO1366_core_53p95M'] = u;
		uppbnd(model)[react_id(model) == 'BIOMASS_Ec_iJO1366_core_53p95M'] = u;
		
		# change the objective to recombinant protein synthesis
		model = changeObjFunc(model,c('BIOMASS_Ec_iJO1366_core_53p95M','RecombinantProtein'),c(0,1));
	}
	
	
	sol = sybil::optimizeProb(model,algorithm="fba",retOptSol=TRUE,poCmd = list("getRedCosts"));

	if ( length(checkSolStat(lp_stat(sol),solver(sol)))!=0 ){## checkSolStat
        print('No feasible solution - nutrients exhausted\n');
        break;
    }
	
	## get reduced costs for each reaction
	if(rcd){
		ppProc_obj = postProc(sol);
		redCosts = pa(ppProc_obj)[[1]];
		redCostMatrix <- c(redCostMatrix,redCosts);
	}
	
	# minimize total flux
	mtf = sybil::optimizeProb(model,algorithm="mtf",wtobj=lp_obj(sol));
	fluxes = getFluxDist(mtf);
	#fluxes = getFluxDist(sol);
	
	# protein production rate
	q_p = fluxes[react_id(model) == 'RecombinantProtein'];
	print(paste("protein synthesis rate: ",q_p))
	# growth rate
	mu = fluxes[react_id(model) == 'BIOMASS_Ec_iJO1366_core_53p95M'];
	print(paste("growth rate: ",mu))
	
	#mu2 = fluxes2[react_id(model) == 'BIOMASS_Ec_iJO1366_core_53p95M'];
	#print(paste("growth rate2: ",mu2))
	
	growth_rate <- c(growth_rate,mu);
	
	uptakeFlux = fluxes[excReactInd];
	if(fld) {
		fluxMatrix <- c(fluxMatrix,fluxes)
	}

	##-----------------------------------------------------------------------------##
	# Update concentrations
	# Note: uptakeFlux: a vector of fluxes of external metabolites; 
	# 					excRxnNames and concentrations[excReactInd] has the same structure or dismension as uptakeFlux
	# 		concentrations: a vector of all metabolites, 
	#						although only external metabolites' concentrations would change along with the time
	# 
	
	### batch phase

	if(C_glc_flag == 0){
		parameters <- list(u = mu,excRxnNames = excRxnNames,
						uptakeRate = uptakeFlux
						);
		state <- c( biomass = biomass,
					S = concentrations[excReactInd],
					pH = pH
					);
		out <- ode(y=state,times=c(0,timeStep),func=batch_bioreactor,parms=parameters);

		biomass = as.numeric(out[2,2]);
		concentrations[excReactInd]  = as.numeric(out[2,3:(length(out[2,]) - 1)]);
		pH = 7.25;

		### Set up a signal meaning batch ending: glucose is consumed out and also oxygen consumption is reduced
		if(concentrations[RxnNames == feedSubstrateRxns[GlucoseInd]] <= 0 && (uptakeFlux[excRxnNames == "EX_o2(e)"] > lowbnd(model)[RxnNames =='EX_o2(e)'] || concentrations[RxnNames == 'EX_ac(e)'] <= 0)){ 
			C_glc_flag = 1;
			concentrations[RxnNames == feedSubstrateRxns[GlucoseInd]] = 0;
			pH = as.numeric(out[2,length(out[2,])]); ## pH rises due to the consumption of acetate
		} 
	}
	
	### when glucose is consumed out, entering into fed-batch phase
	if(C_glc_flag == 1){
		if(i > 0){
			parameters <- list(u = mu,q_p = q_p,excRxnNames = excRxnNames,
							feedSubstrateRxns = feedSubstrateRxns,
							feedConcentrations = feedConcentrations,
							#E = E,
							uptakeRate = uptakeFlux
							);
			state <- c( biomass = biomass,V = V, protein = protein,
						S = concentrations[excReactInd],
						pH = pH
						);
			out <- ode(y=state,times=c(0,timeStep),func=fedbatch_bioreactor,parms=parameters);
	
			#print(out);
			biomass = as.numeric(out[2,2]);
			V = as.numeric(out[2,3]);
			protein = as.numeric(out[2,4]);
			concentrations[excReactInd]  = as.numeric( out[2,5:(length(out[2,])-1)] );
			if(pH > 7.25){
				pH = as.numeric(out[2,length(out[2,])]);
			}
			else{ pH = 7.25; }
		}
		i = i + 1;
	}
	
    biomassVec = c(biomassVec,biomass);

    concentrations[concentrations <= 0] = 0;
	concentrations=ifelse(abs(concentrations) < 1e-9,0,concentrations);
    concentrationMatrix = c(concentrationMatrix,concentrations[excReactInd],V,protein,pH);
    
    # Update bounds for uptake reactions
    uptakeBound[excReactInd] =  concentrations[excReactInd]/(biomass*timeStep);
	### feeding glucose: glc_feed_rates * C_glc_feed = -uptakeFlux[react_id(model)[excReactInd] == "EX_glc(e)"]*biomass*V
	if(C_glc_flag == 1) {
		### for the first time of glucose feeding, set a fixed growth rate: 0.15/h, Yx/s: 0.5g DCW/g glucose
		### glc_feed_rates * C_glc_feed = u/Yx/s*biomass*V
		### and also transfer unit g to mol: 180.16g/mol for glucose
		
		if(i == 1) {
			if(!any(react_id(model) == "PGL")){
				uptakeBound[react_id(model) == feedSubstrateRxns[1] ] = u_fix/0.418/180.16*1000*biomass*timeStep/(biomass*timeStep);
			}
			else{
				#uptakeBound[react_id(model) == feedSubstrateRxns[1] ] = u_fix/0.547/180.16*1000*biomass*timeStep/(biomass*timeStep); #HNC47
				#uptakeBound[react_id(model) == feedSubstrateRxns[1] ] = u_fix/0.4413/180.16*1000*biomass*timeStep/(biomass*timeStep); # HNC47 with Fc
				uptakeBound[react_id(model) == feedSubstrateRxns[1] ] = u_fix/yield_rate/180.16*1000*biomass*timeStep/(biomass*timeStep); # HNC47 with Fc with 25% glucose + 10% YE
			}
			init_feed_rate = uptakeBound[react_id(model) == feedSubstrateRxns[1] ]*biomass*V/feedConcentrations[GlucoseInd]*1000;
			print(paste("Initial feed rate: ", init_feed_rate," ml/h"));
		}
		### init feeding rate: 0.36 ml/h
		#if(i == 1) {uptakeBound[react_id(model) == feedSubstrateRxns[GlucoseInd] ] = 0.36/1000*feedConcentrations[GlucoseInd]/(biomass*V); }
		
		## inducing protein expression when biomass grows to the predefined cut-off
		## feeding rate: 2.9ml/h
		else if (biomass >= x_ind){
			print(paste("Feeding rate at induction: ",feedRate_ind))
			uptakeBound[react_id(model) == feedSubstrateRxns[GlucoseInd]] = feedRate_ind/1000*feedConcentrations[GlucoseInd]/(biomass*V);
		}
		## set to a fixed growth rate by setting the uptake rate
		else{ uptakeBound[react_id(model) == feedSubstrateRxns[GlucoseInd]] = -uptakeFlux[excRxnNames == feedSubstrateRxns[GlucoseInd]]*biomass*timeStep/(biomass*timeStep); }
	}
	
	##print(paste("Glucose uptake rate: ",uptakeBound[react_id(model) == feedSubstrateRxns[GlucoseInd]],sep=""))
	
    # This is to avoid any numerical issues
    uptakeBound[uptakeBound > 1000] = 1000;
	
	### suppose o2 is supplied enough and other N,S,P are also supplied enough
	### noInitConcentration: reactions for which there are no initial concentrations but required for cell's growth
	### e.g., EX_pi(e) is needed for cells, and each run some of pi will be consumed; Many runs later, N won't be enough to support cell's continuing growth,
	### although at first the concentration of EX_nh4(e) is set to 1000.
	uptakeBound[react_id(model) == "EX_o2(e)"] = 15;
	#uptakeBound[noInitConcentration] = 1000;
	
    # Figure out if the computed bounds were above the original bounds
    aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
	
    # Revert to original bounds if the rate was too high
    uptakeBound[aboveOriginal] = originalBound[aboveOriginal];# uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
    uptakeBound=ifelse(abs(uptakeBound) < 1e-9,0,uptakeBound);
	
    ## Change lower bounds according to the result of last step
    lowbnd(model)[excReactInd]  = -uptakeBound[excReactInd];  
    
	#uppb_tmp <- getColsUppBnds(problem(lpmod), which(excReactInd));
    #changeColsBnds(problem(lpmod),which(excReactInd),lb=-uptakeBound[excReactInd],ub=uppb_tmp);
    
    if (verboseMode > 2) print(paste(stepNo,sep="    ",biomass));
    #waitbar(stepNo/nSteps,h);
    timeVec = c(timeVec,stepNo*timeStep);
}# end loop

## Preparing OUTPUT
#concentrationMatrix,excRxnNames,timeVec,biomassVec
return (optsol_dGSMEF(
		      concmat=concentrationMatrix,
		      exRxn=c(excRxnNames,"V","Protein","pH"),
              tmVec=timeVec,  
              bmVec=biomassVec,
			  flux=fluxMatrix,
			  redCost=redCostMatrix,
			  Rxn=RxnNames
        )
  )
}
