// Stata code for paper "Rotten, stale, fresh: three flavors of heritability"
// by Alex Gamma & Michael Liebrenz

// Code by Alex Gamma, 30.10.2017

// We assess the interactionist challenge to behavioral genetics by
// simulating additive and multiplicative regression models 
// for interacting causes in development

// We assume x1/x2 are the levels of two biomolecules
// that are both needed to produce some product, 
// whose levels are denoted y<something>

// Both an additive and a multiplicative relationship between
// y and x1/x2 are tested



// Clear all observations
clear

// Set random seed
set seed 31

// Parameters
local N  = 1000		// Number of simulated observations
local mu = 15 		// mean of x1/x2
local sd = 8 		// standard deviation of x1/x2

// Generate 'space' for the simulated observations
set obs `N'


// Simulate a range of correlations betw x1 and x2:
// r = 0 .25 .5 .75 .9
// Run analysis for each correlation

// Counter for the different analysis runs
local n = 0

// Run loop for each correlation
foreach c of numlist 0 .25 .5 .75 .9 {

	local `++n'
	// Specify desired correlation matrix
	matrix C = (1, `c' \ `c', 1)
	// Draw observations from the normal distribution
	drawnorm x1_`n' x2_`n', means(`mu', `mu') sds(`sd', `sd') corr(C)
	
	
	// Truncate x1/x2 at 0 (negative values not possible)
	replace x1_`n' = 0 if x1_`n' < 0
	replace x2_`n' = 0 if x2_`n' < 0
	
	// Reassess correlations
	qui cor x1_`n' x2_`n'
	local corr = round(`r(rho)', .01)
		
	// Generate additive relationship betw y and x	
	gen yadd_`n' = x1_`n' + x2_`n'
	
	// Encode interactionist constraint: set y = 0 if x1 or x2 are 0
	replace yadd_`n' = 0 if x1_`n' * x2_`n' == 0
	
	// Generate multiplicative relationship betw y and x
	gen ymult_`n' = x1_`n' * x2_`n'
		
	
	// Regression analyses
	// Each y is tested with an additive (main effects only) model
	// and an interaction model (main effects + interaction term)
	
	// yadd: additive model
	reg yadd_`n' x1_`n' x2_`n'
	// Store regression estimates for later plotting
	est store yadd_add_`n'
	
	// yadd: interaction model
	reg yadd_`n' c.x1_`n'##c.x2_`n'
	est store yadd_i_`n'
	
	// ymult: additive model
	reg ymult_`n' x1_`n' x2_`n'
	est store ymult_add_`n'
	
	// ymult: interaction model
	reg ymult_`n' c.x1_`n'##c.x2_`n'
	est store ymult_i_`n'
	
	// Plot coefficients of the 4 models
	// If necessary, install package "coefplot" by typing
	// "ssc install coefplot" on Stata's command line
	coefplot ///
		(yadd_add_`n', label("Additive outcome, additive model") offset(-.1)) 		///
		(yadd_i_`n', label("Additive outcome, interaction model") offset(0)) 		/// 
		(ymult_add_`n', label("Multiplicative outcome, additive model") offset(0)) 	///
		(ymult_i_`n', label("Multiplicative outcome, interaction model") offset(0)) ///
		, drop(_cons) ti("Correlation of x1 with x2 = `corr'", size(medium)) 		///
		  coeflabels( 																///
			x1_`n' 		  = "x1" 													///
			x2_`n' 		  = "x2" 													///
			c.x1_`n'#c.x2_`n' = "Interaction x1 * x2" 								///
		  ) 																		///
		  xlab(0 1 5 10 15) xti(coefficient size)									///
		  legend(cols(1) pos(11))													///
		  name(coefplot_4models_`n', replace)

}



*LF
