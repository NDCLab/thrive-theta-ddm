#include <Rcpp.h>
#include <Rcpp/Rmath.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix simSSP_model_GBnoScale(NumericVector parms, int trialType, int nTrials,
double dt, double vari) {
  
  // This is an rcpp function that simulates trial-level results for the SSP-DDM mode fit, 
  // given a set of input parameters. This is implemented in rcpp to speed up the 
  // simulations (and fits, respectively). Inputs needed to run:
  
  // parms: numeric vector (length 5) corresponding to the SSP params
  // trialType: integer; 1 = congruent, 2 = incongruent (what kind of trials to simulate)
  // nTrials: integer; how many trials to simulate
  // dt: double; scaler diffusion process (i.e. size of each step)
  // var: double; drift variance for each step

  //////////////////////////////////////////////////////////////////////////////
  // set fixed parameters
  double sdRand = sqrt(dt * vari); //get standard deviation from variance

  // get the free parameters from the numeric vector parms
  double A;
    A = parms[0]; // boundary separation = response A boundary
  double B;
    B = 0.0; // response B boundary (essentially 0)
  double tEr;
    tEr = parms[1]; // non-decision time
  double p;
    p = parms[2]; // perceptual strength
  double rd;
    rd = parms[3]; // shrinking rate
  double sda;
    sda = parms[4]; // attentional window width

    // declare random number variable for noise
    double noise;
      noise = 0.0;

  // set empty matrix to store trial data
  int nRow = nTrials; //first set number of rows
  NumericMatrix trialData(nRow, 2); //set the matrix. 2 columns (RT & accuracy)

  // keep track of current evidence (response selection)
  double currEvidenceResp;
    currEvidenceResp = 0.0;
  double t; //to log number of diffusion steps
    t = 0.0;
  double drift; //the drift rate for the current moment
    drift = 0.0;

  // declare trial-specific parameters
  double sd_t; //current sd of spotlight
    sd_t = 0.0;
  double a_target;
    a_target = 0.0;
  double a_flanker;
    a_flanker = 0.0;
  double p_target;
    p_target = 0.0;
  double p_flanker;
    p_flanker = 0.0;
  double mu_target;
    mu_target = 0.0;
  double mu_flanker;
    mu_flanker = 0.0;
  //////////////////////////////////////////////////////////////////////////////

  // start trial loop here

  for (int i=0; i<=nTrials - 1; i++){

    // keep track of the current evidence (response selection)
    currEvidenceResp = A/2; //starting point is halfway between 0 and boundary

    // reset counter to log how many steps taken in diffusion process
    t = 0.0;

      //////////////////////////////////////
      // diffusion simulation starts here

      while((currEvidenceResp <= A) && (currEvidenceResp >= B)){
        
        t = t + 1.0; // update diffusion step number
        
        // calculate current sd of spotlight
        sd_t = sda - (rd * t);
        if(sd_t <= 0.001){ // clip min sd to 0.001 // why limiting 0.001
          sd_t = 0.001;
        }

      // find area of spotlight over target and flanker
      // :Rf_pnorm5(value, mean of dist., sd of dist., 1 = p<x, 0 = p and 1=log(p))
      a_target = ::Rf_pnorm5(0.5, 0.0, sd_t, 1, 0) - ::Rf_pnorm5(-0.5, 0.0, sd_t, 1, 0);
      a_flanker = ::Rf_pnorm5(10.0, 0.0, sd_t, 1, 0) - ::Rf_pnorm5(0.5, 0.0, sd_t, 1, 0); // why 10
  //      a_flanker = 2 * a_flanker;

      // get perceptual input for targets and flanker
      p_target = p;
      p_flanker = p;

        // flip the sign if current trial is incongruent
        if(trialType==2){
          p_flanker = -p_flanker;
        }

      // current evidence
      mu_target = (p_target * a_target) * dt;
      mu_flanker = ((2 * p_flanker) * a_flanker) * dt; // why multiplying by 2

      // current drift rate
      drift = mu_flanker + mu_target;
        // drift = drift * dt;

        // get random noise
        RNGScope scope;
        noise = ::Rf_rnorm(drift, sdRand);

        // Rcout << "trial : " << i << "\n";
        // Rcout << "tstep : " << t << "\n";
        // Rcout << "sd_t : " << sd_t << "\n";
        // Rcout << "drift : " << drift << "\n";
        // Rcout << "sdRand : " << sdRand << "\n";
        // Rcout << "noise : " << noise << "\n";
        // Rcout << "currEvience : " << currEvidenceResp << "\n";

      //update the response selection random walk
      currEvidenceResp = currEvidenceResp + noise;
      
      } // while loop ends here

   trialData(i, 0) = (t * dt) + tEr; // populate simulated RT

    if(currEvidenceResp >= A){
      trialData(i, 1) = 1; // populate simulated accuracy
    }
      if(currEvidenceResp <= B){
        trialData(i, 1) = 0;
      }
  } // trial loop ends here
     return trialData;
}