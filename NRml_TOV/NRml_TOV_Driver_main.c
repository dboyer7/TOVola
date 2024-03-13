#include "NRml_TOV_funcs.c" //nrpy_odiegm itself, for use in NRml_TOV.
#include "NRml_TOV_user_methods.c" //Edited specifically for the TOV solver.
#include "../GRHayL/GRHayL/include/nrpyeos_tabulated.h" //Check that this filepath is correct.
#include <ghl.h> //Need the GRHayL Library

/********************************
// UPDATE: NRml_TOV is a special version of the Odie main driver. It is built to handle the TOV equations specifically, and is built to handle a variety of different types of EOS
//
// Simple Polytrope:
//	-- The Simple Polytrope is a special version of the Piecewise Polytrope with one EOS region relating Pressure to baryon density (P=K*rho^Gamma), where K and Gamma are parameters to be input
//      -- Only one K and Gamma are necessary to run this implementation, and a Gamma_Thermal defaulted to 2.0
//      -- Tested against the TOV solver from NRPy, and got agreement around 10^-13 in relative error.
// Piecewise Polytrope:
//	-- The more general verion of the simple polytrope where you can choose the number of regions for your polytrope, each region denoted by (P=K_i*rho^Gamma_i)
//      -- You can get your necessary parameters by following the method outlined in Read et. al. (2008)
//      -- Since we follow Read's method, you only need to input one K value (The GRHayL library will handle the rest of them), and all the Gammas and rho boundary points for the number of specified regions.
//      -- Test against the TOV solver from NRPy, and got agreement around 10^-13 in relative error.
// Tabulated EOS:
//	-- The tabulated solver reads in a EOS table from your computer.
//      -- Once located, GRHayL slices the table for beta_equilibirum and then uses an interpollator to find the necessary densities through the calculations.
//      -- Tested against the TOV solver from NRPy, and got agreement around 10^-8 in relative error (Adaptive methods do cause some eventual interpolation error down the road)
//
// ALL EDITABLE VARIABLES ARE WITHIN THE FIRST CHUNK OF LINES OF THIS FILE!
// 40 Lines of white space separates what you should and should not edit.
//
// Read the documentation for more help in using this program.
//
// David Boyer (03/13/24)
********************************/
int main(){

  //COMBINED TEMPLATE: All variables present

    printf("Beginning ODE Solver \"Odie\" for NRml_TOV...\n");

    //Choose your EOS type and (if necessary) path
    ghl_eos_parameters ghl_eos={};
    char type = 't'; //s=simple, p=piecewise, t=tabulated
    char* table_path = "/home/boye4060/projects/New_Odie_TOV/NRml_TOV/EOS_Tables/SLy4_3335_rho391_temp163_ye66_adjusted.h5"; //Replace with your own path to your table.

    //These values are for initializing the Simple and Piecewise Polytrope
    int neos = 7;
    double* Gamma = (double *)malloc(neos * sizeof(double));
    double* rhoBound = (double *)malloc(neos * sizeof(double));
    double K0 = 168.57487497864867;
    double GammaTh = 2.0;
    Gamma[0] = 1.58425;
    Gamma[1] = 1.28733;
    Gamma[2] = 0.62223;
    Gamma[3] = 1.35692;
    Gamma[4] = 3.005;
    Gamma[5] = 2.988;
    Gamma[6] = 2.851;
    //...
    //Gamma[neos-1] = final region
    //Add or takeaway Gammas as necessary
    rhoBound[0] = 3.95143746008251e-11;
    rhoBound[1] = 6.126433097526977e-07;
    rhoBound[2] = 4.254975682734708e-06;
    rhoBound[3] = 0.0002367785968890901;
    rhoBound[4] = 0.0008115303644041096;
    rhoBound[5] = 0.001619215953548485;
    //...
    //rhoBound[neos-2] = final region
    //Add or takeaway rhoBounds as necessary
    
    //These values are for initializing a tabulated EOS:
    //Designed to enforce the table limits.                                                                             
    double YeAtm = 1.0e-2;
    double YeMin = -1.0;
    double YeMax = 1.0;
    double TAtm = 1.0e-2;
    double TMin = -1.0;
    double TMax = 1.0e10;
    double T_in = TAtm;

    //These values are for initializing all types of EOS.
    double rhoCentral = 1.58e-3;
    double rhoAtm = 1.0e-13;
    double rhoMin = -1.0;
    double rhoMax = -1.0;

    //Values for ODE steps and Error
    double step = 1e-5;
    double current_position = 0.0;
    const int size = 100000000;
    int adams_bashforth_order = 4; //Only edit if you are using an AB method (NOT RECCOMMENDED)
    bool adaptive_step = true;
    double error_limit = 1e-16;
    
    char file_name[] = "ooData.txt"; 

    //Choose your desired ODE method
    //Recommended: (ARKF: (RK4(5))) or (ADP8: (DP7(8)))
    //Other Methods listed in Documentation
    const nrpy_odiegm_step_type * step_type;
    step_type = nrpy_odiegm_step_ARKF;

    //additional tolerance values you can edit.
    //Only mess with these if you know what you are doing.
    //For more info, see OdieGM and the GSL ODE solver.
    double scale_factor = 0.9;
    double error_safety = 4.0/15.0;
    double ay_error_scaler = 1.0;
    double ady_error_scaler = 1.0;
    double max_step_adjustment = 3.0;
    double min_step_adjustment = 0.1;
    double absolute_max_step = 1e-1;
    double absolute_min_step = 1e-10;
    double error_upper_tolerance = 1.0e0;
    double error_lower_tolerance = 5.0e-1;

    // AFTER THIS POINT THERE SHOULD BE NO NEED FOR USER INPUT, THE CODE SHOULD HANDLE ITSELF.
    



    
























    

    




    



    //Initializing the EOS
    if(type=='s'){
      printf("Simple Polytrope\n");
      if(neos!=1){
        printf("Error: Too many regions for the simple polytrope.\n");
        printf("Check your value for neos, or use a piecewise polytrope (type = 'p')\n");
        printf("Shutting down due to error...\n");
        return 1;}
      ghl_initialize_hybrid_eos_functions_and_params(rhoAtm,rhoMin,rhoMax,neos,rhoBound,Gamma,K0,GammaTh,&ghl_eos);}
    else if(type=='p'){
      printf("Piecewise Polytrope\n");
      ghl_initialize_hybrid_eos_functions_and_params(rhoAtm,rhoMin,rhoMax,neos,rhoBound,Gamma,K0,GammaTh,&ghl_eos);}
    else if(type=='t'){
      printf("Tabulated EOS\n");
      ghl_initialize_tabulated_eos_functions_and_params(table_path,rhoAtm,rhoMin,rhoMax,YeAtm,YeMin,YeMax,TAtm,TMin,TMax,&ghl_eos);
      ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(T_in, &ghl_eos);
    }
    else{
      printf("ERROR: Invalid EOS type. Must be either 's' (simple), 'p' (piecewise), or 't' (tabulated)\n");
      printf("Shutting down due to error...\n");
      return 1;}


    const nrpy_odiegm_step_type * step_type_2;
    step_type_2 = step_type;

    //A few error checks for the ODE steps.
    if(step_type != step_type_2){
    	printf("ERROR: Hybridizing method disabled for the TOV solver.\n");
	printf("Shutting down due to error.\n");
    	return 1;}

    if(adaptive_step==true){
      if(step_type != nrpy_odiegm_step_AHE
	 && step_type != nrpy_odiegm_step_ABS
	 && step_type != nrpy_odiegm_step_ARKF
	 && step_type != nrpy_odiegm_step_ACK
	 && step_type != nrpy_odiegm_step_ADP5
	 && step_type != nrpy_odiegm_step_ADP8){
	printf("ERROR: Incompatible step type with adaptive method.\n");
	printf("Shutting down due to error.\n");
	return 1;}
    }

    if(adaptive_step==false){
      if(step_type == nrpy_odiegm_step_AHE
         || step_type == nrpy_odiegm_step_ABS
         || step_type == nrpy_odiegm_step_ARKF
         || step_type == nrpy_odiegm_step_ACK
         || step_type == nrpy_odiegm_step_ADP5
         || step_type == nrpy_odiegm_step_ADP8){
        printf("ERROR: Incompatible step type with non-adaptive method.\n");
        printf("Shutting down due to error.\n");
        return 1;}
    }

    //GSL requires the use of the variable 'no_adaptive step'
    //However, double negatives are confusing for user interface.
    //We use the postive variable inputted by the user to declare the GSL compatible variable here.
    bool no_adaptive_step;
    if (adaptive_step == true){
      no_adaptive_step=false;}
    else{
      no_adaptive_step=true;}

    //Define an ODE struct
    double absolute_error_limit = error_limit;
    double relative_error_limit = error_limit;
    int number_of_equations = 4;
    int number_of_constants = 2;
    struct constant_parameters cp = { 
      .dimension = number_of_constants,
      .neos = neos,
      .Gamma=(double *)malloc(neos * sizeof(double)),
      .rhoBound=(double *)malloc(neos * sizeof(double)),
    };
    cp.Gamma = Gamma;
    cp.rhoBound = rhoBound;
    cp.type = type;
    cp.rhoCentral = rhoCentral;
    cp.ghl_eos = ghl_eos;
    

    nrpy_odiegm_system system = {diffy_Q_eval,known_Q_eval,number_of_equations,&cp};
    

    nrpy_odiegm_driver *d;
    d = nrpy_odiegm_driver_alloc_y_new(&system, step_type, step, absolute_error_limit, relative_error_limit); 

    //setting the control struct. See OdieGM and the GSL ODE solver for more info.
    d->c->scale_factor = scale_factor;
    d->c->error_safety = error_safety;
    d->c->ay_error_scaler = ay_error_scaler;
    d->c->ady_error_scaler = ady_error_scaler;
    d->c->max_step_adjustment = max_step_adjustment;
    d->c->min_step_adjustment = min_step_adjustment;
    d->c->absolute_max_step = absolute_max_step;
    d->c->absolute_min_step = absolute_min_step;
    d->c->error_upper_tolerance = error_upper_tolerance;
    d->c->error_lower_tolerance = error_lower_tolerance;

    //Some OdieGM checks for AB methods
    int method_type = 1;
    if (step_type->rows == step_type->columns) {
        method_type = 0;  
    } 
    if (step_type->rows == 19) { 
        method_type = 2;
    } else {
        adams_bashforth_order = 0;
    }
    d->s->adams_bashforth_order = adams_bashforth_order;
    d->e->no_adaptive_step = no_adaptive_step;
    

    if (method_type == 2) {
        printf("Method Order: %i.\n",adams_bashforth_order);
    } else {
        printf("Method Order: %i.\n",step_type->order);            
    }

    //Declare ODEs and begin SOLVING!
    double y[number_of_equations];
    

    double c[number_of_constants];
    
    get_initial_condition(y,&cp); 
    //const_eval(current_position, y,&cp);
    assign_constants(c,&cp); 

    FILE *fp2;
    fp2 = fopen(file_name,"w");
    printf("Printing to file '%s'.\n",file_name);

    
    // Print the order: r, e(rho), rho_b, P, M, nu, 2*nu, Iso_r
    printf("INITIAL: Position:,\t%f,\t",current_position);

    fprintf(fp2, "%15.14e",current_position); //r
    fprintf(fp2, " %15.14e", c[0]); //e(or rho)
    fprintf(fp2, " %15.14e", c[1]); //rho_b
    fprintf(fp2, " %15.14e", y[0]); //P
    fprintf(fp2, " %15.14e", y[2]); //M
    fprintf(fp2, " %15.14e", exp(y[1])); //nu
    fprintf(fp2, " %15.14e", exp(2*y[1])); //2*nu
    fprintf(fp2, " %15.14e", y[3]); //Iso_r
    for (int n = 0; n < number_of_equations; n++) {
        printf("Equation %i:,\t%15.14e,\t",n, y[n]);
    }
    for (int n = 0; n < number_of_constants; n++) {
        printf("Constant %i:,\t%15.14e,\t",n, c[n]);
    }
    printf("\n");
    fprintf(fp2,"\n");

    for (int i = 0; i < size; i++){
        
        if (method_type == 2 && i == 0 && step_type_2 != nrpy_odiegm_step_AB) {
            d->s->type = step_type_2;
            d->s->rows = step_type_2->rows;
            d->s->columns = step_type_2->columns;
            d->s->method_type = 0;
            d->s->adams_bashforth_order = adams_bashforth_order;
            d->e->no_adaptive_step = true;
        } else if (step_type != step_type_2 && method_type == 2 && i == adams_bashforth_order) {
            d->s->type = step_type;
            d->s->rows = step_type->rows;
            d->s->columns = step_type->columns;
            d->s->method_type = 2;
            d->s->adams_bashforth_order = adams_bashforth_order;
            d->e->no_adaptive_step = true;
        }

        nrpy_odiegm_evolve_apply(d->e, d->c, d->s, &system, &current_position, current_position+step, &step, y);
        

        exception_handler(current_position,y);
        const_eval(current_position,y,&cp);
        assign_constants(c,&cp);

	// New printing routine for TOV solver
	// Order: r, e, rho_b, P, M, nu, 2*nu, Iso_r

	fprintf(fp2, "%15.14e",current_position); //r
	fprintf(fp2, " %15.14e", c[0]); //e(or rho)
	fprintf(fp2, " %15.14e", c[1]); //rho_b
	fprintf(fp2, " %15.14e", y[0]); //P
	fprintf(fp2, " %15.14e", y[2]); //M
	fprintf(fp2, " %15.14e", y[1]); //nu
	fprintf(fp2, " %15.14e", 2*y[1]); //2*nu
	fprintf(fp2, " %15.14e", y[3]); //Iso_r
    
        fprintf(fp2,"\n");

        if (do_we_terminate(current_position, y, &cp) == 1) {
            i = size-1;
             
        } 
        if (i == size-1) { 
            printf("FINAL: Position:,\t%15.14e,\t",current_position);
            for (int n = 0; n < number_of_equations; n++) {
                printf("Equation %i:,\t%15.14e,\t",n, y[n]);
            }
            for (int n = 0; n < number_of_constants; n++) {
                printf("Constant %i:,\t%15.14e,\t",n, c[n]);
            }
            printf("\n");
        }
    }


    fclose(fp2);

    nrpy_odiegm_driver_free(d);
    free(Gamma);
    free(rhoBound);
    NRPyEOS_free_memory(&ghl_eos);
    printf("ODE Solver \"Odie\" for NRml_TOV Shutting Down...\n");
    return 0;
    
} // -GM, master of dogs (OdieGM code)
  // -David Boyer, Edited for NRml_TOV
