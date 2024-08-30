#include "NRml_TOV_funcs.c" //nrpy_odiegm itself, for use in NRml_TOV.
#include "NRml_TOV_user_methods.c" //Edited specifically for the TOV solver.
#include "../GRHayL/GRHayL/include/nrpyeos_tabulated.h" //Check that this filepath is correct.
#include <ghl.h> //Need the GRHayL Library

/********************************
// UPDATE: NRml_TOV is a special version of the Odie main driver. It is built to handle the TOV equations specifically, and is built to handle a variety of different types of EOS
// Simple Polytrope:
//	-- TODO: Explain SP
// Piecewise Polytrope:
//	-- TODO: Explain PP
// Tabulated EOS:
//	-- TODO: Explain Tab
//
// TODO: Add Comments on how to use the program
//
// TODO: Once this code is complete, you should probably organize some of these comments. I would like it to look a little cleaner.
//
// David Boyer (TODO:Date Completed?) 
********************************/
int main(){
  
  //TABULATED EOS TEMPLATE
  
  
    printf("Beginning ODE Solver \"Odie\" for NRml_TOV...\n");

    //Choose your EOS type and (if necessary) path
    ghl_eos_parameters ghl_eos={};
    char type = 't'; //s=simple, p=piecewise, t=tabulated
    char* table_path = "/home/boye4060/projects/New_Odie_TOV/NRml_TOV/EOS_Tables/SLy4_3335_rho391_temp163_ye66_adjusted.h5"; //Replace with your own path to your table.
    
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
    int adams_bashforth_order = 4; //Only edit if you are using an AB method
    bool adaptive_step = true;
    double error_limit = 1e-16;
    
    char file_name[] = "ooData.txt"; 

    const nrpy_odiegm_step_type * step_type;
    step_type = nrpy_odiegm_step_ARKF;
    const nrpy_odiegm_step_type * step_type_2;
    step_type_2 = step_type;

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
    if(type=='t'){
      printf("Tabulated EOS\n");
      ghl_initialize_tabulated_eos_functions_and_params(table_path,rhoAtm,rhoMin,rhoMax,YeAtm,YeMin,YeMax,TAtm,TMin,TMax,&ghl_eos);
      ghl_tabulated_compute_Ye_P_eps_of_rho_beq_constant_T(T_in, &ghl_eos);
    }
    else{
      printf("ERROR: Invalid EOS type. Must be 't' (tabulated) for the template\n");
      printf("Shutting down due to error...\n");
      return 1;}

    
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
    //However, double negatives are confusing for user interfacing.
    //We use the postive variable to declare the GSL compatible variable here.

    bool no_adaptive_step;
    if (adaptive_step == true){
      no_adaptive_step=false;}
    else{
      no_adaptive_step=true;}
    
    double absolute_error_limit = error_limit;
    double relative_error_limit = error_limit;
    int number_of_equations = 4;
    int number_of_constants = 2;
    struct constant_parameters cp = { 
      .dimension = number_of_constants
    };
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
    
    double y[number_of_equations];
    

    double c[number_of_constants];
    
    get_initial_condition(y,&cp); 
    //const_eval(current_position, y,&cp);
    assign_constants(c,&cp); 

    FILE *fp2;
    fp2 = fopen(file_name,"w");
    printf("Printing to file '%s'.\n",file_name);

    
    // Print the order: r, e(rho), rho_b, P, M, Lapse, Lapse_2, Iso_r
    printf("INITIAL: Position:,\t%f,\t",current_position);

    fprintf(fp2, "%15.14e",current_position); //r
    fprintf(fp2, " %15.14e", c[0]); //e(or rho)
    fprintf(fp2, " %15.14e", c[1]); //rho_b
    fprintf(fp2, " %15.14e", y[0]); //P
    fprintf(fp2, " %15.14e", y[2]); //M
    fprintf(fp2, " %15.14e", exp(y[1])); //TODO: Currently Lapse. Calculate Conformal Factor
    fprintf(fp2, " %15.14e", exp(2*y[1])); //TODO: Same as above for Conformal Factor 2.
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
        

        // Printing section.
        // Uncomment for live updates. Prints to the file automatically.
        //printf("%15.14e\n",current_position);

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
    NRPyEOS_free_memory(&ghl_eos);
    //if(type == 't'){
    //NRPyEOS_tabulated_free_beq_quantities(&ghl_eos);
      //ghl_tabulated_free_beq_quantities(&ghl_eos);
    //}
    printf("ODE Solver \"Odie\" for NRml_TOV Shutting Down...\n");
    return 0;
    
} // -GM, master of dogs (OdieGM code)
  // -David Boyer, Edited for NRml_TOV
