
#include "nrpy_odiegm_proto.c"

//This file actually holds the ODE evolution code. Ripped directly from OdieGM with no changes made.
//For use in NRml_TOV.
//DO NOT TOUCH THIS FILE!!! All edits that need to be made are in the NRml_TOV_Driver.c file.
//For more info, see OdieGM and the GSL ODE solver.






// #include "nrpy_odiegm_proto.c"

// This file contains the actual definitions for the funcitons outlined in nrpy_odiegm_proto.c

// Memory allocation functions.
nrpy_odiegm_step *
nrpy_odiegm_step_alloc (const nrpy_odiegm_step_type * T, size_t dim)
{
  // Allocate the step "object", set all values, even those that may not be used. 
  nrpy_odiegm_step *s = (nrpy_odiegm_step *) malloc (sizeof (nrpy_odiegm_step));
  s->type = T;
  s->method_type = 1;
  s->adams_bashforth_order = 0;
  s->rows = T->rows;
  s->columns = T->columns;
  // these last two assignments might be unecessary, but it will be convenient if this number
  // can be acessed at both levels. 
  if (T->rows == T->columns) {
    s->method_type = 0; // aka, normal RK-type method. 
  }
  if (T->rows == 19) {
    s->method_type = 2; // AB method. 
    s->adams_bashforth_order = 4; // default order chosen, if user wants control they will 
    // specify elsewhere after allocation is run.  
  }

  s->y_values = (double *) malloc ((double)19.0 * dim * sizeof (double));
  // This here is the array used to store past values.
  // Only used for AB methods, but it still needs to be dynamically allocated. 
  // Having an adams_bashforth_order of 0 doesn't throw any errors, which is conveinent.

  return s;
}

nrpy_odiegm_evolve *
nrpy_odiegm_evolve_alloc (size_t dim)
{
  // Allocate the evolve "object" and set all values, even those that may not be used.
  nrpy_odiegm_evolve *e = (nrpy_odiegm_evolve *) malloc (sizeof (nrpy_odiegm_evolve));
  e->y0 = (double *) malloc (dim * sizeof (double));
  e->yerr = (double *) malloc (dim * sizeof (double));
  // Fill these with 0 just in case someone tries to allocate something. 
  for (int n = 0; n < dim; n++) {
    e->y0[n] = 0.0;
    e->yerr[n] = 0.0;
  }
  
  e->count = 0;
  e->last_step = 0.0; // By default we don't use this value. 
  e->bound = 0.0; // This will be adjusted when the first step is taken.
  e->current_position = 0.0; //This will be regularly adjusted as the program goes on. 
  e->no_adaptive_step = false; // We assume adaptive by default. 
  return e;
}

nrpy_odiegm_control *
nrpy_odiegm_control_y_new (double eps_abs, double eps_rel)
{
    // Allocate the control "object." Unusual wording of function name is due to us needing
    // a GSL replacement. 
    nrpy_odiegm_control *c = (nrpy_odiegm_control *) malloc (sizeof (nrpy_odiegm_control));
    c->abs_lim = eps_abs;
    c->rel_lim = eps_rel;
    //defaults
    
    c->scale_factor = 0.9;
    c->error_safety = 4.0/15.0;
    c->ay_error_scaler = 1.0;
    c->ady_error_scaler = 1.0;
    c->max_step_adjustment = 5.0;
    c->min_step_adjustment = 0.2;
    c->absolute_max_step = 0.1;
    c->absolute_min_step = 1e-10;
    c->error_upper_tolerance = 1.1;
    c->error_lower_tolerance = 0.5;
    
    // These are all the default values, virtually all responsible for adaptive timestep and 
    // error estimation.

    return c;
}

nrpy_odiegm_driver * nrpy_odiegm_driver_alloc_y_new (const nrpy_odiegm_system * sys,
                               const nrpy_odiegm_step_type * T,
                               const double hstart,
                               const double epsabs, const double epsrel)
{
    // Initializes an ODE driver "object" which contains all the "objets" above, making a system
    // that is prepared to evaluate a system of differential equations. 

    nrpy_odiegm_driver *state;
    state = (nrpy_odiegm_driver *) calloc (1, sizeof (nrpy_odiegm_driver));
    const size_t dim = sys->dimension; 
    state->sys = sys;
    state->s = nrpy_odiegm_step_alloc (T, dim);

    state->e = nrpy_odiegm_evolve_alloc (dim);
    state->h = hstart; // the step size. 

    state->c = nrpy_odiegm_control_y_new (epsabs, epsrel);

  // There were functions here in GSL that assigned the driver to the objects contained in the driver.
  // We will not be doing that insanity. 

  return state;
}

// Memory freeing functions. 
void nrpy_odiegm_control_free (nrpy_odiegm_control * c)
{
  free (c);
}
void nrpy_odiegm_evolve_free (nrpy_odiegm_evolve * e)
{
  free (e->yerr);
  free (e->y0);
  free (e);
}
void nrpy_odiegm_step_free (nrpy_odiegm_step * s)
{ 
  free (s->y_values);
  free (s);
}
void nrpy_odiegm_driver_free (nrpy_odiegm_driver * state)
{
  // In most cases, this method should be called alone, calling the others would be redundant. 
  if (state->c)
    nrpy_odiegm_control_free (state->c);

  if (state->e)
    nrpy_odiegm_evolve_free (state->e);

  if (state->s)
    nrpy_odiegm_step_free (state->s);

  free (state);
}

// The actual stepping functions follow. 

// The goal is for these functions to be completely agnostic to whatever the user is doing, 
// they should always work regardless of the form of the system passed, the method passed, and even
// if the user does something dumb it shouldn't crash. It will spit out nonsense in those cases, though. 

int nrpy_odiegm_evolve_apply (nrpy_odiegm_evolve * e, nrpy_odiegm_control * c,
                             nrpy_odiegm_step * s,
                             const nrpy_odiegm_system * dydt, double *t,
                             double t1, double *h, double y[]) {
    // This is the big one, the function that ACTUALLY performs the step.

    // First off, check if we're at the desired edge or not. 
    if (*t + *h > t1) {
        *h = t1 - *t;
        // If we're going past an endpoint we want, reduce the step size. 
        // Otherwise continue as normal. 
        // No need to stop the adaptive time step! If we need to increase the size, we
        // Still report the smaller value, so it'll go through. 
        e->last_step = 1.0; // This is generally not used but the user might want it or something
        // to tell that this has been triggered. 
    }

    // Gotta read in several things... improves readability.
    // Don't need a million arrows everywhere if we do this. 
    int number_of_equations = (int)(dydt->dimension);
    double current_position = *t;
    e->current_position = *t;
    double step = *h; 

    unsigned long int i = e->count;
    if (i == 0) {
        e->bound = current_position;
        // If this is our first ever step, record what the starting position was. 
    }

    bool no_adaptive_step = e->no_adaptive_step;

    int method_type = s->method_type; 
    int rows = s->type->rows;
    int columns = s->type->columns;
    int adams_bashforth_order = s->adams_bashforth_order;

    double absolute_error_limit = c->abs_lim;
    double relative_error_limit = c->rel_lim;
    double scale_factor = c->scale_factor;
    double error_safety = c->error_safety;
    double ay_error_scaler = c->ay_error_scaler;
    double ady_error_scaler = c->ady_error_scaler;
    double max_step_adjustment = c-> max_step_adjustment;
    double min_step_adjustment = c->min_step_adjustment;
    double absolute_max_step = c->absolute_max_step;
    double absolute_min_step = c->absolute_min_step;
    double error_upper_tolerance = c->error_upper_tolerance;
    double error_lower_tolerance = c->error_lower_tolerance;

    double y_values[number_of_equations][adams_bashforth_order];

    int counter = 0; // This counter is reused time and time again for sifting through memory
    // Allow me to express my dislike of void pointers. 

    // The following section only runs if we're using an AB method, otherwise it jumps over. 
    if (adams_bashforth_order != 0) {
        if (i == 0) {
            // First time initialization of the y_values array for AB methods. 
            for (int n = 0; n< number_of_equations; n++) {
                y_values[n][0] = y[n];
                for (int m = 1; m < adams_bashforth_order; m++) {
                    y_values[n][m] = 0; // These values shouldn't be used, but zero them anyway. 
                } 
            }
        } else {
            // Load values from known y_values if not first step for AB method. 
            for (int n = 0; n< number_of_equations; n++) {
                for (int m = 0; m < adams_bashforth_order; m++) {
                    y_values[n][m] = *((double *)(*s).y_values+counter); // Gotta fill in an array... joy...
                    counter++;
                    // This has to be done this way due to the array being passed as a void pointer. 
                } 
            }
        }
    }

    // Read in the step type. 
    const nrpy_odiegm_step_type * step_type;
    step_type = s->type;

    counter = 0;
    if (method_type == 2) {
        rows = adams_bashforth_order;
        columns = adams_bashforth_order;
    }
    double butcher[rows][columns];
    // This is the butcher table that actually defines the method we use. 
    if (method_type != 2) { // If we aren't using AB method, just fill it without anything special. 
        for (int k=0; k < rows; k++) {
            for (int j = 0; j < columns; j++) {
                butcher[k][j] = *((double *)(*step_type).butcher+counter);
                counter++;
            }
        }
    } else { // If we ARE using an AB method, we need to construct it a little more carefully. 
        counter = counter + 19*(19-adams_bashforth_order);
        // Every row has 19 elements, and we need to clear 19-order rows, 
        // leaving only the order behind. 
        for (int i=0; i < adams_bashforth_order; i++) {
            counter = counter + 19-adams_bashforth_order; 
            // for every row, clear the unneeded zeroes. 
            for (int j = 0; j < adams_bashforth_order; j++) {
                butcher[i][j] = *((double *)(*step_type).butcher+counter);
                // This slowly counts through the array via complciated void pointer nonsense. 
                counter++;
            }
        }
    }

    if (method_type != 2) {
        // To use adaptive time-step, we need to store data at different step values:
        double y_big_step[number_of_equations];
        double y_smol_steps[number_of_equations];

        // One could argue that since the small steps will become our result 
        // we shouldn't declare it, however we are actually
        // NOT going to assign the results to the actual answer y until we compare and run the adaptive
        // time-step algorithm. We might throw out all the data and need to run it again! 
        double error_estimate[number_of_equations];
        // even if we aren't limiting the constants, we can still report their error. 
        
        double original_step = step;
        // We need to be able to refer to the original step so we can 
        // see if we're adjusting it too much at once. 
        double previous_step = step;
        // if we end up in a situation where the adaptive method wants to oscillate back and forth, 
        // we will occasionally need to know what the step we found before the current step is. 

        // We rather explicitly do not actually take any steps until we confirm the error is below what we want.
        bool error_satisfactory = false;
        bool under_error = false;
        bool over_error = false;
        // It's important to declare these outside the error_satisfactory loop 
        // since to update the stepper we need to know exactly what kind of step change we just did. 

        // This is a slapped together solution for indexing. 
        // Uses multiplication by 1 or 0 instead of an if statement on a bool. 
        int quick_patch = 1;
        if (method_type == 2) {
            quick_patch = 0;
        }
        // This constant removes certain components from consideraiton. 

        bool floored = false;
        // This is for a check hard-coded in for if we hit the *absolute minimum* step size. 
        // We have to make sure to run the loop one more time, so rather than exiting the loop
        // we set this to true and run once more. 

        while (error_satisfactory == false) {
            
            // All of the bellow values start off thinking they are the values from the 
            // previous step or initial conditions. 
            // We must reset them every time we return here.  
            for (int n = 0; n < number_of_equations; n++) {
                y_big_step[n] = y[n];
                y_smol_steps[n] = y[n];
            } 
            for (int iteration = 1; iteration < 4; iteration++) {
                // So, we want to use Adaptive Timestep methodology. 
                // This will involve evaluating each step three times, 
                // In order to compare the evolution of two different 
                // step sizes and get an error estimate. 
                // Iteration 1 performs a normal step. 
                // Iteration 2 perofrms a half step.
                // Iteration 3 performs another half step after the previous one. 
                // Naturally the half-step results are reported as truth, 
                // but we get an error estimate from the difference
                // between the two values. 

                // For inherently adaptive methods we only go through iteration 1 and 2
                // Though instead of doing a half step, we use a second evaluation built
                // into the method. 
                
                // For AB method we only go through once, but do so with some additional operations. 

                if (i == 0 && iteration == 1 && method_type == 0 && adams_bashforth_order == 0) {
                    // Don't take unecessary steps, if we are on the first step 
                    // and have no need for the large step, ignore it.
                    // Since we always want the first step to go through 
                    // don't bother calculating things we don't need. 
                    iteration = 2;
                    // This doesn't actually apply to inherently adaptive methods 
                    // since we cheat and do it in one iteration. 
                }

                double scale = 1.0;
                // This is the number we use to scale. It's either 1 or 1/2, 
                // Depending on what size step we want. 
                int shift = 0;
                // This is the number we set if we want to shift where we are evaluating from. 
                if (iteration == 1.0) {
                    // Scale remains 1
                    // Shift remains 0
                } else if (iteration == 2.0) {
                    scale = 0.5; // Using half-steps.
                    // Shfit remains 0
                } else {
                    scale = 0.5; //Using half-steps.
                    shift = 1; 
                }
                // Every time it's needed, we multiply the step by the scale. 

                double K[rows-method_type*quick_patch][number_of_equations];
                // These are the K-values that are required to evaluate RK-like methods. 
                // They will be determined based on the provided butcher table.
                // This is a 2D matrix since each diffyQ has its own set of K-values. 
                // Note that we subtract the method type from the row: 
                // adaptive RK butcher tables are larger. 

                // Since we'll be calling K while it's empty, 
                // even though there should be no errors due
                // to the way it's set up, let's go ahead and fill it with zeroes.
                for (int j = 0; j<rows-method_type*quick_patch; j++) {
                    for (int n = 0; n<number_of_equations; n++) {
                        K[j][n]=0.0;
                    }
                } 

                double y_insert[number_of_equations];
                //  We also need an array for the inserted y-values for each equation. 

                double dy_out[number_of_equations];
                //  GSL demands that we use two separate arrays for y and y', so here's y'. 

                for (int j = 1; j < rows-method_type*quick_patch; j++) {
                    // Due to the way the butcher table is formatted, 
                    // start our index at 1 and stop at the end. 
                    double x_Insert = current_position+shift*step*scale + butcher[j-1][0]*step*scale;

                    // x_Insert does not change much for different tables, 
                    // just adjust the "step correction" term.
                    // x_Insert is the same for every equation, too.

                    for (int n = 0; n < number_of_equations; n++) {
                        y_insert[n] = y_smol_steps[n];
                    } 
                    // Note that we are setting our buffer value, y_insert, to y_smol_steps. 
                    // This is because y_smol_steps is y at first, but we will need to evolve it
                    // forward two steps, so on the second small step this will be different. 
                    // (If using a method that requires that step, otherwise this is just a formality.)

                    for (int n = 1; n < columns; n++) {
                        // Once again, start at index of 1 rather than 0.
                        for (int q = 0; q < number_of_equations; q++) {
                            y_insert[q] = y_insert[q] + butcher[j-1][n]*K[n][q];
                        }
                        // Each individual y_insert portion is dependent on one of the K values.
                        // K values are initially set to zero even though technically whenever 
                        // we would use an undeclared K-value the butcher table would have a zero.
                        // You know, just in case something goes wrong. 
                    }

                    // Now we actually evaluate the differential equations.
                    dydt->function(x_Insert, y_insert, dy_out, dydt->params);
                    // y_insert goes in, dy_out comes out.

                    for (int n = 0; n < number_of_equations; n++) {
                        K[j][n] = step*scale*dy_out[n];
                        // Fill in the K-values we just calculated. 
                    } 
                }

                // Now that we have all the K-values set, we need to find 
                // the actual result in one final loop.
                for (int n = 0; n< number_of_equations; n++) {
                    K[0][n] = y_smol_steps[n]; // The 0th spot in the K-values is reserved for 
                    // holding the final value while it's being calculated. 
                    for (int j = 1; j < columns; j++) {
                        K[0][n] = K[0][n] + butcher[rows-1-method_type*quick_patch][j]*K[j][n]; 
                        // This is where the actual approximation is finally performed. 
                    }
                    y_smol_steps[n] = K[0][n]; // Set ySmol to the new estimated value. 
                }
                // Note that we specifically set ySmol to the value, not anything else. 
                // This is because we wish to avoid abusing if statements.

                if (iteration == 1) {
                    for (int n = 0; n<number_of_equations; n++) {
                        y_big_step[n] = y_smol_steps[n];
                        y_smol_steps[n] = y[n];
                        // we still need to reset the value for SmolSteps on the first iteration
                        // no matter the type of method. 
                    }
                }
                // This only runs on the first iteration, 
                // setting the big step to the right value
                // and resetting the small steps for when we actually use it. 
                // This odd structure exists purely for efficiency. 
                
                // If we are in an adaptive method situation, 
                // use that method and exit the iterations loop.
                if (method_type == 1) {
                    for (int n = 0; n< number_of_equations; n++) {
                        K[0][n] = y_smol_steps[n]; // The 0th spot in the K-values is reserved 
                        // for holding the final value while it's being calculated. 
                        for (int j = 1; j < columns; j++) {
                            K[0][n] = K[0][n] + butcher[rows-1][j]*K[j][n]; 
                            // This is where the actual approximation is finally performed. 
                            // This time we use the bottom row, not the second to bottom row 
                            // (for adaptive methods)
                        }
                        y_smol_steps[n] = K[0][n]; // Set ySmol to the new estimated value. 
                    }

                        iteration = 4; // Break out after we get to the end, 
                        // we don't need to go any further. 
                }

                if (adams_bashforth_order != 0) {
                    iteration = 4;
                    // We only iterate once for AB. Thus, break out if we are AB. 
                    for (int n = 0; n < number_of_equations; n++) {
                        y_smol_steps[n] = y_big_step[n];
                    }
                }
            }
            // Now that the step and double step have been taken 
            // (or whatever alternative method is being used),
            // time to calculate some errors and see if we move on to the next step. 
            // First, from our parameters declared at the beginning, determine what our error limit is. 
            // Using GSL's version we frist estimate our error based on what we know.
            if (i != 0 && adams_bashforth_order == 0) {
                // Literally none of this is used for the AB method. 
                for (int n = 0; n<number_of_equations; n++) {
                    error_estimate[n] = sqrt((y_big_step[n] - y_smol_steps[n])*(y_big_step[n] - y_smol_steps[n]))* error_safety;
                    // The 4/15 for error_safety is taken from GSL's solver, a 'saftey factor' 
                    // with unknown reasoning. 
                }

                double error_limiter[number_of_equations];
                // Since the definition of the error limiter uses a derivative, 
                // we cannot use it to limit the constant's error. 
                // We originally had the error limiter set its own values. 
                // GSL's formatting requries us to change this. 

                    dydt->function(current_position+step,y_smol_steps, error_limiter, dydt->params);

                // Now SmolSteps is used to set the error_limiter. 
                for (int n = 0; n<number_of_equations; n++) {
                    error_limiter[n] = absolute_error_limit + relative_error_limit*(ay_error_scaler*sqrt(y_smol_steps[n]*y_smol_steps[n]) + ady_error_scaler*step*sqrt(error_limiter[n]*error_limiter[n]));
                }
                // The error limiter is set for every equation. Now we need to perform checks.

                double ratio_ED = 0.0;
                for (int n = 0; n<number_of_equations; n++) { 
                    if (ratio_ED < error_estimate[n]/error_limiter[n]) {
                        ratio_ED = error_estimate[n]/error_limiter[n];
                        // pick out the largest of these ratios for use, every time. 
                    }
                }

                counter = 0;
                for (int n = 0; n< number_of_equations; n++) {
                    *((double *)(*e).yerr+counter) = error_estimate[n]; // Gotta fill in an array... joy...
                    counter++;
                }

                under_error = false;
                over_error = false;
                // Make sure to set our values to false every loop. 
                // These will be set to true when the exit condition is tripped. 

                if (ratio_ED >  error_upper_tolerance) {
                    // If we are 10% (or whatever value is specified) over what the error we want is, adjust. 
                    over_error = true;
                } else if (ratio_ED <= error_lower_tolerance) {
                    // If we are 50% (or whatever value is specified) under what the error we want is, adjust. 
                    under_error = true;
                }
                if (no_adaptive_step == false && step != (min_step_adjustment * original_step)) {
                    // Before adjusting, record what the step size was a second ago. 
                    previous_step = step;
                    
                    // If we have no trouble...
                    if (under_error == false && over_error == false) {
                        error_satisfactory = true;
                    }
                    // ...Say that we're cleared to move to the next step. 
                    // However, if one of them was triggered, we need to adjust. 
                    // In these cases we change the actual step size. 
                    // It is theoretically possible for both to be triggered on different equations. 
                    // In that case, over_error takes prescedent. 
                    // We would rather have more accuracy than less in odd situations like that. 

                    // These if statements perform step adjustment if needed. Based on GSL's algorithm. 
                    else if (over_error == true) {
                        step = step * scale_factor * pow(ratio_ED,-1.0/butcher[rows-1-method_type*quick_patch][0]);
                    } else { // If under_error is true and over_error is false 
                        //is the only way to get here. The true-true situation is skipped.
                        step = step * scale_factor * pow(ratio_ED,-1.0/(butcher[rows-1-method_type*quick_patch][0]+1));
                        error_satisfactory = true;
                    }

                    // Check to see if we're adjusting the step too much at once. 
                    // If we are, declare that we're done. 
                    if (step > max_step_adjustment * original_step) {
                        step = max_step_adjustment * original_step;
                        error_satisfactory = true;
                    } else if (step < min_step_adjustment * original_step){
                        step = min_step_adjustment * original_step;
                        // We still have to go through again to make sure this applies, though. 
                        // Thus there is no errorSatisfacotry = true here. 
                    }

                    if (floored == true) {
                        error_satisfactory = true;
                    } 

                    // We also declare some minium and maximum step conditions. 
                    if (step > absolute_max_step) {
                        step = absolute_max_step;
                        error_satisfactory = true;
                    } else if (step < absolute_min_step){
                        step = absolute_min_step;
                        floored = true;
                        // This is set here since we need to run through one more time, 
                        // not end right here. 
                    }

                } else {
                    error_satisfactory = true;
                    under_error = false;
                    // This area is triggered when we purposefully take single steps.
                    // Or, alternatively, when we hit the minimum step size 
                    // adjustment on the *previous* step
                    // but still needed to go through one more time. 
                }
                // With that, the step size has been changed. If error_satisfactory is still false, 
                // it goes back and performs everything again with the new step size. 
            } else {
                error_satisfactory = true;
                // We always want the *first* step to go through without change, 
                // often the first step is chosen for a specific reason. 
                // In our work this generally came from a need to plot data sets against each other. 
                // Also do this if we are using the AB method, as it has no error checks. 
            }
        }
        
        // Finally, we actually update the real answer. 
        for (int n = 0; n<number_of_equations; n++) {
            if (method_type == 1) {
                y[n]=y_big_step[n];
            } else {
                y[n]=y_smol_steps[n];
            }
            // This check is required due to the way the butcher tables are stored. 
            // There may be a more efficient way to do this. 
        }

        if (under_error == true) {
            current_position = current_position + previous_step;
            // If we had an under_error and increased the step size, 
            // Well, we kept the older points so we use that to update our current location.
            // previous_step rather than original_step since sometimes multiple iterations go through. 
        } else {
            // In any other case we use the new step.
            // Even the case where the step wasn't actually changed. 
            if (no_adaptive_step == true) {
                current_position = e->bound + (i+1)*step;
            } else {
                current_position = current_position + step;
            }
        }

        // Before, the values were Printed here. This method no longer prints, 
        // printing is done outside any method. 

        if (adams_bashforth_order > 0) {
            // At the END of every loop, we "shift" the values in the array "down" one space, 
            // that is, into the "past."
            // Present values are 0, previous step is 1, step before that is 2, etc. 
            for (int n = 0; n < number_of_equations; n++) {
                for (int m = adams_bashforth_order - 1; m > 0; m--) {
                    y_values[n][m] = y_values[n][m-1];
                    // Note that we start at the last column, m, and move the adjacent column to it. 
                    // This pushes off the value at the largest m value, 
                    // since it's far enough in the past we no longer care.
                }
                y_values[n][0] = y[n]; 
                // Present values update to what we just calculated. 
                // We have now completed stepping. 
            }  
        }
    } else {
        // This loop is for the Adams-Bashforth method, which is implemented 
        // entirely differnetly from all RK methods.
        // As such it needs an entirely different algorithm. 

        // This is normally where we would calulate the K values, 
        // but they are entirely unecessary here.

        double y_insert[number_of_equations];
        // We also need an array for the inserted y-values for each equation. 

        double dy_out[number_of_equations];
        // GSL demands that we use two separate arrays for y and y', so here's y'. 

        double x_Insert; // This is generally going to be rather simple. 

        // First, determine which row to use in the AB butcher table. 
        int current_row;
        if (i < adams_bashforth_order-1) {
            current_row = adams_bashforth_order-1-i;
            // Basically, keep track of how many steps we actually have on offer to use. 
        } else {
            current_row = 0;
            // The highest order part of the method is used when we hit a certain step. 
        }

        for (int m = adams_bashforth_order-current_row-1; m >= 0; m--) {
            // We actually need m=0 in this case, the "present" is evaluated. 
            x_Insert = e->bound + step*(i-m);
            // The "current locaiton" depends on how far in the past we are.
            for (int j = 0; j < number_of_equations ; j++) {
                y_insert[j] = y_values[j][m];
            }
            // Grab the correct y_values for the proper time/location. 

            // Now we actually evaluate the differential equations.
            dydt->function(x_Insert, y_insert, dy_out, dydt->params);

            // With that evaluation, we can change the value of y for each equation. 
            for (int n = 0; n< number_of_equations; n++) {
                y[n] = y[n] + step*butcher[current_row][m+current_row]*dy_out[n];

            }
            // Keep in mind this is procedural, y isn't right until all 
            // values of m have been cycled through. 
        }

        // At the END of every loop, we "shift" the values in the array 
        // down one space, that is, into the "past"
        // Present values are 0, previous step is 1, step before that is 2, etc. 
        for (int n = 0; n < number_of_equations; n++) {
            for (int m = adams_bashforth_order-1; m > 0; m--) {
                y_values[n][m] = y_values[n][m-1];
                // Note that we start at the last column, m, and move the adjacent column to it. 
                // This pushes off the value at the largest m value, 
                // since it's far enough in the past we no longer care.
            }
            y_values[n][0] = y[n]; 
            // Present values update to what we just calculated. 
            // We have now completed stepping. 
        }         

        current_position = e->bound+step*(i+1);
            
    }
    
    // Now we adjust any values that changed so everything outside the function can know it. 
    *h = step;
    *t = current_position;
    e->current_position = current_position;
    e->count = i+1;

    // Update y_values, very important. We spent all that time shifting everything, 
    // we need to be able to access it next time this function is called! 
    counter = 0;

    if (adams_bashforth_order != 0) {
        // Put the new y_values back into the stored array. 
        for (int n = 0; n< number_of_equations; n++) {
            for (int m = 0; m < adams_bashforth_order; m++) {
                *((double *)(*s).y_values+counter) = y_values[n][m]; // Gotta fill in an array... joy...
                counter++;
            } 
        }
    }

    // In case the user needs it for some reason we also save the result to the evolve object.
    counter = 0;
    for (int n = 0; n< number_of_equations; n++) {
        *((double *)(*e).y0+counter) = y[n]; // Gotta fill in an array... joy...
        counter++;
    }

    return 0;                      
}

int nrpy_odiegm_evolve_apply_fixed_step (nrpy_odiegm_evolve * e,
                                        nrpy_odiegm_control * con,
                                        nrpy_odiegm_step * step,
                                        const nrpy_odiegm_system * dydt,
                                        double *t, double h0,
                                        double y[]){
    // This method performs a single fixed time step. 
    e->no_adaptive_step = true;
    nrpy_odiegm_evolve_apply(e, con, step, dydt, t, *t+h0, &h0, y);

    return 0;
}

int nrpy_odiegm_driver_apply (nrpy_odiegm_driver * d, double *t,
                             const double t1, double y[]){
    // Takes as many steps as requested at the driver level. 
    // Only really useful if you don't want to report anything until the end. Which. Sure.
    while (*t < t1) {
        nrpy_odiegm_evolve_apply(d->e, d->c, d->s, d->sys, t, t1, &(d->h), y);
    }

    return 0;
}
int nrpy_odiegm_driver_apply_fixed_step (nrpy_odiegm_driver * d, double *t,
                                        const double h,
                                        const unsigned long int n,
                                        double y[]){
    // This just forces a fixed-step extrapolation. 
    d->e->no_adaptive_step = true;
    nrpy_odiegm_driver_apply(d, t, h*(double)n, y);

    return 0;
}

