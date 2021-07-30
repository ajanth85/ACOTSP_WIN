/*

       AAAA    CCCC   OOOO   TTTTTT   SSSSS  PPPPP
      AA  AA  CC     OO  OO    TT    SS      PP  PP
      AAAAAA  CC     OO  OO    TT     SSSS   PPPPP
      AA  AA  CC     OO  OO    TT        SS  PP
      AA  AA   CCCC   OOOO     TT    SSSSS   PP

######################################################
##########    ACO algorithms for the TSP    ##########
######################################################

      Version: 1.0
      File:    main.c
      Author:  Thomas Stuetzle
      Purpose: main routines and control for the ACO algorithms
      Check:   README and gpl.txt
      Copyright (C) 2002  Thomas Stuetzle
*/

/***************************************************************************

    Program's name: acotsp

    Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for the 
    symmetric TSP 

    Copyright (C) 2004  Thomas Stuetzle

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    email: stuetzle no@spam informatik.tu-darmstadt.de
    mail address: Universitaet Darmstadt
                  Fachbereich Informatik
                  Hochschulstr. 10
                  D-64283 Darmstadt
		  Germany

***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "ants.h"
#include "utilities.h"
#include "InOut.h"
#include "TSP.h"
#include "parse.h"
//#include "dos_timer.h"
#include "timer.h"
#include "ls.h"
#include <cfloat>

long int termination_condition( struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper, clock_t* start_time, double* elapsed)
/*    
      FUNCTION:       checks whether termination condition is met 
      INPUT:          none
      OUTPUT:         0 if condition is not met, number neq 0 otherwise
      (SIDE)EFFECTS:  none
*/
{
  return ( ((hyper->n_tours >= hyper->max_tours) && (elapsed_time( elapsed, start_time ) >= hyper->max_time)) || 
	  (ant_hyper->best_so_far_ant->tour_length <= hyper->optimal));  // ((*ant_hyper->best_so_far_ant).tour_length
}



void construct_solutions( struct problem* instance,  int n, struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper, long int *seed)
/*    
      FUNCTION:       manage the solution construction phase
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  when finished, all ants of the colony have constructed a solution  
*/
{
    long int k;        /* counter variable */
    long int step;    /* counter of the number of construction steps */

    TRACE ( printf("construct solutions for all ants\n"); );

    /* Mark all cities as unvisited */
    for ( k = 0 ; k < ant_hyper->n_ants ; k++) {
	ant_empty_memory( &(ant_hyper->ant[k]), n );
    }
    
    step = 0; 
    /* Place the ants on same initial city */
    for (k = 0; k < ant_hyper->n_ants; k++) {
        place_ant(&(ant_hyper->ant[k]), step, n, seed);
    }

    while ( step < n-1 ) {
	step++;

	for ( k = 0 ; k < ant_hyper->n_ants ; k++ ) {
	    neighbour_choose_and_move_to_next( &(ant_hyper->ant[k]), step, instance, n, ant_hyper, seed);
	    if (ant_hyper->acs_flag )
		local_acs_pheromone_update( &(ant_hyper->ant[k]), step, instance, n, ant_hyper);
	}
    }

    step = n;
    for ( k = 0 ; k < ant_hyper->n_ants ; k++ ) {
        ant_hyper->ant[k].tour[n] = ant_hyper->ant[k].tour[0];
        ant_hyper->ant[k].tour_length = compute_tour_length(ant_hyper->ant[k].tour, instance, n );
	if (ant_hyper->acs_flag )
	    local_acs_pheromone_update( &(ant_hyper->ant[k]), step, instance, n, ant_hyper);
    }
    hyper->n_tours += ant_hyper->n_ants;
}



void init_try( long int ntry, struct problem *instance, int n, struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper, struct ls_parameters* ls_params, long int *seed, clock_t *start_time, double *elapsed)
/*    
      FUNCTION: initilialize variables appropriately when starting a trial
      INPUT:    trial number
      OUTPUT:   none
      COMMENTS: none
*/
{

    TRACE ( printf("INITIALIZE TRIAL\n"); );
  
    start_timers(start_time);
    hyper->time_used = elapsed_time(elapsed, start_time);
    hyper->time_passed = hyper->time_used;

    fprintf(hyper->comp_report,"seed %ld\n",*seed);
    fflush(hyper->comp_report);
    /* Initialize variables concerning statistics etc. */
  
    hyper->n_tours      = 1;
    hyper->iteration    = 1;
    hyper->restart_iteration = 1;
    hyper->lambda       = 0.05;            
    ant_hyper->best_so_far_ant->tour_length = INFTY;  // (*ant_hyper->best_so_far_ant).tour_length = INFTY;
    hyper->found_best   = 0;
  
    /* Initialize the Pheromone trails, only if ACS is used, pheromones
       have to be initialized differently */
    if ( !(ant_hyper->acs_flag || ant_hyper->mmas_flag || ant_hyper->bwas_flag) ) {
        ant_hyper->trail_0 = 1. / ( (ant_hyper->rho) * nn_tour(instance, n, hyper, ant_hyper, ls_params, seed) );
	/* in the original papers on Ant System, Elitist Ant System, and
	   Rank-based Ant System it is not exactly defined what the
	   initial value of the pheromones is. Here we set it to some
	   small constant, analogously as done in MAX-MIN Ant System.  
	*/
	init_pheromone_trails(ant_hyper->trail_0, n, ant_hyper);
    } 
    if (ant_hyper->bwas_flag ) {
        ant_hyper->trail_0 = 1. / ( (double) n * (double) nn_tour(instance, n, hyper, ant_hyper, ls_params, seed)) ;
	init_pheromone_trails(ant_hyper->trail_0, n, ant_hyper);
    } 
    if (ant_hyper->mmas_flag ) {
        ant_hyper->trail_max = 1. / ( (ant_hyper->rho) * nn_tour(instance, n, hyper, ant_hyper, ls_params, seed) );
        ant_hyper->trail_min = ant_hyper->trail_max / ( 2. * n );
	init_pheromone_trails(ant_hyper->trail_max, n, ant_hyper);
    }
    if (ant_hyper->acs_flag ) {
        ant_hyper->trail_0 = 1. / ( (double) n * (double) nn_tour( instance, n, hyper, ant_hyper, ls_params, seed) ) ;
	init_pheromone_trails(ant_hyper->trail_0, n, ant_hyper);
    }
  
    /* Calculate combined information pheromone times heuristic information */
    compute_total_information(instance, n, ant_hyper);
    
    fprintf(hyper->comp_report,"begin try %li \n",ntry);
    fprintf(hyper->stat_report,"begin try %li \n",ntry);
}



void local_search(struct problem *instance, int n, struct ant_hyper_parameters* ant_hyper, struct ls_parameters* ls_params, long int *seed)
/*    
      FUNCTION:       manage the local search phase; apply local search to ALL ants; in 
                      dependence of ls_flag one of 2-opt, 2.5-opt, and 3-opt local search
		      is chosen.
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants of the colony have locally optimal tours
      COMMENTS:       typically, best performance is obtained by applying local search 
                      to all ants. It is known that some improvements (e.g. convergence 
		      speed towards high quality solutions) may be obtained for some 
		      ACO algorithms by applying local search to only some of the ants.
		      Overall best performance is typcially obtained by using 3-opt.
*/
{
    long int k;

    TRACE ( printf("apply local search to all ants\n"); );

    for ( k = 0 ; k < ant_hyper->n_ants ; k++ ) {
	if ( ls_params->ls_flag == 1 )
	    two_opt_first(ant_hyper->ant[k].tour, instance, n, ls_params, seed);    /* 2-opt local search */
	else if ( ls_params->ls_flag == 2 )
	    two_h_opt_first(ant_hyper->ant[k].tour, instance, n, ls_params, seed );  /* 2.5-opt local search */
	else if ( ls_params->ls_flag == 3 )
	    three_opt_first(ant_hyper->ant[k].tour, instance, n, ls_params, seed);  /* 3-opt local search */
	else {
	    fprintf(stderr,"type of local search procedure not correctly specified\n");
	    exit(1);
	}
    ant_hyper->ant[k].tour_length = compute_tour_length(ant_hyper->ant[k].tour, instance, n );
    }
}



void update_statistics(struct problem *instance, int n, struct hyper_parameters* hyper, struct ant_hyper_parameters *ant_hyper, struct ls_parameters* ls_params, clock_t* start_time, double* elapsed)
/*    
      FUNCTION:       manage some statistical information about the trial, especially
                      if a new best solution (best-so-far or restart-best) is found and
                      adjust some parameters if a new best solution is found
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated; trail_min 
                      and trail_max used by MMAS may be updated
*/
{

    long int iteration_best_ant;
    double p_x; /* only used by MMAS */

    iteration_best_ant = find_best(ant_hyper); /* iteration_best_ant is a global variable */

    if ( ant_hyper->ant[iteration_best_ant].tour_length < ant_hyper->best_so_far_ant->tour_length ) {

        hyper->time_used = elapsed_time( elapsed, start_time ); /* best sol found after time_used */
	copy_from_to( &(ant_hyper->ant[iteration_best_ant]), ant_hyper->best_so_far_ant, n );
	copy_from_to( &(ant_hyper->ant[iteration_best_ant]), ant_hyper->restart_best_ant, n );

    hyper->found_best = hyper->iteration;
    hyper->restart_found_best = hyper->iteration;
    hyper->found_branching = node_branching(hyper->lambda, instance, n, ant_hyper);
    hyper->branching_factor = hyper->found_branching;
	if (ant_hyper->mmas_flag ) {
	    if ( !ls_params->ls_flag ) {
		p_x = exp(log(0.05)/n); 
        ant_hyper->trail_min = 1. * (1. - p_x) / (p_x * (double)((ant_hyper->nn_ants + 1) / 2));
        ant_hyper->trail_max = 1. / ( (ant_hyper->rho) * ant_hyper->best_so_far_ant->tour_length );
        ant_hyper->trail_0 = ant_hyper->trail_max;
        ant_hyper->trail_min = ant_hyper->trail_max * ant_hyper->trail_min;
	    } else {
            ant_hyper->trail_max = 1. / ( (ant_hyper->rho) * ant_hyper->best_so_far_ant->tour_length );
            ant_hyper->trail_min = ant_hyper->trail_max / ( 2. * n );
            ant_hyper->trail_0 = ant_hyper->trail_max;
	    }
	}
	write_report(hyper, ant_hyper, start_time, elapsed);
    }
    if (ant_hyper->ant[iteration_best_ant].tour_length < ant_hyper->restart_best_ant->tour_length ) {
	copy_from_to( &(ant_hyper->ant[iteration_best_ant]), ant_hyper->restart_best_ant, n );
    hyper->restart_found_best = hyper->iteration;
	printf("restart best: %ld, restart_found_best %ld, time %.2f\n",ant_hyper->restart_best_ant->tour_length, hyper->restart_found_best, elapsed_time ( elapsed, start_time ));
    }
}



void search_control_and_statistics(struct problem *instance, int n, struct hyper_parameters *hyper, struct ant_hyper_parameters *ant_hyper, clock_t* start_time, double* elapsed)
/*    
      FUNCTION:       occasionally compute some statistics and check whether algorithm 
                      is converged 
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  restart-best and best-so-far ant may be updated; trail_min 
                      and trail_max used by MMAS may be updated
*/
{
    TRACE ( printf("SEARCH CONTROL AND STATISTICS\n"); );

    if (!(hyper->iteration % 100)) {
	population_statistics(instance, n, hyper, ant_hyper);
	hyper->branching_factor = node_branching(hyper->lambda, instance, n, ant_hyper);
	printf("\nbest so far %ld, iteration: %ld, time %.2f, b_fac %.5f\n",(*ant_hyper->best_so_far_ant).tour_length,hyper->iteration,elapsed_time( elapsed, start_time),hyper->branching_factor);
    
	if ( ant_hyper->mmas_flag && (hyper->branching_factor < hyper->branch_fac) && (hyper->iteration - hyper->restart_found_best > 250) ) {
	    /* MAX-MIN Ant System was the first ACO algorithm to use
	       pheromone trail re-initialisation as implemented
	       here. Other ACO algorithms may also profit from this mechanism.
	    */
	    printf("INIT TRAILS!!!\n"); (*ant_hyper->restart_best_ant).tour_length = INFTY; 
	    init_pheromone_trails( ant_hyper->trail_max, n, ant_hyper);
	    compute_total_information(instance, n, ant_hyper);
	    hyper->restart_iteration = hyper->iteration;
	    hyper->restart_time = elapsed_time( elapsed, start_time );
	}
	printf("try %li, iteration %li, b-fac %f \n\n", hyper->n_try,hyper->iteration,hyper->branching_factor);
    }
}



void as_update( int n, struct ant_hyper_parameters *ant_hyper)
/*    
      FUNCTION:       manage global pheromone deposit for Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants deposit pheromones on matrix "pheromone"
*/
{
    long int   k;

    TRACE ( printf("Ant System pheromone deposit\n"); );

    for ( k = 0 ; k < ant_hyper->n_ants ; k++ )
	global_update_pheromone( &(ant_hyper->ant[k]), n, ant_hyper);
}



void eas_update( int n, struct ant_hyper_parameters *ant_hyper)
/*    
      FUNCTION:       manage global pheromone deposit for Elitist Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants plus elitist ant deposit pheromones on matrix "pheromone"
*/
{
    long int   k;

    TRACE ( printf("Elitist Ant System pheromone deposit\n"); );

    for ( k = 0 ; k < ant_hyper->n_ants ; k++ )
	global_update_pheromone( &(ant_hyper->ant[k]) , n, ant_hyper);
    global_update_pheromone_weighted(ant_hyper->best_so_far_ant, ant_hyper->elitist_ants, n, ant_hyper);
}



void ras_update( int n, struct ant_hyper_parameters *ant_hyper)
/*    
      FUNCTION:       manage global pheromone deposit for Rank-based Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  the ras_ranks-1 best ants plus the best-so-far ant deposit pheromone 
                      on matrix "pheromone"
      COMMENTS:       this procedure could be implemented slightly faster, but it is 
                      anyway not critical w.r.t. CPU time given that ras_ranks is 
		      typically very small.
*/
{
    long int i, k, b, target;
    long int *help_b;

    TRACE ( printf("Rank-based Ant System pheromone deposit\n"); );

    help_b = malloc(ant_hyper->n_ants  * sizeof(long int) );
    for ( k = 0 ; k < ant_hyper->n_ants ; k++ )
	help_b[k] = ant_hyper->ant[k].tour_length;

    for ( i = 0 ; i < ant_hyper->ras_ranks-1 ; i++ ) {
	b = help_b[0]; target = 0;
	for ( k = 0 ; k < ant_hyper->n_ants ; k++ ) {
	    if ( help_b[k] < b ) {
		b = help_b[k]; target = k;
	    }
	}
	help_b[target] = LONG_MAX;
	global_update_pheromone_weighted( &(ant_hyper->ant[target]), ant_hyper->ras_ranks-i-1 , n, ant_hyper);
    }
    global_update_pheromone_weighted(ant_hyper->best_so_far_ant, ant_hyper->ras_ranks, n, ant_hyper);
    free ( help_b );
}



void mmas_update( int n, struct hyper_parameters *hyper, struct ant_hyper_parameters *ant_hyper, struct ls_parameters* ls_params)
/*    
      FUNCTION:       manage global pheromone deposit for MAX-MIN Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  either the iteration-best or the best-so-far ant deposit pheromone 
                      on matrix "pheromone"
*/
{
    /* we use default upper pheromone trail limit for MMAS and hence we
       do not have to worry regarding keeping the upper limit */

    long int iteration_best_ant;

    TRACE ( printf("MAX-MIN Ant System pheromone deposit\n"); );

    if ( hyper->iteration % ant_hyper->u_gb ) {
	iteration_best_ant = find_best(ant_hyper);
	global_update_pheromone( &(ant_hyper->ant[iteration_best_ant]), n, ant_hyper);
    }
    else {
	if (ant_hyper->u_gb == 1 && (hyper->restart_found_best - hyper->iteration > 50))
	    global_update_pheromone(ant_hyper->best_so_far_ant , n, ant_hyper);
	else
	    global_update_pheromone(ant_hyper->restart_best_ant , n, ant_hyper);
    }

    if ( ls_params->ls_flag ) {
	/* implement the schedule for u_gb as defined in the 
	   Future Generation Computer Systems article or in Stuetzle's PhD thesis.
	   This schedule is only applied if local search is used.
	*/
	if ( (hyper->iteration - hyper->restart_iteration ) < 25 )
        ant_hyper->u_gb = 25;
	else if ( (hyper->iteration - hyper->restart_iteration) < 75 )
        ant_hyper->u_gb = 5;
	else if ( (hyper->iteration - hyper->restart_iteration) < 125 )
        ant_hyper->u_gb = 3;
	else if ( (hyper->iteration - hyper->restart_iteration) < 250 )
        ant_hyper->u_gb = 2;
	else 
        ant_hyper->u_gb = 1;
    } else
        ant_hyper->u_gb = 25;
  
}



void bwas_update( int n, struct hyper_parameters *hyper, struct ant_hyper_parameters *ant_hyper, long int *seed, clock_t* start_time, double* elapsed)
/*    
      FUNCTION:       manage global pheromone deposit for Best-Worst Ant System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  either the iteration-best or the best-so-far ant deposit pheromone 
                      on matrix "pheromone"
*/
{
    long int   iteration_worst_ant, distance_best_worst;

    TRACE ( printf("Best-worst Ant System pheromone deposit\n"); );

    global_update_pheromone(ant_hyper->best_so_far_ant, n, ant_hyper);
    iteration_worst_ant = find_worst(ant_hyper);
    bwas_worst_ant_update(ant_hyper->best_so_far_ant, &(ant_hyper->ant[iteration_worst_ant]), n, ant_hyper);  //Used to be &ant[iteration_worst_ant]
    distance_best_worst = distance_between_ants(ant_hyper->best_so_far_ant, &(ant_hyper->ant[iteration_worst_ant]), n);  // Used to be &ant[iteration_worst_ant]
/*    printf("distance_best_worst %ld, tour length worst %ld\n",distance_best_worst,ant[iteration_worst_ant].tour_length); */
    if ( distance_best_worst < (long int) (0.05 * (double) n) ) {
        ant_hyper->restart_best_ant->tour_length = INFTY; //Used to be (*restart_best_ant).tour_length = INFTY;
	init_pheromone_trails(ant_hyper->trail_0, n, ant_hyper);
    hyper->restart_iteration = hyper->iteration;
    hyper->restart_time = elapsed_time( elapsed, start_time );
	printf("init pheromone trails with %.15f, iteration %ld\n", ant_hyper->trail_0, hyper->iteration);
    }
    else 
	bwas_pheromone_mutation(n, hyper, ant_hyper, seed);
}



void acs_global_update(struct problem *instance, int n, struct ant_hyper_parameters *ant_hyper)
/*    
      FUNCTION:       manage global pheromone deposit for Ant Colony System
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  the best-so-far ant deposits pheromone on matrix "pheromone"
      COMMENTS:       global pheromone deposit in ACS is done per default using 
                      the best-so-far ant; Gambardella & Dorigo examined also iteration-best
		      update (see their IEEE Trans. on Evolutionary Computation article), 
		      but did not use it for the published computational results.
*/
{
    TRACE ( printf("Ant colony System global pheromone deposit\n"); );

    global_acs_pheromone_update(ant_hyper->best_so_far_ant, instance, n, ant_hyper);
}



void pheromone_trail_update(struct problem *instance, int n, struct hyper_parameters *hyper, struct ant_hyper_parameters *ant_hyper, struct ls_parameters* ls_params, long int *seed, clock_t* start_time, double* elapsed)
/*    
      FUNCTION:       manage global pheromone trail update for the ACO algorithms
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  pheromone trails are evaporated and pheromones are deposited 
                      according to the rules defined by the various ACO algorithms.
*/
{
    /* Simulate the pheromone evaporation of all pheromones; this is not necessary 
       for ACS (see also ACO Book) */
    if ( ant_hyper->as_flag || ant_hyper->eas_flag || ant_hyper->ras_flag || ant_hyper->bwas_flag || ant_hyper->mmas_flag ) {
	if ( ls_params->ls_flag ) {
	    if (ant_hyper->mmas_flag )
		mmas_evaporation_nn_list(instance, n, ant_hyper);
	    else
		evaporation_nn_list(instance, n, ant_hyper);
	    /* evaporate only pheromones on arcs of candidate list to make the 
	       pheromone evaporation faster for being able to tackle large TSP 
	       instances. For MMAS additionally check lower pheromone trail limits.
	    */
	} else {
	    /* if no local search is used, evaporate all pheromone trails */
	    evaporation(n, ant_hyper);
	}
    }

    /* Next, apply the pheromone deposit for the various ACO algorithms */
    if (ant_hyper->as_flag )
	as_update(n, ant_hyper); 
    else if (ant_hyper->eas_flag )
	eas_update(n, ant_hyper);
    else if (ant_hyper->ras_flag )
	ras_update(n, ant_hyper);
    else if (ant_hyper->mmas_flag )
	mmas_update(n, hyper, ant_hyper, ls_params);
    else if (ant_hyper->bwas_flag )
	bwas_update(n, hyper, ant_hyper, seed, start_time, elapsed);
    else if (ant_hyper->acs_flag )
	acs_global_update(instance, n, ant_hyper);

  /* check pheromone trail limits for MMAS; not necessary if local
     search is used, because in the local search case lower pheromone trail
     limits are checked in procedure mmas_evaporation_nn_list */
    if (ant_hyper->mmas_flag && !ls_params->ls_flag )
	check_pheromone_trail_limits(n, ant_hyper);

  /* Compute combined information pheromone times heuristic info after
     the pheromone update for all ACO algorithms except ACS; in the ACS case 
     this is already done in the pheromone update procedures of ACS */
    if (ant_hyper->as_flag || ant_hyper->eas_flag || ant_hyper->ras_flag || ant_hyper->mmas_flag || ant_hyper->bwas_flag ) {
	if ( ls_params->ls_flag ) {
	    compute_nn_list_total_information(instance, n, ant_hyper);
	} else {
	    compute_total_information(instance, n, ant_hyper);
	}
    }
}


void internal_main(int argc, char* argv[]) {
    /*
      FUNCTION:       main control for running the ACO algorithms
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  none
      COMMENTS:       this function controls the run of "max_tries" independent trials

*/

    const char* progname;
    struct options input_options;

    progname = argv[0] != NULL && *(argv[0]) != '\0'
        ? argv[0]
        : "acotsp";

    parse_options(&input_options, progname, argc, argv);
    omp_set_num_threads(atoi(input_options.arg_threads));

    char thesis_report_stats[100];
    
    char file_type[10];
    sprintf(file_type, "%.*s", strlen(input_options.arg_tsplibfile) - 4, input_options.arg_tsplibfile);
    sprintf(thesis_report_stats, "thesis_report_stats_%d.%s", omp_get_max_threads(), file_type);
    FILE* thesis_stats = fopen(thesis_report_stats, "w");
    fprintf(thesis_stats, "threads,avg_best_time,avg_best_iteration,avg_best_tour,std_time,std_iteration,std_tour\n");

    struct all_threads_result* all_threads_result = malloc(sizeof(*all_threads_result));
    all_threads_result->avg_best_iteration = INT_MAX;
    all_threads_result->avg_best_time = DBL_MAX;
    all_threads_result->avg_best_tour = INT_MAX;
    all_threads_result->avg_std_iteration = 0;
    all_threads_result->avg_std_time = 0;
    all_threads_result->avg_std_tour = 0;

    // Open MP directives
#pragma omp parallel
    {
    long int i;

    
    // Previosuly global variables
    clock_t* start_time  = malloc(sizeof(*start_time));
    double elapsed;
    long int seed;
    long int n = 0;          /* number of cities in the instance to be solved */
    struct thread_result* thread_result = malloc(sizeof(*thread_result));

    struct problem *instance =  malloc(sizeof(*instance));
    struct hyper_parameters *hyperparams =  malloc (sizeof(*hyperparams)) ;
    hyperparams->name_buf = calloc(LINE_BUF_LEN, sizeof(char));

    struct ant_hyper_parameters* ant_hyperparams = malloc(sizeof (*ant_hyperparams));
    struct ls_parameters* ls_params =malloc(sizeof (*ls_params));

    start_timers(start_time);
    // Open MP directives
#pragma critical 
        {
            init_program(argc, argv, instance, &n, hyperparams, ant_hyperparams, ls_params, &seed);

            seed = (long int)time(NULL);

            instance->nn_list = compute_nn_lists(instance, n, ant_hyperparams, ls_params);
            ant_hyperparams->pheromone = generate_double_matrix(n, n);
            ant_hyperparams->total = generate_double_matrix(n, n);

            hyperparams->time_used = elapsed_time(&elapsed, start_time);
            printf("Initialization took %.10f seconds\n", hyperparams->time_used);
        }

        for (hyperparams->n_try = 0; hyperparams->n_try < hyperparams->max_tries; hyperparams->n_try++) {

            init_try(hyperparams->n_try, instance, n, hyperparams, ant_hyperparams, ls_params, &seed, start_time, &elapsed);

            while (!termination_condition(hyperparams, ant_hyperparams, start_time, &elapsed)) {

                construct_solutions(instance, n, hyperparams, ant_hyperparams, &seed);

                if (ls_params->ls_flag > 0)
                    local_search(instance, n, ant_hyperparams, ls_params, &seed);

                update_statistics(instance, n, hyperparams, ant_hyperparams, ls_params, start_time, &elapsed);

                pheromone_trail_update(instance, n, hyperparams, ant_hyperparams, ls_params, &seed, start_time, &elapsed);

                search_control_and_statistics(instance, n, hyperparams, ant_hyperparams, start_time, &elapsed);

                hyperparams->iteration++;
            }
            exit_try(hyperparams->n_try, instance, n, hyperparams, ant_hyperparams, start_time, &elapsed, thread_result);
        }
        exit_program(instance, hyperparams, thread_result);


#pragma critical
        {
            all_threads_result->avg_std_iteration += thread_result->std_iteration;
            all_threads_result->avg_std_time += thread_result->std_time;
            all_threads_result->avg_std_tour += thread_result->std_tour;

            if (thread_result->avg_best_tour < all_threads_result->avg_best_tour)
            {
                all_threads_result->avg_best_iteration = thread_result->avg_best_iteration;
                all_threads_result->avg_best_time = thread_result->avg_best_time;
                all_threads_result->avg_best_tour = thread_result->avg_best_tour;
                all_threads_result->number_of_threads = omp_get_max_threads();
            }
            else if (thread_result->avg_best_tour == all_threads_result->avg_best_tour)
            {
                if (thread_result->avg_best_time < all_threads_result->avg_best_time)
                {
                    all_threads_result->avg_best_iteration = thread_result->avg_best_iteration;
                    all_threads_result->avg_best_time = thread_result->avg_best_time;
                    all_threads_result->avg_best_tour = thread_result->avg_best_tour;
                    all_threads_result->number_of_threads = omp_get_max_threads();
                }
            }
        }

        free(instance->distance);
        free(instance->nn_list);
        free(ant_hyperparams->pheromone);
        free(ant_hyperparams->total);
        free(hyperparams->best_in_try);
        free(hyperparams->best_found_at);
        free(hyperparams->time_best_found);
        free(hyperparams->time_total_run);
        for (i = 0; i < ant_hyperparams->n_ants; i++) {
            free(ant_hyperparams->ant[i].tour);
            free(ant_hyperparams->ant[i].visited);
        }
        free(ant_hyperparams->ant);
        free(ant_hyperparams->best_so_far_ant->tour);
        free(ant_hyperparams->best_so_far_ant->visited);
        free(ant_hyperparams->prob_of_selection);

        free(instance);
        free(hyperparams);
        free(ant_hyperparams);
        free(ls_params);
        free(thread_result);
    }

    all_threads_result->avg_std_iteration /= all_threads_result->number_of_threads;
    all_threads_result->avg_std_time /= all_threads_result->number_of_threads;
    all_threads_result->avg_std_tour /= all_threads_result->number_of_threads;

    printf("Best tour: %d\tBest iteration: %d\tBest time: %f\tNumber of threads: %d\n", all_threads_result->avg_best_tour, all_threads_result->avg_best_iteration, all_threads_result->avg_best_time, all_threads_result->number_of_threads);
    fprintf(thesis_stats, "%d,%f,%d,%d,%f,%f,%f\n", 
        all_threads_result->number_of_threads, 
        all_threads_result->avg_best_time,
        all_threads_result->avg_best_iteration, 
        all_threads_result->avg_best_tour, 
        all_threads_result->avg_std_time, 
        all_threads_result->avg_std_iteration,
        all_threads_result->avg_std_tour);

    free(all_threads_result);
    fclose(thesis_stats);
}

/* --- main program ------------------------------------------------------ */

int main(int argc, char *argv[]) {
/*    
      FUNCTION:       main control for running the ACO algorithms
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  none
      COMMENTS:       this function controls the run of "max_tries" independent trials
*/
    internal_main(argc, argv);

    return(0);
}
