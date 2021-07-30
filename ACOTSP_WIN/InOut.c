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
      File:    InOut.c
      Author:  Thomas Stuetzle
      Purpose: mainly input / output / statistic routines
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <time.h>

#include "InOut.h"
#include "TSP.h"
#include "dos_timer.h"
#include "utilities.h"
#include "ants.h"
#include "ls.h"
#include "parse.h"



struct point * read_etsp(const char *tsp_file_name, struct problem *instance, int *n) 
/*    
      FUNCTION: parse and read instance file
      INPUT:    instance name
      OUTPUT:   list of coordinates for all nodes
      COMMENTS: Instance files have to be in TSPLIB format, otherwise procedure fails
*/
{
    FILE         *tsp_file;
    char         buf[LINE_BUF_LEN];
    long int     i, j;
    struct point *nodeptr;

    tsp_file = fopen(tsp_file_name, "r");
    if ( tsp_file == NULL ) {
	fprintf(stderr,"No instance file specified, abort\n");
	exit(1);
    }
    assert(tsp_file != NULL);
    printf("\nreading tsp-file %s ... \n\n", tsp_file_name);

    fscanf(tsp_file,"%s", buf);
    while ( strcmp("NODE_COORD_SECTION", buf) != 0 ) {
	if ( strcmp("NAME", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s ", buf); )
	    fscanf(tsp_file, "%s", buf);
	    strcpy(instance->name, buf);
	    TRACE ( printf("%s \n", instance->name); )
	    buf[0]=0;
	}
	else if ( strcmp("NAME:", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    strcpy(instance->name, buf);
	    TRACE ( printf("%s \n", instance->name); )
	    buf[0]=0;
	}
	else if ( strcmp("COMMENT", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    TRACE ( printf("%s", buf); )
	    buf[0]=0;
	}
	else if ( strcmp("COMMENT:", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    TRACE ( printf("%s", buf); )
	    buf[0]=0;
	}
	else if ( strcmp("TYPE", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s ", buf); )
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s\n", buf); )
	    if( strcmp("TSP", buf) != 0 ) {
		fprintf(stderr,"\n Not a TSP instance in TSPLIB format !!\n");
		exit(1);
	    }
	    buf[0]=0;
	}
	else if ( strcmp("TYPE:", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s\n", buf); )
	    if( strcmp("TSP", buf) != 0 ) {
		fprintf(stderr,"\n Not a TSP instance in TSPLIB format !!\n");
		exit(1);
	    }
	    buf[0]=0;
	}
	else if( strcmp("DIMENSION", buf) == 0 ){
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s ", buf); );
	    fscanf(tsp_file, "%ld", n);
	    instance->n = *n;
	    TRACE ( printf("%ld\n", n); );
	    assert ( *n > 2 && *n < 6000);
	    buf[0]=0;
	}
	else if ( strcmp("DIMENSION:", buf) == 0 ) {
	    fscanf(tsp_file, "%ld", n);
	    instance->n = *n;
	    TRACE ( printf("%ld\n", *n); );
	    assert ( *n > 2 && *n < 6000);
	    buf[0]=0;
	}
	else if( strcmp("DISPLAY_DATA_TYPE", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    TRACE ( printf("%s", buf); );
	    buf[0]=0;
	}
	else if ( strcmp("DISPLAY_DATA_TYPE:", buf) == 0 ) {
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    TRACE ( printf("%s", buf); );
	    buf[0]=0;
	}
	else if( strcmp("EDGE_WEIGHT_TYPE", buf) == 0 ){
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s ", buf); );
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s\n", buf); );
	    if ( strcmp("EUC_2D", buf) == 0 ) {
		distance = round_distance;
	    }
	    else if ( strcmp("CEIL_2D", buf) == 0 ) {
		distance = ceil_distance;
	    }
	    else if ( strcmp("GEO", buf) == 0 ) {
		distance = geo_distance;
	    }
	    else if ( strcmp("ATT", buf) == 0 ) {
		distance = att_distance;
	    }
	    else if ( strcmp("EUC_TOROID", buf) == 0 ) {
	        distance = toroid_distance;
	   }
	    else
		fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
	    strcpy(instance->edge_weight_type, buf);
	    buf[0]=0;
	}
	else if( strcmp("EDGE_WEIGHT_TYPE:", buf) == 0 ){
	    /* set pointer to appropriate distance function; has to be one of 
	       EUC_2D, CEIL_2D, GEO, ATT or EUC_TOROID. Everything else fails */
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s\n", buf); )
		printf("%s\n", buf);
	    printf("%s\n", buf);
	    if ( strcmp("EUC_2D", buf) == 0 ) {
		distance = round_distance;
	    }
	    else if ( strcmp("CEIL_2D", buf) == 0 ) {
		distance = ceil_distance;
	    }
	    else if ( strcmp("GEO", buf) == 0 ) {
		distance = geo_distance;
	    }
	    else if ( strcmp("ATT", buf) == 0 ) {
		distance = att_distance;
	    }
	    else if ( strcmp("EUC_TOROID", buf) == 0 ) {
	        distance = toroid_distance;
	   }
	    else {
		fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
		exit(1);
	    }
	    strcpy(instance->edge_weight_type, buf);
	    buf[0]=0;
	}
	buf[0]=0;
	fscanf(tsp_file,"%s", buf);
    }


    if( strcmp("NODE_COORD_SECTION", buf) == 0 ){
	TRACE ( printf("found section contaning the node coordinates\n"); )
	    }
    else{
	fprintf(stderr,"\n\nSome error ocurred finding start of coordinates from tsp file !!\n");
	exit(1);
    }

    if( (nodeptr = malloc(sizeof(struct point) * *n)) == NULL )
	exit(EXIT_FAILURE);
    else {
	for ( i = 0 ; i < *n ; i++ ) {
	    fscanf(tsp_file,"%ld %lf %lf", &j, &nodeptr[i].x, &nodeptr[i].y );
	}
    }
    TRACE ( printf("number of cities is %ld\n",*n); )
    TRACE ( printf("\n... done\n"); )
	return (nodeptr);
}



void write_report( struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper, clock_t* start_time, double* elapsed)
/*    
      FUNCTION: output some info about trial (best-so-far solution quality, time)
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
*/
{
    printf("best %ld, iteration: %ld, time %.2f\n",(*ant_hyper->best_so_far_ant).tour_length, hyper->iteration,elapsed_time(elapsed, start_time));
    fprintf(hyper->comp_report,"best %ld\t iteration %ld\t tours %ld\t time %.3f\n",(*ant_hyper->best_so_far_ant).tour_length, hyper->iteration,hyper->n_tours,hyper->time_used);
}



void print_default_parameters(struct hyper_parameters* hyper, struct ant_hyper_parameters* ant_hyper, struct ls_parameters* ls_params)
/*    
      FUNCTION: output default parameter settings
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
*/
{
    fprintf(stderr,"\nDefault parameter settings are:\n\n");
    fprintf(stderr,"max_tries\t\t %ld\n", hyper->max_tries);
    fprintf(stderr,"max_tours\t\t %ld\n", hyper->max_tours);
    fprintf(stderr,"max_time\t\t %.2f\n", hyper->max_time);
    fprintf(stderr,"optimum\t\t\t %ld\n", hyper->optimal);
    fprintf(stderr,"n_ants\t\t\t %ld\n", ant_hyper->n_ants);
    fprintf(stderr,"nn_ants\t\t\t %ld\n", ant_hyper->nn_ants);
    fprintf(stderr,"alpha\t\t\t %.2f\n", ant_hyper->alpha);
    fprintf(stderr,"beta\t\t\t %.2f\n", ant_hyper->beta);
    fprintf(stderr,"rho\t\t\t %.2f\n", ant_hyper->rho);
    fprintf(stderr,"q_0\t\t\t %.2f\n", ant_hyper->q_0);
    fprintf(stderr,"elitist_ants\t\t 0\n");
    fprintf(stderr,"ras_ranks\t\t 6\n");
    fprintf(stderr,"ls_flag\t\t\t %ld\n", ls_params->ls_flag);
    fprintf(stderr,"nn_ls\t\t\t %ld\n", ls_params->nn_ls);
    fprintf(stderr,"dlb_flag\t\t %ld\n", ls_params->dlb_flag);
    fprintf(stderr,"as_flag\t\t\t %ld\n", ant_hyper->as_flag);
    fprintf(stderr,"eas_flag\t\t %ld\n", ant_hyper->eas_flag);
    fprintf(stderr,"ras_flag\t\t %ld\n", ant_hyper->ras_flag);
    fprintf(stderr,"mmas_flag\t\t %ld\n", ant_hyper->mmas_flag);
    fprintf(stderr,"bwas_flag\t\t %ld\n", ant_hyper->bwas_flag);
    fprintf(stderr,"acs_flag\t\t %ld\n", ant_hyper->acs_flag);
}



void set_default_parameters(struct hyper_parameters* hyper, struct ant_hyper_parameters* ant_hyper, struct ls_parameters* ls_params)
/*    
      FUNCTION: set default parameter settings
      INPUT:    none
      OUTPUT:   none
      COMMENTS: none
*/
{
    ls_params->ls_flag        = 3;     /* per default run 3-opt*/
    ls_params->dlb_flag       = TRUE;  /* apply don't look bits in local search */
    ls_params->nn_ls          = 20;    /* use fixed radius search in the 20 nearest neighbours */
    ant_hyper->n_ants         = 25;    /* number of ants */
    ant_hyper->nn_ants        = 20;    /* number of nearest neighbours in tour construction */
    ant_hyper->alpha          = 1.0;
    ant_hyper->beta           = 2.0;
    ant_hyper->rho            = 0.5;
    ant_hyper->q_0            = 0.0;
    hyper->max_tries      = 10;
    hyper->max_tours      = 100;
    hyper->max_time       = 10.0;
    hyper->optimal        = 1;
    hyper->branch_fac     = 1.00001;
    ant_hyper->as_flag        = FALSE;
    ant_hyper->eas_flag       = FALSE;
    ant_hyper->ras_flag       = FALSE;
    ant_hyper->mmas_flag      = TRUE;
    ant_hyper->u_gb           = INFTY;
    ant_hyper->bwas_flag      = FALSE;
    ant_hyper->acs_flag       = FALSE;
    ant_hyper->ras_ranks      = 6;
    ant_hyper->elitist_ants   = 100;
}



void population_statistics (struct problem *instance, int n, struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper)
/*    
      FUNCTION:       compute some population statistics like average tour length, 
                      standard deviations, average distance, branching-factor and 
		      output to a file gathering statistics
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  none
*/
{
    long int j, k;
    long int *l;
    double   pop_mean, pop_stddev, avg_distance = 0.0;
    
    l = calloc(ant_hyper->n_ants, sizeof(long int));
    for( k = 0 ; k < ant_hyper->n_ants ; k++ ) {
	l[k] = ant_hyper->ant[k].tour_length;
    }
    
    pop_mean = mean( l, ant_hyper->n_ants);
    pop_stddev = std_deviation( l, ant_hyper->n_ants, pop_mean );
    hyper->branching_factor = node_branching(hyper->lambda, instance, n, ant_hyper);
    
    for ( k = 0 ; k < ant_hyper->n_ants-1 ; k++ )
	for ( j = k+1 ; j < ant_hyper->n_ants ; j++) {
	    avg_distance += (double)distance_between_ants( &(ant_hyper->ant[k]), &(ant_hyper->ant[j]), n);
	}
    avg_distance /= ((double)ant_hyper->n_ants * (double)(ant_hyper->n_ants-1) / 2.);
    

    fprintf(hyper->stat_report,"%ld\t%.1f\t%.5f\t%.7f\t%.5f\t%.1f\t%.1f\t%.5f\n", hyper->iteration, pop_mean, pop_stddev, pop_stddev / pop_mean, hyper->branching_factor,  (hyper->branching_factor - 1.) * (double)n, avg_distance, avg_distance / (double)n);
}



double node_branching(double l, struct problem *instance, int n, struct ant_hyper_parameters* ant_hyper)
/*    
      FUNCTION:       compute the average node lambda-branching factor 
      INPUT:          lambda value 
      OUTPUT:         average node branching factor 
      (SIDE)EFFECTS:  none
      COMMENTS:       see the ACO book for a definition of the average node 
                      lambda-branching factor 
*/
{
  long int  i, m;
  double    min, max, cutoff;
  double    avg;
  double    *num_branches;

  num_branches = calloc(n, sizeof(double));

  for ( m = 0 ; m < n ; m++ ) {
    /* determine max, min to calculate the cutoff value */
    min = ant_hyper->pheromone[m][instance->nn_list[m][1]];
    max = ant_hyper->pheromone[m][instance->nn_list[m][1]];
    for ( i = 1 ; i < ant_hyper->nn_ants ; i++ ) {
      if (ant_hyper->pheromone[m][instance->nn_list[m][i]] > max )
	max = ant_hyper->pheromone[m][instance->nn_list[m][i]];
      if (ant_hyper->pheromone[m][instance->nn_list[m][i]] < min )
	min = ant_hyper->pheromone[m][instance->nn_list[m][i]];
    }
    cutoff = min + l * (max - min);
    
    for ( i = 0 ; i < ant_hyper->nn_ants ; i++ ) {
      if (ant_hyper->pheromone[m][instance->nn_list[m][i]] > cutoff )
	num_branches[m] += 1.;
    }
  }
  avg = 0.;
  for ( m = 0 ; m < n ; m++ ) {
    avg += num_branches[m];
  }
  free ( num_branches );
  /* Norm branching factor to minimal value 1 */
  return ( avg / (double) (n * 2)  );
}



void output_solution( struct problem *instance, int n, struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper)
/*    
      FUNCTION:       output a solution together with node coordinates
      INPUT:          none
      OUTPUT:         none
      COMMENTS:       not used in the default implementation but may be useful anyway
*/
{

    long int i;
    
    for ( i = 0 ; i < n ; i++ ) {
	fprintf(hyper->stat_report," %ld %f %f\n",(*ant_hyper->best_so_far_ant).tour[i],instance->nodeptr[(*ant_hyper->best_so_far_ant).tour[i]].x, instance->nodeptr[(*ant_hyper->best_so_far_ant).tour[i]].y);
    }
    printf("\n"); 
}



void exit_try( long int ntry, struct problem *instance, int n, struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper, clock_t* start_time, double* elapsed)
/*    
      FUNCTION:       save some statistical information on a trial once it finishes
      INPUT:          trial number
      OUTPUT:         none
      COMMENTS:       
*/
{
  checkTour( (*ant_hyper->best_so_far_ant).tour, instance, n);
/*    printTourFile( (*best_so_far_ant).tour ); */

  printf("\n Best Solution in try %ld is %ld\n",ntry, (*ant_hyper->best_so_far_ant).tour_length);
  fprintf(hyper->report,"Best: %ld\t Iterations: %6ld\t B-Fac %.5f\t Time %.2f\t Tot.time %.2f\n", (*ant_hyper->best_so_far_ant).tour_length, hyper->found_best, hyper->found_branching, hyper->time_used, elapsed_time( elapsed, start_time ));
  printf(" Best Solution was found after %ld iterations\n", hyper->found_best);

  hyper->best_in_try[ntry] = (*ant_hyper->best_so_far_ant).tour_length;
  hyper->best_found_at[ntry] = hyper->found_best;
  hyper->time_best_found[ntry] = hyper->time_used;
  hyper->time_total_run[ntry] = elapsed_time(elapsed, start_time);
  printf("\ntry %ld, Best %ld, found at iteration %ld, found at time %f\n",ntry, hyper->best_in_try[ntry], hyper->best_found_at[ntry], hyper->time_best_found[ntry]);

  fprintf(hyper->comp_report,"end try %ld\n\n",ntry);
  fprintf(hyper->stat_report,"end try %ld\n\n",ntry);
  output_solution(instance, n, hyper, ant_hyper);
  fflush(hyper->report);
  fflush(hyper->comp_report);
  fflush(hyper->stat_report);

}



void exit_program( struct problem *instance, struct hyper_parameters *hyper, struct thread_result* thread_result)
/*    
      FUNCTION:       save some final statistical information on a trial once it finishes
      INPUT:          none
      OUTPUT:         none
      COMMENTS:       
*/
{
  long int best_tour_length, worst_tour_length;
  double   t_avgbest, t_stdbest, t_avgtotal, t_stdtotal;
  double   avg_sol_quality = 0., avg_cyc_to_bst = 0., stddev_best, stddev_iterations;

  best_tour_length = best_of_vector(hyper->best_in_try , hyper->max_tries );
  worst_tour_length = worst_of_vector(hyper->best_in_try , hyper->max_tries );

  avg_cyc_to_bst = mean(hyper->best_found_at , hyper->max_tries );
  stddev_iterations = std_deviation(hyper->best_found_at, hyper->max_tries, avg_cyc_to_bst );

  avg_sol_quality = mean(hyper->best_in_try , hyper->max_tries );
  stddev_best = std_deviation(hyper->best_in_try, hyper->max_tries, avg_sol_quality);

  t_avgbest = meanr(hyper->time_best_found, hyper->max_tries );
  printf(" t_avgbest = %f\n", t_avgbest );
  t_stdbest = std_deviationr(hyper->time_best_found, hyper->max_tries, t_avgbest);

  t_avgtotal = meanr(hyper->time_total_run, hyper->max_tries );
  printf(" t_avgtotal = %f\n", t_avgtotal );
  t_stdtotal = std_deviationr(hyper->time_total_run, hyper->max_tries, t_avgtotal);

  fprintf(hyper->report,"\nAverage-Best: %.2f\t Average-Iterations: %.2f", avg_sol_quality, avg_cyc_to_bst);
  fprintf(hyper->report,"\nStddev-Best: %.2f \t Stddev Iterations: %.2f", stddev_best, stddev_iterations);
  fprintf(hyper->report,"\nBest try: %ld\t\t Worst try: %ld\n", best_tour_length, worst_tour_length);
  fprintf(hyper->report,"\nAvg.time-best: %.2f stddev.time-best: %.2f\n", t_avgbest, t_stdbest);
  fprintf(hyper->report,"\nAvg.time-total: %.2f stddev.time-total: %.2f\n", t_avgtotal, t_stdtotal);

  thread_result->avg_best_iteration = avg_cyc_to_bst;
  thread_result->avg_best_time = t_avgbest;
  thread_result->avg_best_tour = avg_sol_quality;
  thread_result->std_iteration = stddev_iterations;
  thread_result->std_tour = stddev_best;
  thread_result->std_time = t_stdbest;

  if (hyper->optimal > 0 ) {
    fprintf(hyper->report," excess best = %f, excess average = %f, excess worst = %f\n",(double)(best_tour_length - hyper->optimal) / (double)hyper->optimal,(double)(avg_sol_quality - hyper->optimal) / (double)hyper->optimal,(double)(worst_tour_length - hyper->optimal) / (double)hyper->optimal);
  }

  fprintf(hyper->comp_report,"end problem %s\n",instance->name);
}



void init_program(long int argc, char* argv[], struct problem *instance, int *n, struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper, struct ls_parameters* ls_params, long int *seed)
/*
      FUNCTION:       initialize the program,
      INPUT:          program arguments, needed for parsing commandline
      OUTPUT:         none
      COMMENTS:
*/
{
    char temp_buffer[LINE_BUF_LEN];

    printf(PROG_ID_STR);
    set_default_parameters(hyper, ant_hyper, ls_params);
    setbuf(stdout, NULL);
    parse_commandline(argc, argv, hyper, ant_hyper, ls_params);

    assert(ant_hyper->n_ants < MAX_ANTS - 1);
    assert(ant_hyper->nn_ants < MAX_NEIGHBOURS);
    assert(ant_hyper->nn_ants > 0);
    assert(ls_params->nn_ls > 0);
    assert(hyper->max_tries <= MAXIMUM_NO_TRIES);

    hyper->best_in_try = calloc(hyper->max_tries, sizeof(long int));
    hyper->best_found_at = calloc(hyper->max_tries, sizeof(long int));
    hyper->time_best_found = calloc(hyper->max_tries, sizeof(double));
    hyper->time_total_run = calloc(hyper->max_tries, sizeof(double));

    *seed = (long int)time(NULL);


    TRACE(printf("read problem data  ..\n\n"); )
        instance->nodeptr = read_etsp(hyper->name_buf, instance, n);
    TRACE(printf("\n .. done\n\n"); )

        ls_params->nn_ls = MIN(*n - 1, ls_params->nn_ls);

    sprintf(temp_buffer, "thead_%d_best.%s", omp_get_thread_num(), instance->name);
    TRACE(printf("%s\n", temp_buffer); )
        hyper->report = fopen(temp_buffer, "w");
    sprintf(temp_buffer, "cmp.%s", instance->name);
    TRACE(printf("%s\n", temp_buffer); )
        hyper->comp_report = fopen(temp_buffer, "w");
    sprintf(temp_buffer, "stat.%s", instance->name);
    TRACE(printf("%s\n", temp_buffer); )
        hyper->stat_report = fopen(temp_buffer, "w");

    printf("calculating distance matrix ..\n\n");
    instance->distance = compute_distances(instance, *n);
    printf(" .. done\n");

    write_params(hyper, ant_hyper, ls_params);
    fprintf(hyper->comp_report, "begin problem %s\n", hyper->name_buf);
    printf("allocate ants' memory ..\n\n");
    allocate_ants(*n, ant_hyper);
    printf(" .. done\n");
    assert(ls_params->nn_ls < *n && 0 < ls_params->nn_ls);

        printf("\nFinally set ACO algorithm specific parameters, typically done as proposed in literature\n\n");
    /* default setting for elitist_ants is 0; if EAS is applied and
       option elitist_ants is not used, we set the default to
       elitist_ants = n */
    if (ant_hyper->eas_flag && ant_hyper->elitist_ants == 0)
        ant_hyper->elitist_ants = *n;
}



void printDist(struct problem *instance, int n) 
/*    
      FUNCTION:       print distance matrix 
      INPUT:          none
      OUTPUT:         none
*/
{
  long int i,j;

  printf("Distance Matrix:\n");
  for ( i = 0 ; i < n ; i++) {
    printf("From %ld:  ",i);
    for ( j = 0 ; j < n - 1 ; j++ ) {
      printf(" %ld", instance->distance[i][j]);
    }
    printf(" %ld\n", instance->distance[i][n-1]);
    printf("\n");
  }
  printf("\n");
}



void printHeur(struct problem *instance, int n) 
/*    
      FUNCTION:       print heuristic information 
      INPUT:          none
      OUTPUT:         none
*/
{
  long int i, j;

  printf("Heuristic information:\n");
  for ( i = 0 ; i < n ; i++) {
    printf("From %ld:  ",i);
    for ( j = 0 ; j < n - 1 ; j++ ) {
      printf(" %.3f ", HEURISTIC(i,j));
    }
    printf(" %.3f\n", HEURISTIC(i,j));
    printf("\n");
  }
  printf("\n");
}



void printTrail(int n, struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper)
/*    
      FUNCTION:       print pheromone trail values 
      INPUT:          none
      OUTPUT:         none
*/
{
  long int i,j;

  printf("pheromone Trail matrix, iteration: %ld\n\n",hyper->iteration);
  for ( i = 0 ; i < n ; i++) {
    printf("From %ld:  ",i);
    for ( j = 0 ; j < n ; j++ ) {
      printf(" %.10f ", ant_hyper->pheromone[i][j]);
      if (ant_hyper->pheromone[i][j] > 1.0 )
	printf("XXXXX\n");
    }
    printf("\n");
  }
  printf("\n");
}



void printTotal(int n, struct ant_hyper_parameters* ant_hyper)
/*    
      FUNCTION:       print values of pheromone times heuristic information 
      INPUT:          none
      OUTPUT:         none
*/
{
  long int i, j;

  printf("combined pheromone and heuristic info\n\n");
  for (i=0; i < n; i++) {
    for (j = 0; j < n - 1 ; j++) {
      printf(" %.15f &", ant_hyper->total[i][j]);
      if (ant_hyper->total[i][j] > 1.0 )
	printf("XXXXX\n");
    }
    printf(" %.15f\n", ant_hyper->total[i][n-1]);
    if (ant_hyper->total[i][n-1] > 1.0 )
      printf("XXXXX\n");
  }
  printf("\n");
}



void printProbabilities(int n, struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper)
/*    
      FUNCTION:       prints the selection probabilities as encountered by an ant 
      INPUT:          none
      OUTPUT:         none
      COMMENTS:       this computation assumes that no choice has been made yet. 
*/
{
  long int i, j;
  double   *p;
  double   sum_prob;

  printf("Selection Probabilities, iteration: %ld\n",hyper->iteration);
  p = calloc( n, sizeof(double) );

  for (i=0; i < n; i++) {
    printf("From %ld:  ",i);
    sum_prob = 0.;
    for ( j = 0 ; j < n ; j++) {
      if ( i == j )
	p[j] = 0.;
      else
	p[j] = ant_hyper->total[i][j];
      sum_prob += p[j];
    }
    for ( j = 0 ; j < n ; j++) {
      p[j] = p[j] / sum_prob;
    }
    for ( j = 0 ; j < n-1 ; j++) {
      printf(" %.5f ", p[j]);
    }
    printf(" %.5f\n", p[n-1]);
    if (!(j % 26)) {
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");
  free ( p );
}



void printTour( long int *t, struct problem *instance, int n)
/*    
      FUNCTION:       print the tour *t
      INPUT:          pointer to a tour
      OUTPUT:         none
*/
{
    long int   i;

    printf("\n");
    for( i = 0 ; i <= n ; i++ ) {
	if (!i%25)
	    printf("\n");
	printf("%ld ", t[i]);
    }
    printf("\n");
    printf("Tour length = %ld\n\n",compute_tour_length( t, instance, n ));
}



void printTourFile( long int *t, struct problem *instance, int n, struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper)
/*    
      FUNCTION:       print the tour *t to cmp.tsplibfile
      INPUT:          pointer to a tour
      OUTPUT:         none
*/
{
    long int   i;

    fprintf(hyper->comp_report,"begin solution\n");
    for( i = 0 ; i < n ; i++ ) {
	fprintf(hyper->comp_report,"%ld ", t[i]);
    }
    fprintf(hyper->comp_report,"\n");
    fprintf(hyper->comp_report,"Tour length %ld\n",compute_tour_length( t, instance, n));
    fprintf(hyper->comp_report,"end solution\n");
}



void checkTour( long int *t, struct problem *instance, int n)
/*    
      FUNCTION:       make a simple check whether tour *t can be feasible
      INPUT:          pointer to a tour
      OUTPUT:         none
*/
{
    long int   i, sum=0;

    for( i = 0 ; i < n ; i++ ) {
	sum += t[i];
    }
    if ( sum != (n-1) * n / 2 ) {
	fprintf(stderr,"Next tour must be flawed !!\n");
	printTour( t, instance, n );
	exit(1);
    }
}



void write_params( struct hyper_parameters *hyper, struct ant_hyper_parameters *ant_hyper, struct ls_parameters* ls_params)
/*    
      FUNCTION:       writes chosen parameter settings in standard output and in 
                      report files 
      INPUT:          none
      OUTPUT:         none
*/
{
  printf("\nParameter-settings: \n\n");
  printf("max-tries %ld\n", hyper->max_tries);
  printf("max-tours %ld\n", hyper->max_tours);
  printf("optimum %ld\n", hyper->optimal);
  printf("time %f\n", hyper->max_time);
  printf("num-ants %ld\n", ant_hyper->n_ants);
  printf("num-neigh %ld\n", ant_hyper->nn_ants);
  printf("alpha %f\n", ant_hyper->alpha);
  printf("beta %f\n", ant_hyper->beta);
  printf("rho %f\n", ant_hyper->rho);
  printf("q_0 %f\n", ant_hyper->q_0);
  printf("branch-up %f\n", hyper->branch_fac);
  printf("ls_flag %ld\n", ls_params->ls_flag);
  printf("nn_ls %ld\n", ls_params->nn_ls);
  printf("dlb_flag %ld\n", ls_params->dlb_flag);
  printf("as_flag %ld\n", ant_hyper->as_flag);
  printf("eas_flag %ld\n", ant_hyper->eas_flag);
  printf("ras_flag %ld\n", ant_hyper->ras_flag);
  printf("mmas_flag %ld\n", ant_hyper->mmas_flag);
  printf("bwas_flag %ld\n", ant_hyper->bwas_flag);
  printf("acs_flag %ld\n", ant_hyper->acs_flag);
  printf("\n");
  fprintf(hyper->report,"\nParameter-settings: \n\n");
  fprintf(hyper->report,"max-tries %ld\n", hyper->max_tries);
  fprintf(hyper->report,"max-tours %ld\n", hyper->max_tours);
  fprintf(hyper->report,"optimum %ld\n", hyper->optimal);
  fprintf(hyper->report,"time %f\n", hyper->max_time);
  fprintf(hyper->report,"num-ants %ld\n", ant_hyper->n_ants);
  fprintf(hyper->report,"num-neigh %ld\n", ant_hyper->nn_ants);
  fprintf(hyper->report,"alpha %f\n", ant_hyper->alpha);
  fprintf(hyper->report,"beta %f\n", ant_hyper->beta);
  fprintf(hyper->report,"rho %f\n", ant_hyper->rho);
  fprintf(hyper->report,"q_0 %f\n", ant_hyper->q_0);
  fprintf(hyper->report,"branch-up %f\n", hyper->branch_fac);
  fprintf(hyper->report,"ls_flag %ld\n", ls_params->ls_flag);
  fprintf(hyper->report,"nn_ls %ld\n", ls_params->nn_ls);
  fprintf(hyper->report,"dlb_flag %ld\n", ls_params->dlb_flag);
  fprintf(hyper->report,"as_flag %ld\n", ant_hyper->as_flag);
  fprintf(hyper->report,"eas_flag %ld\n", ant_hyper->eas_flag);
  fprintf(hyper->report,"ras_flag %ld\n", ant_hyper->ras_flag);
  fprintf(hyper->report,"mmas_flag %ld\n", ant_hyper->mmas_flag);
  fprintf(hyper->report,"bwas_flag %ld\n", ant_hyper->bwas_flag);
  fprintf(hyper->report,"acs_flag %ld\n", ant_hyper->acs_flag);
  fprintf(hyper->report,"\n");
  fprintf(hyper->comp_report,"%s",PROG_ID_STR);
  fprintf(hyper->comp_report,"\nParameter-settings: \n\n");
  fprintf(hyper->comp_report,"max-tries %ld\n", hyper->max_tries);
  fprintf(hyper->comp_report,"max-tours %ld\n", hyper->max_tours);
  fprintf(hyper->comp_report,"optimum %ld\n", hyper->optimal);
  fprintf(hyper->comp_report,"time %f\n", hyper->max_time);
  fprintf(hyper->comp_report,"num-ants %ld\n", ant_hyper->n_ants);
  fprintf(hyper->comp_report,"num-neigh %ld\n", ant_hyper->nn_ants);
  fprintf(hyper->comp_report,"alpha %f\n", ant_hyper->alpha);
  fprintf(hyper->comp_report,"beta %f\n", ant_hyper->beta);
  fprintf(hyper->comp_report,"rho %f\n", ant_hyper->rho);
  fprintf(hyper->comp_report,"q_0 %f\n", ant_hyper->q_0);
  fprintf(hyper->comp_report,"branch-up %f\n", hyper->branch_fac);
  fprintf(hyper->comp_report,"ls_flag %ld\n", ls_params->ls_flag);
  fprintf(hyper->comp_report,"nn_ls %ld\n", ls_params->nn_ls);
  fprintf(hyper->comp_report,"dlb_flag %ld\n", ls_params->dlb_flag);
  fprintf(hyper->comp_report,"as_flag %ld\n", ant_hyper->as_flag);
  fprintf(hyper->comp_report,"eas_flag %ld\n", ant_hyper->eas_flag);
  fprintf(hyper->comp_report,"ras_flag %ld\n", ant_hyper->ras_flag);
  fprintf(hyper->comp_report,"mmas_flag %ld\n", ant_hyper->mmas_flag);
  fprintf(hyper->comp_report,"bwas_flag %ld\n", ant_hyper->bwas_flag);
  fprintf(hyper->comp_report,"acs_flag %ld\n", ant_hyper->acs_flag);
  fprintf(hyper->comp_report,"\n");
}











