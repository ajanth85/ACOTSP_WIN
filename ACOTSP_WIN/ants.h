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
      File:    ants.h
      Author:  Thomas Stuetzle
      Purpose: implementation of procedures for ants' behaviour
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



#define HEURISTIC(m,n)     (1.0 / ((double) instance->distance[m][n] + 0.1))
/* add a small constant to avoid division by zero if a distance is 
zero */

#define EPSILON            0.00000000000000000000000000000001

#define MAX_ANTS       1024    /* max no. of ants */
#define MAX_NEIGHBOURS 512     /* max. no. of nearest neighbours in candidate set */

/* Note that *tour needs to be allocated for length n+1 since the first city of 
a tour (at position 0) is repeated at position n. This is done to make the 
computation of the tour length easier 
*/
typedef struct {
  long int  *tour;
  char      *visited;
  long int  tour_length;
} ant_struct;

struct ant_hyper_parameters {
     ant_struct* ant;      /* this (array of) struct will hold the colony */
     ant_struct* best_so_far_ant;   /* struct that contains the best-so-far ant */
     ant_struct* restart_best_ant;  /* struct that contains the restart-best ant */

     double** pheromone; /* pheromone matrix, one entry for each arc */
     double** total;     /* combination of pheromone times heuristic information */

     double* prob_of_selection;


     long int n_ants;      /* number of ants */
     long int nn_ants;     /* length of nearest neighbor lists for the ants'
                                    solution construction */

     double rho;           /* parameter for evaporation */
     double alpha;         /* importance of trail */
     double beta;          /* importance of heuristic evaluate */
     double q_0;           /* probability of best choice in tour construction */


     long int as_flag;     /* = 1, run ant system */
     long int eas_flag;    /* = 1, run elitist ant system */
     long int ras_flag;    /* = 1, run rank-based version of ant system */
     long int mmas_flag;   /* = 1, run MAX-MIN ant system */
     long int bwas_flag;   /* = 1, run best-worst ant system */
     long int acs_flag;    /* = 1, run ant colony system */

     long int elitist_ants;    /* additional parameter for elitist ant system,
                        it defines the number of elitist ants */

     long int ras_ranks;       /* additional parameter for rank-based version of ant
                        system */

     double   trail_max;       /* maximum pheromone trail in MMAS */
     double   trail_min;       /* minimum pheromone trail in MMAS */
     long int u_gb;            /* every u_gb iterations update with best-so-far ant;
                        parameter used by MMAS for scheduling best-so-far update
                     */

     double   trail_0;         /* initial pheromone trail level in ACS  and BWAS */
} ;

/* Pheromone manipulation etc. */

void init_pheromone_trails ( double initial_trail, int n, struct ant_hyper_parameters* ant_hyper);

void evaporation ( int n, struct ant_hyper_parameters* ant_hyper);

void evaporation_nn_list (struct problem *instance, int n, struct ant_hyper_parameters* ant_hyper);

void global_update_pheromone ( ant_struct *a, int n, struct ant_hyper_parameters* ant_hyper);

void global_update_pheromone_weighted ( ant_struct *a, long int weight, int n, struct ant_hyper_parameters* ant_hyper);

void compute_total_information(struct problem *instance, int n, struct ant_hyper_parameters* ant_hyper);

void compute_nn_list_total_information(struct problem *instance, int n, struct ant_hyper_parameters* ant_hyper);

/* Ants' solution construction */

void ant_empty_memory( ant_struct *a, int n );

void place_ant( ant_struct *a , long int phase, int n, long int *seed);

void choose_best_next( ant_struct *a, long int phase, int n, struct ant_hyper_parameters* ant_hyper);

void neighbour_choose_best_next( ant_struct *a, long int phase, struct problem *instance, int n, struct ant_hyper_parameters* ant_hyper);

void choose_closest_next( ant_struct *a, long int phase, struct problem* instance, int n);

void neighbour_choose_and_move_to_next( ant_struct *a , long int phase, struct problem *instance, int n, struct ant_hyper_parameters* ant_hyper, long int* seed);

/* Auxiliary procedures related to ants */

long int find_best ( struct ant_hyper_parameters* ant_hyper);

long int find_worst( struct ant_hyper_parameters* ant_hyper);

void copy_from_to(ant_struct *a1, ant_struct *a2, int n);

void allocate_ants ( int n, struct ant_hyper_parameters* ant_hyper);

long int nn_tour(struct problem *instance, int n, struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper, struct ls_parameters* ls_params, long int* seed);

long int distance_between_ants( ant_struct *a1, ant_struct *a2, int n);

/* Procedures specific to MAX-MIN Ant System */

void mmas_evaporation_nn_list(struct problem* instance, int n, struct ant_hyper_parameters* ant_hyper);

void check_nn_list_pheromone_trail_limits(struct problem *instance, int n, struct ant_hyper_parameters* ant_hyper);

void check_pheromone_trail_limits( int n, struct ant_hyper_parameters* ant_hyper);

/* Procedures specific to Ant Colony System */

void global_acs_pheromone_update( ant_struct *a, struct problem* instance, int n, struct ant_hyper_parameters* ant_hyper);

void local_acs_pheromone_update( ant_struct *a, long int phase, struct problem* instance, int n, struct ant_hyper_parameters* ant_hyper);

/* Procedures specific to Best Worst Ant System */

void bwas_worst_ant_update( ant_struct *a1, ant_struct *a2, int n, struct ant_hyper_parameters* ant_hyper);

void bwas_pheromone_mutation( int n, struct hyper_parameters* hyper, struct ant_hyper_parameters *ant_hyper, long int* seed);
