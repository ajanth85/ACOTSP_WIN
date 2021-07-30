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
      File:    InOut.h
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



#define PROG_ID_STR      "\nACO algorithms for the TSP, v1.0\n"
#define CALL_SYNTAX_STR  "call syntax:  acotsp <param-list>\n"

#define LINE_BUF_LEN     100

struct point * read_etsp(const char *tsp_file_name, struct problem instance, int n);

struct point * read_etsp(const char *tsp_file_name, struct problem instance, int n);

struct thread_result {

    double std_iteration;
    double std_time;
    double std_tour;

    long int avg_best_iteration;
    double avg_best_time;
    long int avg_best_tour;
};

struct all_threads_result {
    double avg_std_iteration;
    double avg_std_time;
    double avg_std_tour;

    long int avg_best_iteration;
    double avg_best_time;
    long int avg_best_tour;
    long int number_of_threads;

};

struct hyper_parameters {
    long int* best_in_try;
    long int* best_found_at;
    double* time_best_found;
    double* time_total_run;

    long int n_try;        /* number of try */
    long int n_tours;      /* number of constructed tours */
    long int iteration;    /* iteration counter */
    long int restart_iteration;    /* iteration counter */
    double   restart_time; /* remember when restart was done if any */

    long int max_number_of_threads;
    long int max_tries;    /* maximum number of independent tries */
    long int max_tours;    /* maximum number of tour constructions in one try */

    double   lambda;       /* Parameter to determine branching factor */
    double   branch_fac;   /* If branching factor < branch_fac => update trails */

    double   max_time;     /* maximal allowed run time of a try  */
    double   time_used;    /* time used until some given event */
    double   time_passed;  /* time passed until some moment*/
    long int optimal;      /* optimal solution value or bound to find */

    double mean_ants;      /* average tour length */
    double stddev_ants;    /* stddev of tour lengths */
    double branching_factor; /* average node branching factor when searching */
    double found_branching;  /* branching factor when best solution is found */

    long int found_best;         /* iteration in which best solution is found */
    long int restart_found_best; /* iteration in which restart-best solution is found */

    FILE* report, * comp_report, * stat_report;

    char *name_buf;
    int  opt;
};


//extern long int *best_in_try;
//extern long int *best_found_at;
//extern double   *time_best_found;
//extern double   *time_total_run;   
//
//extern long int n_try;        /* number of try */
//extern long int n_tours;      /* number of constructed tours */
//extern long int iteration;    /* iteration counter */
//extern long int restart_iteration;    /* iteration counter */
//extern double   restart_time; /* remember when restart was done if any */
//
//extern long int max_tries;    /* maximum number of independent tries */
//extern long int max_tours;    /* maximum number of tour constructions in one try */
//
//extern double   lambda;       /* Parameter to determine branching factor */
//extern double   branch_fac;   /* If branching factor < branch_fac => update trails */
//
//extern double   max_time;     /* maximal allowed run time of a try  */
//extern double   time_used;    /* time used until some given event */
//extern double   time_passed;  /* time passed until some moment*/
//extern long int optimal;      /* optimal solution value or bound to find */
//
//extern double mean_ants;      /* average tour length */
//extern double stddev_ants;    /* stddev of tour lengths */
//extern double branching_factor; /* average node branching factor when searching */
//extern double found_branching;  /* branching factor when best solution is found */
//
//extern long int found_best;         /* iteration in which best solution is found */
//extern long int restart_found_best; /* iteration in which restart-best solution is found */
//
//extern FILE *report, *comp_report, *stat_report;
//
//extern char name_buf[LINE_BUF_LEN];
//extern int  opt;

void write_report( struct hyper_parameters *hyper, struct ant_hyper_parameters* ant_hyper );

void print_default_parameters(struct hyper_parameters* hyper, struct ant_hyper_parameters* ant_hyper, struct ls_parameters *ls_params);

void set_default_parameters(struct hyper_parameters* hyper, struct ant_hyper_parameters* ant_hyper, struct ls_parameters* ls_params);

void init_try( long int ntry, long int* seed);

void output_solution( struct problem *instance, int n, struct hyper_parameters* hyper, struct ant_hyper_parameters* ant_hyper);

void exit_try( long int ntry, struct problem *instance, int n, struct hyper_parameters* hyper, struct ant_hyper_parameters* ant_hyper, struct thread_result* thread_result);

void exit_program( struct problem *instance, struct hyper_parameters* hyper);

void init_program( long int argc, char *argv[], struct problem *instance, int *n, struct hyper_parameters* hyper, struct ant_hyper_parameters* ant_hyper, struct ls_parameters* ls_params, long int* seed);

void printDist(struct problem *instance, int n);

void printHeur(struct problem *instance, int n);

void printTrail(int n, struct hyper_parameters* hyper, struct ant_hyper_parameters* ant_hyper);

void printTotal(int n, struct ant_hyper_parameters* ant_hyper);

void printProbabilities(int n, struct hyper_parameters* hyper, struct ant_hyper_parameters* ant_hyper);

void printTour( long int *t, struct problem *instance, int n);

void printTourFile( long int *t , struct problem *instance, int n, struct hyper_parameters* hyper, struct ant_hyper_parameters* ant_hyper);

void checkTour( long int *t, struct problem *instance, int n);

void population_statistics ( struct problem *instance, int n, struct hyper_parameters* hyper, struct ant_hyper_parameters* ant_hyper);

double node_branching(double l, struct problem *instance, int n, struct ant_hyper_parameters* ant_hyper);

void write_params (struct hyper_parameters* hyper, struct ant_hyper_parameters* ant_hyper, struct ls_parameters* ls_params);
