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
      File:    ls.h
      Author:  Thomas Stuetzle
      Purpose: header file for local search routines
      Check:   README and gpl.txt
      Copyright (C) 1999  Thomas Stuetzle
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

struct ls_parameters {
    long int ls_flag; /* indicates whether and which local search is used */ 

    long int nn_ls; /* maximal depth of nearest neighbour lists used in the   local search */ 

    long int dlb_flag;  /* flag indicating whether don't look bits are used. I recommend to always use it if local search is applied */
} ;


void two_opt_first( long int *tour, struct problem *instance, int n, struct ls_parameters* ls_params, long int* seed);

void two_h_opt_first( long int *tour, struct problem *instance, int n, long int* seed);

void three_opt_first( long int *tour, struct problem *instance, int n, struct ls_parameters* ls_params, long int* seed);
