/*-----------------------------------------------------------------------------
 * BarnesHut.c
 *
 * Implementation of a Barnes-Hut Quadtree algorithm
 * 
 * Assignment 1 in the course 4EM30 group 9:
 *   Scientific Computing for Mechanical Engineering
 *   2017-2018
 *
 * (c) 2018 Joris Remmers, Michael Geurtsen and Tim Ramirez TU/e
-----------------------------------------------------------------------------*/

#include "mylib.h"
#include <time.h>

int main( void )

{  
  time_t       t0,t1,t2,t3,t4,t5;
  double  theta = 0.5;
     
  BodyList     *blist = malloc(sizeof( BodyList ) + 100000*sizeof(Body));
  QuadTree     *qt    = malloc(sizeof( QuadTree ) + 200000*sizeof(Node));

  readInput( "input50000.txt" , blist );
  
//  printBodies( blist );
  
  t0 = clock();

  initQuadTree( qt , blist->domainSize );

  for ( int i = 0 ; i < blist->nBod ; i++ )
  {
    addBodyToNode( qt , &blist->body[i] , 0 );
  }

  t1 = clock();

  printf("Time needed to build the QuadTree: %f seconds.\n", 
           (double)(t1 - t0)/CLOCKS_PER_SEC );

  clearBruteForces( blist );

  t2 = clock();

  bruteForces( blist );

  t3 = clock();

  printf("Time needed to calculate brute forces: %f seconds.\n", 
           (double)(t3 - t2)/CLOCKS_PER_SEC );

//  printQuadTree(qt);

  GenerateXMLfile(qt, blist);

  t4 = clock();

  clearBarnesHut( blist );

  barnesHut(qt, blist, theta);

  t5 = clock();

  printf("Time needed to calculate Barnes-Hut forces for theta %e: %f seconds.\n", 
         theta, (double)(t5 - t4)/CLOCKS_PER_SEC );

  error ( blist, theta );

  printf("Average error for theta %e: %f.\n", theta, blist->forceError);

//  printForces ( blist );

  return 0;
}

