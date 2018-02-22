/*-----------------------------------------------------------------------------
 * BarnesHut.c
 *
 * Implementation of a Barnes-Hut Quadtree algorithm
 * 
 * Part of assignment 1 in the course 4EM30:
 *   Scientific Computing for Mechanical Engineering
 *   2017-2018
 *
 * (c) 2018 Joris Remmers TU/e
-----------------------------------------------------------------------------*/

#include "mylib.h"
#include <time.h>

int main( void )

{  
  time_t       t0,t1,t2,t3,t4,t5;
     
  BodyList     *blist = malloc(sizeof( BodyList ) + 100000*sizeof(Body));
  QuadTree     *qt    = malloc(sizeof( QuadTree ) + 200000*sizeof(Node));

  readInput( "input5000.txt" , blist );
  
//  printBodies( blist );
  
  t4 = clock();

  initQuadTree( qt , blist->domainSize );

  for ( int i = 0 ; i < blist->nBod ; i++ )
  {
    addBodyToNode( qt , &blist->body[i] , 0 );
  }

  t5 = clock();

  printf("Time needed to build the QuadTree:   %f  seconds.\n", 
           (double)(t5 - t4)/CLOCKS_PER_SEC );

  clearBruteForces( blist );

  t0 = clock();

  bruteForces( blist );

  t1 = clock();

  printf("Time needed to calculate brute forces:   %f  seconds.\n", 
           (double)(t1 - t0)/CLOCKS_PER_SEC );

//  printQuadTree(qt);

  GenerateXMLfile(qt, blist);

  double  theta;

  for ( theta = 0.0 ; theta <= 1.0 ; theta+=0.1 )
  {
    t2 = clock();

    clearBarnesHut( blist );

    barnesHut(qt, blist, theta);

    t3 = clock();

    printf("Time needed to calculate Barnes-Hut forces:   %f  seconds.\n", 
           (double)(t3 - t2)/CLOCKS_PER_SEC );

    error ( blist );
  }

//  printForces ( blist );

  writeForcesToText ( blist );

//  testPrint( qt, blist );

  return 0;
}

