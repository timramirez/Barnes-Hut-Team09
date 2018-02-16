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
//bla
{  
  time_t       t0,t1;
     
  BodyList     *blist = malloc(sizeof( BodyList ) + 100000*sizeof(Body));
  QuadTree     *qt    = malloc(sizeof( QuadTree ) + 200000*sizeof(Node));

  readInput( "input50.txt" , blist );
  
  printBodies( blist );

  t0 = clock();
  
  initQuadTree( qt , blist->domainSize );

  for ( int i = 0 ; i < blist->nBod ; i++ )
  {
    addBodyToNode( qt , &blist->body[i] , 0 );
  }

  t1 = clock();

  printf("Time needed to fill the QuadTree:   %f  seconds.\n", 
           (double)(t1 - t0)/CLOCKS_PER_SEC );
 
  return 0;
}

