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

  printQuadTree(qt);









// Plot of the position of bodies and the structure of the Barnes-Hut nodes (SVG)

FILE *f = fopen("file.xml", "w");
if (f == NULL)
{
    printf("Error opening file!\n");
    exit(1);
}

//print initialisation lines
fprintf(f, "<?xml version=\"1.0\" standalone=\"no\"?> \n");
fprintf(f, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n");
fprintf(f, "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\"> \n");
fprintf(f, "<svg width=\"1000px\" height=\"1000px\" version=\"1.1\" \n");
fprintf(f, "xmlns=\"http://www.w3.org/2000/svg\"> \n");

//Generating position plot
fprintf(f, "<line x1='220' y1='20' x2='220' y2='420' style='stroke:rgb(0,0,0);stroke-width:0.25'/> \n");
fprintf(f, "<line x1='20' y1='220' x2='420' y2='220' style='stroke:rgb(0,0,0);stroke-width:0.25'/> \n");

  for ( int iBod = 0 ; iBod < blist->nBod ; iBod++ )
  {
    fprintf(f, "<circle cx=\"%f\" cy=\"%f\" r=\"1\" fill=\"red\"/> \n" , (blist->body[iBod].pos.x+20)*10, (blist->body[iBod].pos.y+20)*10);
  }
   



//end document
fprintf(f, "</svg> \n");

fclose(f);

  return 0;
}

