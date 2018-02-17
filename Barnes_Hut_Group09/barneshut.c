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

  readInput( "input500.txt" , blist );
  
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
fprintf(f, "<line x1='200' y1='0' x2='200' y2='400' style='stroke:rgb(0,0,0);stroke-width:0.25'/> \n");
fprintf(f, "<line x1='0' y1='200' x2='400' y2='200' style='stroke:rgb(0,0,0);stroke-width:0.25'/> \n");

    for (int i=0 ; i <= qt->node[0].box.point2.x*2 ; i++)
        {
        fprintf(f, "<line x1='%d' y1='0' x2='%d' y2='400' style='stroke:rgb(0,0,0);stroke-width:0.10'/> \n", i*10, i*10);
        fprintf(f, "<line x1='0' y1='%d' x2='400' y2='%d' style='stroke:rgb(0,0,0);stroke-width:0.10'/> \n", i*10, i*10);
        }

    for (int i=0 ; i <= qt->node[0].box.point2.x*2 ; i++)
        {
        fprintf(f, "<text font-size=\"4\"  x=\"%d\" y=\"206\"> \n", i*10);
        fprintf(f, " %d \n", i-20);
        fprintf(f, "</text>\n");
        }

    for (int i=0 ; i <= qt->node[0].box.point2.x*2 ; i++)
        {
        fprintf(f, "<text font-size=\"4\"  x=\"194\" y=\"%d\"> \n", 400-i*10);
        fprintf(f, " %d \n", i-20);
        fprintf(f, "</text>\n");
        }

    for ( int iBod = 0 ; iBod < blist->nBod ; iBod++ )
        {
        fprintf(f, "<circle cx=\"%f\" cy=\"%f\" r=\"0.9\" fill=\"red\"/> \n" , (blist->body[iBod].pos.x+qt->node[0].box.point2.x)*10, (((blist->body[iBod].pos.y)*-1)+qt->node[0].box.point2.y)*10);
        }

//Generating Barnes-Hut plot
    for ( int iNod = 0 ; iNod < qt->nNod ; iNod++ )
        {
        fprintf(f, "<line x1='%e' y1='%e' x2='%e' y2='%e' style='stroke:rgb(0,0,0);stroke-width:0.10'/> \n", 10+((qt->node[iNod].box.point1.x)+20)*10 , 500+(((qt->node[iNod].box.point1.y)*-1)+20)*10 , 10+((qt->node[iNod].box.point1.x)+20)*10 , 500+(((qt->node[iNod].box.point2.y)*-1)+20)*10 );
        fprintf(f, "<line x1='%e' y1='%e' x2='%e' y2='%e' style='stroke:rgb(0,0,0);stroke-width:0.10'/> \n", 10+((qt->node[iNod].box.point1.x)+20)*10 , 500+(((qt->node[iNod].box.point1.y)*-1)+20)*10 , 10+((qt->node[iNod].box.point2.x)+20)*10 , 500+(((qt->node[iNod].box.point1.y)*-1)+20)*10 );
        fprintf(f, "<line x1='%e' y1='%e' x2='%e' y2='%e' style='stroke:rgb(0,0,0);stroke-width:0.10'/> \n", 10+((qt->node[iNod].box.point2.x)+20)*10 , 500+(((qt->node[iNod].box.point2.y)*-1)+20)*10 , 10+((qt->node[iNod].box.point1.x)+20)*10 , 500+(((qt->node[iNod].box.point2.y)*-1)+20)*10 );
        fprintf(f, "<line x1='%e' y1='%e' x2='%e' y2='%e' style='stroke:rgb(0,0,0);stroke-width:0.10'/> \n", 10+((qt->node[iNod].box.point2.x)+20)*10 , 500+(((qt->node[iNod].box.point2.y)*-1)+20)*10 , 10+((qt->node[iNod].box.point2.x)+20)*10 , 500+(((qt->node[iNod].box.point1.y)*-1)+20)*10 );
        }


   
//end document
fprintf(f, "</svg> \n");

fclose(f);

  return 0;
}

