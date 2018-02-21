/*-----------------------------------------------------------------------------
 * Implementation of library mylib.h for BarnesHut.c
 *
 * Implementation of a Barnes-Hut Quadtree algorithm
 * 
 * Part of assignment 1 in the course 4EM30:
 *   Scientific Computing for Mechanical Engineering
 *   2017-2018
 *
 * (c) 2018 Joris Remmers TU/e
-----------------------------------------------------------------------------*/


#include <stdlib.h>
#include "mylib.h"

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


void readInput
  
  ( char*           name  ,
    BodyList*       blist )
  
{
  FILE         *fp;
  float        x,y,vx,vy,mass,domainSize;
  int          iBod,nBod;
  Vector       pos,velo;

  clearBodyList( blist );
  
  if ( ( fp=fopen(name,"r") ) == NULL) 
  {
    printf("Cannot open file.\n");
  }
  
  fscanf( fp, "%d", &nBod );
  fscanf( fp, "%e", &domainSize );

  blist->domainSize = domainSize;
  
  for ( iBod = 0 ; iBod < nBod ; iBod++ )
  {    
    fscanf(fp,"%e %e %e %e %e",&x,&y,&vx,&vy,&mass);
    
    pos.x  = x;
    pos.y  = y;
    velo.x = vx;
    velo.y = vy;

    addBody( blist , pos , velo , mass );
  }

  fclose( fp );
}


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


int addBody

  ( BodyList   *bl  ,
    Vector     pos  ,
    Vector     velo ,
    double     mass )

{
  bl->body[bl->nBod].pos.x  = pos.x;
  bl->body[bl->nBod].pos.y  = pos.y;

  bl->body[bl->nBod].velo.x = velo.x;
  bl->body[bl->nBod].velo.y = velo.y;

  bl->body[bl->nBod].mass   = mass;
  bl->body[bl->nBod].idx    = bl->nBod;

  bl->nBod += 1;

  return bl->nBod-1;
}


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


void clearBodyList
  
  ( BodyList   *bl )

{
  bl->nBod = 0;
}


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


void printBodies
  
  ( BodyList*       blist )

{
  int iBod;

  for ( iBod = 0 ; iBod < blist->nBod ; iBod++ )
  {
    printf("Body %d; Pos.: %e %e; Vel.: %e %e; Mass: %e\n",iBod,
                                                           blist->body[iBod].pos.x,
                                                           blist->body[iBod].pos.y,
                                                           blist->body[iBod].velo.x,
                                                           blist->body[iBod].velo.y,
                                                           blist->body[iBod].mass );
  }
}


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


void resetQuadTree

  ( QuadTree*       quadtree )

{
  quadtree->nNod = 0;
}


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


void updateCoM

  ( Node*           node ,
    Body*           body )

{
  node->com.x = node->com.x*node->mass + body->pos.x*body->mass;
  node->com.y = node->com.y*node->mass + body->pos.y*body->mass;
  node->mass  = node->mass + body->mass;

  node->com.x = node->com.x/node->mass;
  node->com.y = node->com.y/node->mass;
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


int getChild

  ( QuadTree*       qt   ,
    Vector          pos  ,
    Node*           node )

{
  double cx = 0.5*( node->box.point2.x+node->box.point1.x );
  double cy = 0.5*( node->box.point2.y+node->box.point1.y );

  int childIdx = -1;
  int idx      = -1;

  if ( pos.x < cx )
  {
    if ( pos.y < cy )
    {
      childIdx = 0;
    }
    else
    {
      childIdx = 3;
    }
  }
  else
  {
    if ( pos.y < cy )
    {
      childIdx = 1;
    }
    else
    {
      childIdx = 2;
    }
  }

  idx = node->child[childIdx];

  if ( idx == -1 )   // Node does not exist
  {
    idx = initNode( qt , childIdx , &node->box );
    node->child[childIdx] = idx;
  }

  return idx;
}


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


void addBodyToNode

  ( QuadTree*       quadtree ,
    Body*           body     ,
    int             idx      )

{ 
  int newidx = -1;
   
  if ( quadtree->node[idx].nBod == 0 ) // Store body in this node
  {
    updateCoM( &quadtree->node[idx] , body );
    quadtree->node[idx].nBod  = 1;
    quadtree->node[idx].idx   = body->idx;
  }
  else if( quadtree->node[idx].nBod == 1 ) // Node was leaf, will be internal
  {
    Body body0;

    body0.pos.x = quadtree->node[idx].com.x;
    body0.pos.y = quadtree->node[idx].com.y;
    body0.mass  = quadtree->node[idx].mass;
    body0.idx   = quadtree->node[idx].idx;

    newidx = getChild( quadtree , body0.pos , &quadtree->node[idx] );
    addBodyToNode( quadtree , &body0 , newidx );

    newidx = getChild( quadtree , body->pos , &quadtree->node[idx] );
    addBodyToNode( quadtree , body , newidx );
    
    updateCoM( &quadtree->node[idx] , body );

    quadtree->node[idx].nBod++;
  }
  else // Node is internal
  {
    newidx = getChild( quadtree , body->pos , &quadtree->node[idx] );
    addBodyToNode( quadtree , body , newidx );
    
    updateCoM( &quadtree->node[idx] , body );
  }
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


void initQuadTree

  ( QuadTree*    qt , 
    double       length )

{
  qt->nNod = 1;

  qt->node[0].box.point1.x = -0.5*length;
  qt->node[0].box.point1.y = -0.5*length;
  qt->node[0].box.point2.x =  0.5*length;
  qt->node[0].box.point2.y =  0.5*length;
  
  qt->node[0].nBod = 0;

  qt->node[0].child[0] = -1;
  qt->node[0].child[1] = -1;
  qt->node[0].child[2] = -1;
  qt->node[0].child[3] = -1;
}


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


int initNode

  ( QuadTree*    qt       , 
    int          childIdx , 
    Box*         box      )

{
  const int idx = qt->nNod;

  qt->node[idx].box.point1.x = box->point1.x;
  qt->node[idx].box.point1.y = box->point1.y;
  qt->node[idx].box.point2.x = box->point2.x;
  qt->node[idx].box.point2.y = box->point2.y;

  Vector center;

  center.x = 0.5f*(box->point1.x+box->point2.x);
  center.y = 0.5f*(box->point1.y+box->point2.y);

  if ( childIdx == 0 )
  {
    qt->node[idx].box.point2.x = center.x;
    qt->node[idx].box.point2.y = center.y;
  }
  else if ( childIdx == 1 )
  {
    qt->node[idx].box.point1.x = center.x;
    qt->node[idx].box.point2.y = center.y;
  }
  else if ( childIdx == 2 )
  {
    qt->node[idx].box.point1.x = center.x;
    qt->node[idx].box.point1.y = center.y;
  }
  else if ( childIdx == 3 )
  {
    qt->node[idx].box.point1.y = center.y;
    qt->node[idx].box.point2.x = center.x;
  }

  qt->node[idx].child[0] = -1;
  qt->node[idx].child[1] = -1;
  qt->node[idx].child[2] = -1;
  qt->node[idx].child[3] = -1;
  qt->node[idx].nBod     =  0;

  qt->node[idx].com.x    = 0.;
  qt->node[idx].com.y    = 0.;
  qt->node[idx].mass     = 0.;
  qt->node[idx].idx      = idx;
 
  qt->nNod += 1;

  return idx;
}


//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------


void printQuadTree

  ( QuadTree*   qt )

{
  int iNod;

  for ( iNod = 0 ; iNod < qt->nNod ; iNod++ )
  {
    printf("Node %d, Box %e %e %e %e\n ",iNod,qt->node[iNod].box.point1.x ,
                                              qt->node[iNod].box.point1.y ,
                                              qt->node[iNod].box.point2.x ,
                                              qt->node[iNod].box.point2.y );
  }

  printf("The number of nodes is: %d\n",qt->nNod);
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

void GenerateXMLfile

( BodyList*       blist,
  QuadTree*       qt )

{

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

    for ( int iBod = 0 ; iBod < blist->nBod ; iBod++ )
        {
        fprintf(f, "<circle cx=\"%f\" cy=\"%f\" r=\"0.9\" fill=\"red\"/> \n" , 10+(blist->body[iBod].pos.x+qt->node[0].box.point2.x)*10, 500+(((blist->body[iBod].pos.y)*-1)+qt->node[0].box.point2.y)*10);
        }

   
//end document
fprintf(f, "</svg> \n");

fclose(f);

}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

void clearForces

  ( BodyList*   blist)

{
  int iBod;

  for ( iBod = 0 ; iBod < blist->nBod ; iBod++ )
  {
    blist->body[iBod].force.x = 0;
    blist->body[iBod].force.y = 0;
  }
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

void bruteForces

  ( BodyList*   blist)

{
  int iBod;

  for ( iBod = 0 ; iBod < blist->nBod ; iBod++ )
  {
    bruteForceBody( blist, iBod );
  }
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

void bruteForceBody

  ( BodyList*   blist ,
    int         iBod  )

{
  Vector  posi;
  Vector  posj;
  Vector  unitij;
  Vector  forceVec;
  Vector  newForcei;
  Vector  newForcej;
  double  mi;
  double  mj;
  int     jBod;
  double  distance;
  double  force;

  posi    = blist->body[iBod].pos;
  mi      = blist->body[iBod].mass;

  for ( jBod = (iBod + 1) ; jBod < blist->nBod ; jBod++ )
  {
    posj          = blist->body[jBod].pos;
    mj            = blist->body[jBod].mass;

    distance      = sqrt(pow((posi.x - posj.x), 2) + pow((posi.y - posj.y), 2));
    force         = GRAV_CONSTANT * mi * mj / pow(distance, 2);

    unitij.x      = (posj.x - posi.x) / distance;
    unitij.y      = (posj.y - posi.y) / distance;
    forceVec.x    = force * unitij.x;
    forceVec.y    = force * unitij.y;

    newForcei.x   = blist->body[iBod].force.x + forceVec.x;
    newForcei.y   = blist->body[iBod].force.y + forceVec.y;
    newForcej.x   = blist->body[jBod].force.x - forceVec.x;
    newForcej.y   = blist->body[jBod].force.y - forceVec.y;

    blist->body[iBod].force = newForcei;
    blist->body[jBod].force = newForcej;
  }
}

//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------

void printBruteForces

  ( BodyList* blist)

{
  int iBod;

  for ( iBod = 0 ; iBod < blist->nBod ; iBod++ )
  {
    printf("The resulting brute force on body %d is %e in x direction and %e in y direction\n", iBod, blist->body[iBod].force.x, blist->body[iBod].force.y);
  }
}
