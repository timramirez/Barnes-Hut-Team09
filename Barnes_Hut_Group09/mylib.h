#ifndef MYLIB_H
#define MYLIB_H

/*-----------------------------------------------------------------------------
 * Header file for Library mylib.h for BarnesHut.c
 *
 * Implementation of a Barnes-Hut Quadtree algorithm
 * 
 * Assignment 1 in the course 4EM30 group 9:
 *   Scientific Computing for Mechanical Engineering
 *   2017-2018
 *
 * (c) 2018 Joris Remmers, Michael Geurtsen and Tim Ramirez TU/e
-----------------------------------------------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#define GRAV_CONSTANT 1

//=============================================================================
//  Structures
//=============================================================================

//-----------------------------------------------------------------------------
//  Vector(2d)
//-----------------------------------------------------------------------------

typedef struct
{
  double      x;
  double      y;
} Vector;

//-----------------------------------------------------------------------------
//  Body, data type that contains the current position, velocity, force, mass and
//    index of a celestial body.
//-----------------------------------------------------------------------------

typedef struct 
{
  Vector      pos;
  Vector      velo;
  Vector      bruteForce;
  Vector      barnesHutForce;
  double      mass;
  double      forceErrori;
  int         idx;
} Body;

//-----------------------------------------------------------------------------
//  BodyList: datatype that contains a list of bodies, the length of the 
//    original domain and the number of bodys. Note that this structure has 
//    a dynamic length.
//-----------------------------------------------------------------------------

typedef struct
{
  double      domainSize;
  double      forceError;
  int         nBod;
  Body        body[];
} BodyList;

//-----------------------------------------------------------------------------
//  Box: data structure that describes the position of a box by its lef bottom
//    corner (point1) and the top right corner (point2).
//-----------------------------------------------------------------------------

typedef struct
{
  Vector      point1;
  Vector      point2;
} Box;

//-----------------------------------------------------------------------------
//  Node: data structure that contains the data of a single node: number of 
//    bodies in the node, the centre of mass, the total mass, the node index
//    the children of this node and the coordinates.
//-----------------------------------------------------------------------------

typedef struct
{
  int         nBod;
  Vector      com;
  double      mass;
  int         idx;
  int         child[4];
  Box         box;
} Node;

//-----------------------------------------------------------------------------
//  QuadTree: datastructure that is basically an array of nodes and a
//   counter to determine the number of active nodes.
//-----------------------------------------------------------------------------


typedef struct
{
  int         nNod;
  Node        node[];
} QuadTree;

//=============================================================================
//  Functions
//=============================================================================

//-----------------------------------------------------------------------------
//  readInput  : Function that reads the body data, position, velocities and 
//               mass from the input file.
//             The file has the following format
//             50   (Number of bodies)
//             1.0  (dimension of domain)
//             2.0 1.0 0.1 0.2 5.0  (x and y position, x and y velocity and mass of the body
//             .....etc.
//    pre    : the file should exist.
//    post   : the bodylist blist is filled.
//    result : -
//-----------------------------------------------------------------------------

void readInput
  
  ( char*           name  ,
    BodyList*       blist );

//-----------------------------------------------------------------------------
//  addBody  : Function that stores a new body in the BodyList.
//    pre    : the structure BodyList is not full.
//    post   : a new body is stored in the next available position.
//    result : the body index number is returned. 
//-----------------------------------------------------------------------------

int addBody

  ( BodyList   *bl  ,
    Vector     pos  ,
    Vector     velo ,
    double     mass );

//-----------------------------------------------------------------------------
//  clearBodyList  : Function to empty the BodyList.
//    pre    : -
//    post   : the counter nBod is set to zero.
//    result : -
//-----------------------------------------------------------------------------

void clearBodyList
  
  ( BodyList   *bl );

//-----------------------------------------------------------------------------
//  printBodies  : Function to print the bodylist to the screen.
//    pre    : -
//    post   : The bodylist is printed to the screen.
//    result : -
//-----------------------------------------------------------------------------

void printBodies
  
  ( BodyList*       blist );

//-----------------------------------------------------------------------------
//  resetQuadTree  : Function that empties the quadtree.
//    pre    : -
//    post   : the counter nNod is set to zero.
//    result : -
//-----------------------------------------------------------------------------

void resetQuadTree

  ( QuadTree*       quadtree );

//-----------------------------------------------------------------------------
//  updateCoM  : Function to update the centre of mass in a node.
//    pre    : node and body should be existing. Node may be empty.
//    post   : the new centre of mass is store in node.
//    result : -
//-----------------------------------------------------------------------------

void updateCoM

  ( Node*           node ,
    Body*           body );

//-----------------------------------------------------------------------------
//  getChild  : Function that determines in which child a body is.
//    pre    : The QuadTRee must be initialised. pos is the psoition of the 
//               body under investigation and node is the current node.
//    post   : the exact child node id is stored in node.child[ ]. If the 
//               child ndoe does not exist, it is created.
//    result : the global node index number is returned. 
//-----------------------------------------------------------------------------

int getChild

  ( QuadTree*       qt   ,
    Vector          pos  ,
    Node*           node );

//-----------------------------------------------------------------------------
//  addBodyToNode : Function that adds a body to node idx.
//    pre    : The QuadTRee must be initialised. 
//    post   : the body is added to the node idx. IF this node is an
//               internal node, this function is called recursively.
//    result : -
//-----------------------------------------------------------------------------

void addBodyToNode

  ( QuadTree*       quadtree ,
    Body*           body     ,
    int             idx      );

//-----------------------------------------------------------------------------
//  initQuadTree : initialises the QuadTree.
//    pre    : The QuadTRee must be initialised. 
//    post   : the body is added to the node idx. IF this node is an
//               internal node, this function is called recursively.
//    result : -
//-----------------------------------------------------------------------------

void initQuadTree

  ( QuadTree*    qt , 
    double       length );

//-----------------------------------------------------------------------------
//  initNode : Function that initialises a new node in the quadtree.
//    pre    : Box is the size of the mother node, childIdx the child index.
//    post   : the node is added to the quadtree. 
//    result : the global node idx number.
//-----------------------------------------------------------------------------

int initNode

  ( QuadTree*    qt       , 
    int          childIdx , 
    Box*         box      );

//-----------------------------------------------------------------------------
//  printQuadTree: Function that prints the quadtree.
//    pre    : The QuadTRee must be initialised. 
//    post   : The quadtree is printed to stdout.
//    result : -
//-----------------------------------------------------------------------------

void printQuadTree

  ( QuadTree*   qt );

//-----------------------------------------------------------------------------
//  GenerateXMLfile: Function that creates a XML file containing the position. 
//                   of the bodies and the structure of the Barnes-Hut nodes.
//    pre    : The QuadTRee and the BodyList must be initialised. 
//    post   : The XML file named 'BarnesHutTree.xml' is created in the main directory of the project.
//    result : -
//-----------------------------------------------------------------------------


void GenerateXMLfile

  ( QuadTree*       qt    , 
    BodyList*       blist );

//-----------------------------------------------------------------------------
//  clearBruteForces: Function that clears the forces on all bodies.
//    pre    : -
//    post   : The forces in all bodies is reset to zero.
//    result : -
//-----------------------------------------------------------------------------

void clearBruteForces

  ( BodyList*  blist );

//-----------------------------------------------------------------------------
//  bruteForces: Function that calculates the brute forces on all bodies.
//    pre    : The forces have to be cleared first.
//    post   : The brute forces are calculated and saved in all bodies.
//    result : -
//-----------------------------------------------------------------------------

void bruteForces

  ( BodyList*  blist );

//-----------------------------------------------------------------------------
//  bruteForceBody: Function that calculates the brute forces on a single body.
//    pre    : -
//    post   : The brute forces are calculated and saved in the body that is has
//             been called for.
//    result : -
//-----------------------------------------------------------------------------

void bruteForceBody

  ( BodyList*   blist ,
    int         iBod  );

//-----------------------------------------------------------------------------
//  printBruteForces: Function that prints the brute forces in x and y on all bodies.
//    pre    : The brute forces have to be initialised first.
//    post   : The brute forces on all bodies are printed.
//    result : -
//-----------------------------------------------------------------------------

void printForces

  ( BodyList*   blist );

//-----------------------------------------------------------------------------
//  barnesHut: Function that calculates the BarnesHut forces on the system.
//    pre    : The foreces have to be cleared first.
//    post   : The forces on the system are calculated via the BarnesHut algorithm.
//    result : -
//-----------------------------------------------------------------------------

void barnesHut

  ( QuadTree*       qt    , 
    BodyList*       blist ,
    double          theta );

//-----------------------------------------------------------------------------
//  barnesHutBody: Function that makes the dicision if the forces on a body caused by
//                 a node may be calculated or the children of the node have to be evaluated.
//    pre    : The quadtree has to initialised.
//    post   : The oby and node are passed on to the calculating function or the children of 
//             the node are evaluated in this function again.
//    result : -
//-----------------------------------------------------------------------------

void barnesHutBody

  ( QuadTree*       qt    , 
    BodyList*       blist ,
    int             iBod  ,
    double          theta , 
    int             iNod  );

//-----------------------------------------------------------------------------
//  checkCalculationMethod: Function that determines what step in the BarnesHut 
//                          algorithm should be made next.
//    pre    : -
//    post   : A step in the BarnesHut algorithm is selected.
//    result : 0 for a leaf, 1 if the center of mass of the node can be used for
//             the calculation and 2 if the children of the node has to be evaluated.
//-----------------------------------------------------------------------------

int checkCalculationMethod

  ( QuadTree*       qt    , 
    BodyList*       blist ,
    int             iBod  ,
    double          theta ,
    int             iNod  );

//-----------------------------------------------------------------------------
//  forceBarnesHut: Function that calculates and updates the Barnes-Hut force on a 
//                  single body casued by a single node in the quad tree.
//    pre    : -
//    post   : The force on a body is updated for the interaction with the node.
//    result : -
//-----------------------------------------------------------------------------

void forceBarnesHut

  ( QuadTree*       qt    ,
    BodyList*       blist ,
    int             iNod  ,
    int             iBod  );

//-----------------------------------------------------------------------------
//  clearBarnesHut: A function that clears the Barnes-Hut forces on all bodies.
//    pre    : -
//    post   : The force in x and y direction of all bodies is set to zero.
//    result : -
//-----------------------------------------------------------------------------

void clearBarnesHut

  ( BodyList*       blist );

//-----------------------------------------------------------------------------
//  error: A function that determines the average error of the Barnes-Hut forces 
//         compared to the brute forces.
//    pre    : The brute forces and Barnes-Hut forces have to be initialised. 
//    post   : The error of the Barnes-Hut force is printed to the screen.
//    result : -
//-----------------------------------------------------------------------------

void error

  ( BodyList*       blist ,
    double          theta );

//-----------------------------------------------------------------------------
//  resetError: A function that resets the average error to zero.
//    pre    : -
//    post   : The average error is set to zero.
//    result : -
//-----------------------------------------------------------------------------

void resetError

  ( BodyList*     blist );

//-----------------------------------------------------------------------------
//  maxError: A function that determines the maximum error on a single body and 
//            prints it to the screen.
//    pre    : -
//    post   : The maximum error on a body is printed to the screen.
//    result : -
//-----------------------------------------------------------------------------

void maxError

  ( BodyList*     blist ,
    double        theta );

#endif
