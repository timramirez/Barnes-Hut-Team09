#ifndef MYLIB_H
#define MYLIB_H

/*-----------------------------------------------------------------------------
 * Header file for Library mylib.h for BarnesHut.c
 *
 * Implementation of a Barnes-Hut Quadtree algorithm
 * 
 * Part of assignment 1 in the course 4EM30:
 *   Scientific Computing for Mechanical Engineering
 *   2017-2018
 *
 * (c) 2018 Joris Remmers TU/e
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
  Vector      force;
  double      mass;
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
//    pre    : the file should exist
//    post   : the bodylist blist is filled
//    result : -
//-----------------------------------------------------------------------------

void readInput
  
  ( char*           name  ,
    BodyList*       blist );

//-----------------------------------------------------------------------------
//  addBody  : Function that stores a new body in the BodyList
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
//  clearBodyList  : Function to empty the BodyList
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
//  addBodyToNode : Function that adds a body to node idx
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
//  initQuadTree : initialises the QuadTree
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
//  printQuadTree: Function that prints the quadtree
//    pre    : The QuadTRee must be initialised. 
//    post   : The quadtree is printed to stdout.
//    result : -
//-----------------------------------------------------------------------------

void printQuadTree

  ( QuadTree*   qt );

//-----------------------------------------------------------------------------
//  GenerateXMLfile: Function that creates a XML file containing the position 
//                   of the bodies and the structure of the Barnes-Hut nodes
//    pre    : The QuadTRee and the BodyList must be initialised. 
//    post   : The XML file is created in the main directory.
//    result : A file called 'file.xml' (stored in the directory of all other files
//             which can be opened in e.g. a Chrome Webbroser.
//-----------------------------------------------------------------------------


void GenerateXMLfile

( BodyList*       blist,
  QuadTree*       qt );

//  clearForces: Function that clears the forces on all bodies 
//    pre    : -
//    post   : The forces in all bodies is reset to zero.
//    result : -
//-----------------------------------------------------------------------------

void clearForces

  ( BodyList*  blist );

//-----------------------------------------------------------------------------
//  bruteForces: Function that calculates the brute forces on all bodies 
//    pre    : -
//    post   : The brute forces are calculated and saved in all bodies.
//    result : -
//-----------------------------------------------------------------------------

void bruteForces

  ( BodyList*  blist );

//-----------------------------------------------------------------------------
//  bruteForceBody: Function that calculates the brute forces on a single body
//    pre    : -
//    post   : The brute forces are calculated and saved in one body.
//    result : -
//-----------------------------------------------------------------------------

void bruteForceBody

  ( BodyList*   blist ,
    int         iBod  );

//-----------------------------------------------------------------------------
//  printBruteForces: Function that prints the brute forces on all bodies
//    pre    : -
//    post   : The brute forces on all bodies are printed.
//    result : -
//-----------------------------------------------------------------------------

void printBruteForces

  ( BodyList*   blist );


#endif
