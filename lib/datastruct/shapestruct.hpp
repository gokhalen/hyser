/*
  A shape function data structure. This is passed to the integration and, shape function routines which then store shape function values, jacobian etc in it. 
  Nachiket Gokhale  gokhalen@bu.edu 
 */

#ifndef HAVE_SHAPESTRUCT_HPP
#define HAVE_SHAPESTRUCT_HPP
#endif

#ifndef HAVE_MEMUTILS_HPP
#include "../utillib/memutils.hpp"
#endif

#include<iostream>
#include<new>

using namespace std;
class shapestruct{
public:

  int eltype;           // element type
  long int elemno;      // global element number 

  int ndime;            // dimension of the problem
  int selfdime;         // selfdime is the dimension of the calling element 
                        // for example a line can live in 3D but it itself is 1D                        // this is for the integration rules
                        // default == ndime
  int numnodes;         // number of nodes

  int *ninte;           // number of integration points in each direction size:(ndime);
  int tinte;            // total number of integration points

  double *jdet;         // jacobian determinant (tinte)
  double ***jacoinv;    //
  double ***jacomat;    // jacobian matrix: linear TJRH pg.147 eqn 3.9.6
                        // (tinte,ndime(global coord),ndime(parent coord der)) 
  double **coord;       // coordinates of nodes (numnodes,ndime)
  double **intpts;      // natural coordinates of integration points (tinte,ndime)
  double **glbcrdint;   // global coordinates of integration points (tinte,ndime)
  double **shape;       // value of shape functions: a 2D array : (tinte,numnodes)
  double ***glbder;     // global derivative of shape function: a 2D array : (tinte,numnodes,ndime)
  double ***natder;     // natural derivative of shape function: (tinte,numnodes,ndime)
  double *wts;          // weights of integration points: (tinte)

  double **wishape;     // wilson's incompatible shape functions (iinte,dime) 
  double ***winatder;   // natural derivatives of wilsons shape functions (iinte,dime, dime)
  double ***wiglbder;      // global derivatives of wilsons shape functions (iinte,dime,dime)


public:
  shapestruct();
  ~shapestruct();
  shapestruct(int cndime,int selfdime,int *cinte,int cnumnodes,int celtype,long int  celemno); // constructor
  void viewgauss();   //view the data set by integration library
  void viewshape();   //view the data set by shapefunction library
  void shapealloc(int *cinte); //allocates memory for shape datastructure, sets number of integration pts in each direction

};


