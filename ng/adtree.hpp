#pragma once
/* *************************************************************************/
/* File:   adtree.hh                                                       */
/* Author: Joachim Schoeberl                                               */
/* Date:   16. Feb. 98                                                     */
/* Redesigned by Wolfram Muehlhuber, May 1998                              */
/* *************************************************************************/

namespace netgen
{

/**
  Alternating Digital Tree
 */

class ADTreeNode
{
public:
  ADTreeNode *left, *right, *father;
  int dim;
  float sep;
  float *data;
  float *boxmin;
  float *boxmax;
  int pi;
  int nchilds;

  ADTreeNode (int adim);
  ~ADTreeNode ();

  friend class ADTree;
};

class ADTreeCriterion
{
public:
  ADTreeCriterion() { }
  virtual int Eval (const ADTreeNode * node) const = 0;
};

class ADTree
{
  int dim;
  ADTreeNode * root;
  float *cmin, *cmax;
  Array<ADTreeNode*> ela;
  const ADTreeCriterion * criterion; 

  Array<ADTreeNode*> stack;
  Array<int> stackdir;
  int stackindex;

public:
  ADTree (int adim, const float * acmin, 
	   const float * acmax);
  ~ADTree ();

  void Insert (const float * p, int pi);
  void SetCriterion (ADTreeCriterion & acriterion);
  void Reset ();
  int Next ();
  void GetMatch (Array<int> & matches);

  void DeleteElement (int pi);

};

class ADTreeNode3
{
public:
  ADTreeNode3 *left, *right, *father;
  float sep;
  float data[3];
  int pi;
  int nchilds;

  ADTreeNode3 ();
  void DeleteChilds ();
  friend class ADTree3;

  static BlockAllocator ball;
  void * operator new(size_t);
  void operator delete (void *);
};

class ADTree3
{
  ADTreeNode3 * root;
  float cmin[3], cmax[3];
  Array<ADTreeNode3*> ela;

public:
  ADTree3 (const float * acmin, 
	   const float * acmax);
  ~ADTree3 ();

  void Insert (const float * p, int pi);
  void GetIntersecting (const float * bmin, const float * bmax,
			Array<int> & pis) const;
  
  void DeleteElement (int pi);

};

class ADTreeNode6
{
public:
  ADTreeNode6 *left, *right, *father;
  float sep;
  float data[6];
  int pi;
  int nchilds;

  ADTreeNode6 ();
  void DeleteChilds ();
  friend class ADTree6;

  static BlockAllocator ball;
  void * operator new(size_t);
  void operator delete (void *);
};

class ADTree6
{
  ADTreeNode6 * root;
  float cmin[6], cmax[6];
  Array<ADTreeNode6*> ela;

public:
  ADTree6 (const float * acmin, 
	   const float * acmax);
  ~ADTree6 ();

  void Insert (const float * p, int pi);
  void GetIntersecting (const float * bmin, const float * bmax,
			Array<int> & pis) const;
  
  void DeleteElement (int pi);

  int Depth () const
  { return DepthRec (root); }
  int Elements () const
  { return ElementsRec (root); }

  int DepthRec (const ADTreeNode6 * node) const;
  int ElementsRec (const ADTreeNode6 * node) const;
};

class Point3dTree 
{
  ADTree3 * tree;

public:
   Point3dTree (const Point<3> & pmin, const Point<3> & pmax);
   ~Point3dTree ();
   void Insert (const Point<3> & p, int pi);
  void DeleteElement (int pi) 
    { tree->DeleteElement(pi); }
   void GetIntersecting (const Point<3> & pmin, const Point<3> & pmax, 
			Array<int> & pis) const;
  const ADTree3 & Tree() const { return *tree; };
};

class Box3dTree
{
  ADTree6 * tree;
  Point<3> boxpmin, boxpmax;
public:
  Box3dTree (const Box<3> & abox);
  Box3dTree (const Point<3> & apmin, const Point<3> & apmax);
  ~Box3dTree ();
  void Insert (const Point<3> & bmin, const Point<3> & bmax, int pi);
  void Insert (const Box<3> & box, int pi)
  {
    Insert (box.PMin(), box.PMax(), pi);
  }
  void DeleteElement (int pi) 
    { tree->DeleteElement(pi); }
  void GetIntersecting (const Point<3> & pmin, const Point<3> & pmax, 
			Array<int> & pis) const;

  const ADTree6 & Tree() const { return *tree; };
};
}