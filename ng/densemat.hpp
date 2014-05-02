#pragma once
/**************************************************************************/
/* File:   densemat.hpp                                                    */
/* Author: Joachim Schoeberl                                              */
/* Date:   01. Oct. 94                                                    */
/**************************************************************************/

/** 
    Data type dense matrix
*/

class DenseMatrix
{
protected:
  int height;
  int width;
  double * data;

public:
  ///
   DenseMatrix ();
  ///
   DenseMatrix (int h, int w = 0);
  ///
   DenseMatrix (const DenseMatrix & m2);
  ///
   ~DenseMatrix ();

  ///
   void SetSize (int h, int w = 0);

  int Height() const { return height; }
  int Width() const {return width; }

  double & operator() (int i, int j) { return data[i*width+j]; }
  double operator() (int i, int j) const { return data[i*width+j]; }
  double & operator() (int i) { return data[i]; }
  double operator() (int i) const { return data[i]; }

  ///
   DenseMatrix & operator= (const DenseMatrix & m2);
  ///
   DenseMatrix & operator+= (const DenseMatrix & m2);
  ///
   DenseMatrix & operator-= (const DenseMatrix & m2);

  ///
   DenseMatrix & operator= (double v);
  ///
   DenseMatrix & operator*= (double v);

  ///
   void Mult (const FlatVector & v, FlatVector & prod) const
  {
    const double *mp = data;
    double *dp = &prod(0);
    for (int i = 0; i < height; i++)
    {
        double sum = 0;
        const double *sp = &v(0);

        for (int j = 0; j < width; j++)
        {
            sum += (*mp) * (*sp);
            mp++;
            sp++;
        }

        *dp = sum;
        dp++;
    }
  }

  ///
   void MultTrans (const Vector & v, Vector & prod) const;
  ///
   void Residuum (const Vector & x, const Vector & b, Vector & res) const;
  ///
   double Det () const;
  ///
  friend DenseMatrix operator* (const DenseMatrix & m1, const DenseMatrix & m2);
  ///
  friend DenseMatrix operator+ (const DenseMatrix & m1, const DenseMatrix & m2);

  /// 
  friend void Transpose (const DenseMatrix & m1, DenseMatrix & m2);
  ///
  friend void Mult (const DenseMatrix & m1, const DenseMatrix & m2, DenseMatrix & m3);
  ///
  friend void CalcAAt (const DenseMatrix & a, DenseMatrix & m2);
  ///
  friend void CalcABt (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2);
  ///
  friend void CalcAtB (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2);
  ///
   void Solve (const Vector & b, Vector & x) const;
  ///
  void SolveDestroy (const Vector & b, Vector & x);
  ///
  const double & Get(int i, int j) const { return data[(i-1)*width+j-1]; }
  ///
  const double & Get(int i) const { return data[i-1]; }
  ///
  void Set(int i, int j, double v) { data[(i-1)*width+j-1] = v; }
  ///
  double & Elem(int i, int j) { return data[(i-1)*width+j-1]; }
  ///
  const double & ConstElem(int i, int j) const { return data[(i-1)*width+j-1]; }
};

extern ostream & operator<< (ostream & ost, const DenseMatrix & m);

template <int WIDTH>
class MatrixFixWidth
{
protected:
  int height;
  double * data;

public:
  ///
  MatrixFixWidth () 
  { height = 0; data = 0; }
  ///
  MatrixFixWidth (int h)
  { height = h; data = new double[WIDTH*height]; }
  ///
  ~MatrixFixWidth ()
  { delete [] data; }

  void SetSize (int h)
  {
    if (h != height)
      {
	delete data;
	height = h;
	data = new double[WIDTH*height]; 
      }
  }

  ///
  int Height() const { return height; }

  ///
  int Width() const { return WIDTH; }

  ///
  MatrixFixWidth & operator= (double v)
  {
    for (int i = 0; i < height*WIDTH; i++)
      data[i] = v; 
    return *this;
  }

  ///
  void Mult (const FlatVector & v, FlatVector & prod) const
  {
    double sum;
    const double * mp, * sp;
    double * dp;

    mp = data;
    dp = &prod[0];
    for (int i = 0; i < height; i++)
      {
	sum = 0;
	sp = &v[0];
	
	for (int j = 0; j < WIDTH; j++)
	  {
	    sum += *mp * *sp;
	    mp++;
	    sp++;
	  }
	    
	*dp = sum;
	dp++;
      }
  }

  double & operator() (int i, int j)
  { return data[i*WIDTH+j]; }

  const double & operator() (int i, int j) const
  { return data[i*WIDTH+j]; }


  MatrixFixWidth & operator*= (double v)
  {
    if (data)
      for (int i = 0; i < height*WIDTH; i++)
        data[i] *= v;
    return *this;
  }

  const double & Get(int i, int j) const { return data[(i-1)*WIDTH+j-1]; }
  ///
  const double & Get(int i) const { return data[i-1]; }
  ///
  void Set(int i, int j, double v) { data[(i-1)*WIDTH+j-1] = v; }
  ///
  double & Elem(int i, int j) { return data[(i-1)*WIDTH+j-1]; }
  ///
  const double & ConstElem(int i, int j) const { return data[(i-1)*WIDTH+j-1]; }
};

template <int WIDTH>
extern ostream & operator<< (ostream & ost, const MatrixFixWidth<WIDTH> & m)
{
  for (int i = 0; i < m.Height(); i++)
    {
      for (int j = 0; j < m.Width(); j++)
	ost << m.Get(i+1,j+1) << " ";
      ost << endl;
    }
  return ost;
};

extern  void CalcAtA (const DenseMatrix & a, DenseMatrix & m2);
extern  void CalcInverse (const DenseMatrix & m1, DenseMatrix & m2);