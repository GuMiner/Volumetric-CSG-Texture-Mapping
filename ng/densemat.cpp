#include "stdafx.h" 
#include "mystdlib.h"

#include "linalg.hpp"

namespace netgen
{
  DenseMatrix :: DenseMatrix () 
  {
    data = NULL;
    height = 0;
    width = 0;
  }

  DenseMatrix :: DenseMatrix (int h, int w)
  {
    if (!w) w = h;
    width = w;
    height = h;
    if (h*w)
      data = new double[h*w];
    else 
      data = 0;

    for (int i = 0 ; i < (h * w); i++)
      data[i] = 0;
  }

  DenseMatrix :: DenseMatrix (const DenseMatrix & m2)
  {
    data = NULL; height = width = 0;
    SetSize (m2.Height(), m2.Width());
    memcpy (data, m2.data, sizeof(double) * Height() * Width());
  }

  DenseMatrix :: ~DenseMatrix ()
  {
    delete [] data;
  }
  
  void DenseMatrix :: SetSize (int h, int w)
  {
    if (!w) w = h;
    if (height == h && width == w)
      return;
          
    height = h;
    width = w;
    
    delete[] data;
    
    if (h*w)  
      data = new double[h*w];
    else
      data = NULL;
  }

  DenseMatrix & DenseMatrix :: operator= (const DenseMatrix & m2)
  {
    SetSize (m2.Height(), m2.Width());
    
    if (data) memcpy (data, m2.data, sizeof(double) * m2.Height() * m2.Width());
    return *this;
  }

  DenseMatrix & DenseMatrix :: operator+= (const DenseMatrix & m2)
  {
    int i;
    double * p, * q;
    
    if (Height() != m2.Height() || Width() != m2.Width())
      {
        cout << "DenseMatrix::Operator+=: Sizes don't fit" << endl;
        return *this;
      }
    
    if (data)
      {
	p = data;
	q = m2.data;
	for (i = Width() * Height(); i > 0; i--)
          {
            *p += *q;
            p++;
            q++;
          }
      }
    else
      cout << "DenseMatrix::Operator+=: Matrix not allocated" << endl;

    return *this;
  }

  DenseMatrix & DenseMatrix :: operator-= (const DenseMatrix & m2)
  {
    int i;
    double * p, * q;

    if (Height() != m2.Height() || Width() != m2.Width())
      {
        cout << "DenseMatrix::Operator-=: Sizes don't fit" << endl;
        return *this;
      }

    if (data)
      {
        p = data;
        q = m2.data;
        for (i = Width() * Height(); i > 0; i--)
          {
            *p -= *q;
            p++;
            q++;
          }
      }
    else
      cout << "DenseMatrix::Operator-=: Matrix not allocated" << endl;

    return *this;
  }

  DenseMatrix & DenseMatrix :: operator= (double v)
  {
    double * p = data;

    if (data)
      for (int i = width*height; i > 0; i--, p++)
        *p = v;

    return *this;
  }

  DenseMatrix & DenseMatrix :: operator*= (double v)
  {
    double * p = data;

    if (data)
      for (int i = width*height; i > 0; i--, p++)
        *p *= v;

    return *this;
  }

  double DenseMatrix :: Det () const
  {
    if (width != height)
      {
        cout << "DenseMatrix :: Det: width != height" << endl;
        return 0;
      }

    switch (width)
      {
      case 1: return data[0];
      case 2: return data[0] * data[3] - data[1] * data[2];
      case 3: return data[0] * data[4] * data[8]
          + data[1] * data[5] * data[6]
          + data[2] * data[3] * data[7] 
          - data[0] * data[5] * data[7]
          - data[1] * data[3] * data[8]
          - data[2] * data[4] * data[6];
      default:
        {
          cout << "Matrix :: Det:  general size not implemented (size=" << width << ")" << endl;
          return 0;
        }
      }
  }

  void CalcInverse (const DenseMatrix & m1, DenseMatrix & m2)
  {
    double det;

    if (m1.Width() != m1.Height())
      {
        cout << "CalcInverse: matrix not symmetric" << endl;
        return;
      }
    if (m1.Width() != m2.Width() || m1.Height() != m2.Height())
      {
        cout << "CalcInverse: dim(m2) != dim(m1)" << endl;
        return;
      }

    if (m1.Width() <= 3)
      {
        det = m1.Det();
        if (det == 0)
          {
            cout << "CalcInverse: Matrix singular" << endl;
            return;
          }

        det = 1.0f / det;
        switch (m1.Width())
          {
          case 1:
            {
              m2(0,0) = det;
              return;
            }
          case 2:
            {
              m2(0,0) = det * m1(3);
              m2(1,1) = det * m1(0);  
              m2(0,1) = -det * m1(1);
              m2(1,0) = - det * m1(2);
              return;
            }
          case 3:
            {
              m2(0, 0) =  det * (m1(4) * m1(8) - m1(5) * m1(7));
              m2(1, 0) = -det * (m1(3) * m1(8) - m1(5) * m1(6));
              m2(2, 0) =  det * (m1(3) * m1(7) - m1(4) * m1(6));

              m2(0, 1) = -det * (m1(1) * m1(8) - m1(2) * m1(7));
              m2(1, 1) =  det * (m1(0) * m1(8) - m1(2) * m1(6));
              m2(2, 1) = -det * (m1(0) * m1(7) - m1(1) * m1(6));

              m2(0, 2) =  det * (m1(1) * m1(5) - m1(2) * m1(4));
              m2(1, 2) = -det * (m1(0) * m1(5) - m1(2) * m1(3));
              m2(2, 2) =  det * (m1(0) * m1(4) - m1(1) * m1(3));
              return;
            }
          }
      }
    
    else
      {
        int i, j, k, n;
        n = m1.Height();
  
        // Gauss - Jordan - algorithm
        int r, hi;
        double max, hr;
      

        Array<int> p(n);   // pivot-permutation
        Vector hv(n);
    
      
        m2 = m1;
      
        // Algorithm of Stoer, Einf. i. d. Num. Math, S 145
        for (j = 1; j <= n; j++)
          p.Set(j, j);
      
        for (j = 1; j <= n; j++)
          {
            // pivot search
            max = fabs(m2.Get(j, j));
            r = j;
	  
            for (i = j+1; i <= n ;i++)
              if (fabs (m2.Get(i, j)) > max)
                {
                  r = i;
                  max = fabs (m2.Get(i, j));
                }
	  
            if (max < 1e-20)
              {
                cerr << "Inverse matrix: matrix singular" << endl;
                return;
              }
	  
            r = j;
	  
            // exchange rows
            if (r > j)
              {
                for (k = 1; k <= n; k++)
                  {
                    hr = m2.Get(j, k);
                    m2.Elem(j, k) = m2.Get(r, k);
                    m2.Elem(r, k) = hr;
                  }
                hi = p.Get(j);
                p.Elem(j) = p.Get(r);
                p.Elem(r) = hi;
              }
	  
	  
            // transformation
            hr = 1 / m2.Get(j, j);
            for (i = 1; i <= n; i++)
              m2.Elem(i, j) *= hr;
            m2.Elem(j, j) = hr;
	  
            for (k = 1; k <= n; k++)
              if (k != j)
                {
                  for (i = 1; i <= n; i++)
                    if (i != j)
                      m2.Elem(i, k) -= m2.Elem(i, j) * m2.Elem(j, k);
                  m2.Elem(j, k) *= -hr;
                }
          }
      
        // col exchange
        for (i = 1; i <= n; i++)
          {
            for (k = 1; k <= n; k++)
              hv(p.Get(k)-1) = m2.Get(i, k);
            for (k = 1; k <= n; k++)
              m2.Elem(i, k) = hv(k-1);
          }
      }
  }

  void CalcAAt (const DenseMatrix & a, DenseMatrix & m2)
  {
    int n1 = a.Height();
    int n2 = a.Width();
    int i, j, k;
    double sum;
    const double *p, *q, *p0;

    if (m2.Height() != n1 || m2.Width() != n1)
      {
        cout << "CalcAAt: sizes don't fit" << endl;
        return;
      }

    for (i = 1; i <= n1; i++)
      {
        sum = 0;
        p = &a.ConstElem(i, 1);
        for (k = 1; k <= n2; k++)
          {
            sum += *p * *p;
            p++;
          }
        m2.Set(i, i, sum);

        p0 = &a.ConstElem(i, 1);
        q = a.data;
        for (j = 1; j < i; j++)
          {
            sum = 0;
            p = p0;

            for (k = 1; k <= n2; k++)
              {
                sum += *p * *q;
                p++;
                q++;
              }
            m2.Set(i, j, sum);
            m2.Set(j, i, sum);
          }
      }
  }

  void CalcAtA (const DenseMatrix & a, DenseMatrix & m2)
  {
    int n1 = a.Height();
    int n2 = a.Width();
    int i, j, k;
    double sum;

    if (m2.Height() != n2 || m2.Width() != n2)
      {
        cout << "CalcAtA: sizes don't fit" << endl;
        return;
      }

    for (i = 1; i <= n2; i++)
      for (j = 1; j <= n2; j++)
        {
          sum = 0;
          for (k = 1; k <= n1; k++)
            sum += a.Get(k, i) * a.Get(k, j);
          m2.Elem(i, j) = sum;
        }
  }

  void CalcABt (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2)
  {
    int n1 = a.Height();
    int n2 = a.Width();
    int n3 = b.Height();
    int i, j, k;
    double sum;

    if (m2.Height() != n1 || m2.Width() != n3 || b.Width() != n2)
      {
        cout << "CalcABt: sizes don't fit" << endl;
        return;
      }

    double * pm2 = &m2.Elem(1, 1);
    const double * pa1 = &a.Get(1, 1);

    for (i = 1; i <= n1; i++)
      {
        const double * pb = &b.Get(1, 1);
        for (j = 1; j <= n3; j++)
          {
            sum = 0;
            const double * pa = pa1;
	  
            for (k = 1; k <= n2; k++)
              {
                sum += *pa * *pb;
                pa++; pb++;
              }
	  
            *pm2 = sum;
            pm2++;
          }
        pa1 += n2;
      }
  }

  void CalcAtB (const DenseMatrix & a, const DenseMatrix & b, DenseMatrix & m2)
  {
    int n1 = a.Height();
    int n2 = a.Width();
    int n3 = b.Width();
    int i, j, k;

    if (m2.Height() != n2 || m2.Width() != n3 || b.Height() != n1)
      {
        cout << "CalcAtB: sizes don't fit" << endl;
        return;
      }

    for (i = 1; i <= n2 * n3; i++)
      m2.data[i-1] = 0;

    for (i = 1; i <= n1; i++)
      for (j = 1; j <= n2; j++)
        {
          const double va = a.Get(i, j);
          double * pm2 = &m2.Elem(j, 1);
          const double * pb = &b.Get(i, 1);

          for (k = 1; k <= n3; ++k, ++pm2, ++pb)
            *pm2 += va * *pb;
        }
  }

  DenseMatrix operator* (const DenseMatrix & m1, const DenseMatrix & m2)
  {
    DenseMatrix temp (m1.Height(), m2.Width());

    if (m1.Width() != m2.Height())
      {
        cout << "DenseMatrix :: operator*: Matrix Size does not fit" << endl;
      }
    else if (temp.Height() != m1.Height())
      {
        cout << "DenseMatrix :: operator*: temp not allocated" << endl;
      }
    else
      {
        Mult (m1, m2, temp);
      }
    return temp;
  }

  void Mult (const DenseMatrix & m1, const DenseMatrix & m2, DenseMatrix & m3)
  {
    double sum;
    double *p1, *p1s, *p1sn, *p1snn, *p2, *p2s, *p2sn, *p3;

    if (m1.Width() != m2.Height() || m1.Height() != m3.Height() ||
        m2.Width() != m3.Width() )
      {
        cout << "DenseMatrix :: Mult: Matrix Size does not fit" << endl;
        cout << "m1: " << m1.Height() << " x " << m1.Width() << endl;
        cout << "m2: " << m2.Height() << " x " << m2.Width() << endl;
        cout << "m3: " << m3.Height() << " x " << m3.Width() << endl;
        return;
      }
    else
      {
        //      int i, j, k;
        int n1 = m1.Height();
        int n2 = m2.Width();
        int n3 = m1.Width();

        p3 = m3.data;
        p1s = m1.data;
        p2sn = m2.data + n2;
        p1snn = p1s + n1 * n3;

        while (p1s != p1snn)
          {
            p1sn = p1s + n3;
            p2s = m2.data;
	  
            while (p2s != p2sn)
              {
                sum = 0;
                p1 = p1s;
                p2 = p2s;
                p2s++;

                while (p1 != p1sn)
                  {
                    sum += *p1 * *p2;
                    p1++;
                    p2 += n2;
                  }
                *p3++ = sum;
              }
            p1s = p1sn;
          }
      }
  }  

  DenseMatrix operator+ (const DenseMatrix & m1, const DenseMatrix & m2)
  {
    DenseMatrix temp (m1.Height(), m1.Width());
    int i, j;

    if (m1.Width() != m2.Width() || m1.Height() != m2.Height())
      {
        cout << "BaseMatrix :: operator+: Matrix Size does not fit" << endl;
      }
    else if (temp.Height() != m1.Height())
      {
        cout << "BaseMatrix :: operator+: temp not allocated" << endl;
      }
    else
      {
        for (i = 1; i <= m1.Height(); i++)
          for (j = 1; j <= m1.Width(); j++)
            {
              temp.Set(i, j, m1.Get(i, j) + m2.Get(i, j));
            }
      }
    return temp;
  }

  void Transpose (const DenseMatrix & m1, DenseMatrix & m2)
  {
    int w = m1.Width();
    int h = m1.Height();
    int i, j;

    m2.SetSize (w, h);

    double * pm2 = &m2.Elem(1, 1);
    for (j = 1; j <= w; j++)
      {
        const double * pm1 = &m1.Get(1, j);
        for (i = 1; i <= h; i++)
          {
            *pm2 = *pm1;
            pm2 ++;
            pm1 += w;
          }
      }
  }

  void DenseMatrix :: MultTrans (const Vector & v, Vector & prod) const
  {
      int i, j;
      int w = Width(), h = Height();
      if (prod.Size() != w)
	prod.SetSize (w);

      const double * pmat = &Get(1, 1);
      const double * pv = &v(0);

      prod = 0;

      for (i = 1; i <= h; i++)
	{
	  double val = *pv;
	  ++pv;

	  double * pprod = &prod(0);

	  for (j = w-1; j >= 0; --j, ++pmat, ++pprod)
	    {
	      *pprod += val * *pmat;
	    }
	}
  }

  void DenseMatrix :: Residuum (const Vector & x, const Vector & b,
                                Vector & res) const
  {
    double sum;

    res.SetSize (Height());

    if (Width() != x.Size() || Height() != b.Size())
      {
        cout << "\nMatrix and Vector don't fit" << endl;
      }
    else if (Height() != res.Size())
      {
        cout << "Base_Matrix::operator*(Vector): prod vector not ok" << endl;
      }
    else
      {
        int h = Height(); 
        int w = Width();
        const double * mp = &Get(1, 1);

        for (int i = 1; i <= h; i++)
          {
            sum = b(i-1);
            const double * xp = &x(0);

            for (int j = 1; j <= w; ++j, ++mp, ++xp)
              sum -= *mp * *xp;
	  
            res(i-1) = sum;
          }
      }
  }

  void DenseMatrix :: Solve (const Vector & v, Vector & sol) const
  {
    DenseMatrix temp (*this);
    temp.SolveDestroy (v, sol);
  }

  void DenseMatrix :: SolveDestroy (const Vector & v, Vector & sol)
  {
    double q;

    if (Width() != Height())
      {
        cout << "SolveDestroy: Matrix not square";
        return;
      }
    if (Width() != v.Size())
      {
        cout << "SolveDestroy: Matrix and Vector don't fit";
        return;
      }

    sol = v;
    if (Height() != sol.Size())
      {
        cout << "SolveDestroy: Solution Vector not ok";
        return;
      }

        int n = Height();
        for (int i = 1; i <= n; i++)
          {
            for (int j = i+1; j <= n; j++)
              {
                q = Get(j,i) / Get(i,i);
                if (q)
                  {
                    const double * pik = &Get(i, i+1);
                    double * pjk = &Elem(j, i+1);

                    for (int k = i+1; k <= n; ++k, ++pik, ++pjk)
                      *pjk -= q * *pik;

                    //  for (k = i+1; k <= Height(); k++)
                    //	Elem(j, k) -= q * Get(i,k);

                    sol(j-1) -= q * sol(i-1);
                  }
              }
          }
      
        for (int i = n; i >= 1; i--)
          {
            q = sol(i-1);
            for (int j = i+1; j <= n; j++)
	      q -= Get(i,j) * sol(j-1);

            sol(i-1) = q / Get(i,i);
          }
  }

  ostream & operator<< (ostream & ost, const DenseMatrix & m)
  {
    for (int i = 0; i < m.Height(); i++)
      {
        for (int j = 0; j < m.Width(); j++)
          ost << m.Get(i+1,j+1) << " ";
        ost << endl;
      }
    return ost;
  }
}