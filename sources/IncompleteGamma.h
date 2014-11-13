
/********************

PhyloBayes MPI. Copyright 2010-2013 Nicolas Lartillot, Nicolas Rodrigue, Daniel Stubbs, Jacques Richer.

PhyloBayes is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
PhyloBayes is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

See the GNU General Public License for more details. You should have received a copy of the GNU General Public License
along with PhyloBayes. If not, see <http://www.gnu.org/licenses/>.

**********************/


# include <cstdlib>
# include <cstdio>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

double PointNormal(double prob);

//yan modify 12.Dec.2004 (from Mr.Bayes)

//------------------------------------------------------------------------------
//
//  Returns the incomplete gamma ratio I(x,alpha) where x is the upper
//  limit of the integration and alpha is the shape parameter.  Returns (-1)
//  if in error.
// LnGamma_alpha = ln(Gamma(alpha)), is almost redundant.
//  (1) series expansion     if (alpha>x || x<=1)
//  (2) continued fraction   otherwise
//
//  RATNEST FORTRAN by
//  Bhattacharjee, G. P.  1970.  The incomplete gamma integral.  Applied
//     Statistics, 19:285-287 (AS32)
//
//------------------------------------------------------------------------------
double IncompleteGamma (double x, double alpha, double LnGamma_alpha)   {
	int 		i;
	double		p = alpha;
        double          g = LnGamma_alpha;
        double		accurate = 1e-8;
        double          overflow = 1e30;
        double		factor;
        double          gin = 0.0;
        double          rn = 0.0;
        double          a = 0.0;
        double          b = 0.0;
        double          an = 0.0;
        double          dif = 0.0;
        double          term = 0.0;
        double          pn[6];

	if (x == 0.0)
		return (0.0);
	if (x < 0 || p <= 0)
		return (-1.0);

	factor = exp(p*log(x)-x-g);
	if (x>1 && x>=p)
		goto l30;
	gin = 1.0;
	term = 1.0;
	rn = p;
	l20:
		rn++;
		term *= x/rn;
		gin += term;
		if (term > accurate)
			goto l20;
		gin *= factor/p;
		goto l50;
	l30:
		a = 1.0-p;
		b = a+x+1.0;
		term = 0.0;
		pn[0] = 1.0;
		pn[1] = x;
		pn[2] = x+1;
		pn[3] = x*b;
		gin = pn[2]/pn[3];
	l32:
		a++;
		b += 2.0;
		term++;
		an = a*term;
		for (i=0; i<2; i++)
			pn[i+4] = b*pn[i+2]-an*pn[i];
		if (pn[5] == 0)
			goto l35;
		rn = pn[4]/pn[5];
		dif = fabs(gin-rn);
		if (dif>accurate)
			goto l34;
		if (dif<=accurate*rn)
			goto l42;
	l34:
		gin = rn;
	l35:
		for (i=0; i<4; i++)
			pn[i] = pn[i+2];
		if (fabs(pn[4]) < overflow)
			goto l32;
		for (i=0; i<4; i++)
			pn[i] /= overflow;
		goto l32;
	l42:
		gin = 1.0-factor*gin;
	l50:
		return (gin);

}
//yan 20. Dec 2004 (from Mr.Bayes)
//------------------------------------------------------------------------------
//
//  Returns ln(gamma(alpha)) for alpha > 0, accurate to 10 decimal places.
// Stirling's formula is used for the central polynomial part of the procedure.
//
//  Pike, M. C. and I. D. Hill.  1966.  Algorithm 291: Logarithm of the gamma
//    function.  Communications of the Association for Computing
//    Machinery, 9:684.
//
//-----------------------------------------------------------------------------
double LnGamma (double alpha)   {
	double	x = alpha;
        double  f = 0.0;
        double  z;

	if (x < 7)
		{
		f = 1.0;
		z = x-1.0;
		while (++z < 7.0)
			f *= z;
		x = z;
		f = -log(f);
		}
	z = 1.0/(x*x);
	return  f + (x-0.5)*log(x) - x + 0.918938533204673 +
		(((-0.000595238095238*z+0.000793650793651)*z-0.002777777777778)*z +
		0.083333333333333)/x;
}

//http://www.csit.fsu.edu/~burkardt/cpp_src/dcdflib/dcdflib.html
//---------------------------------------------------------------------------
//              double PointChi2 (double prob, double v)
//---------------------------------------------------------------------------
//------------------------------------------------------------------------------
//
//  Returns z so That Prob{x<z} = prob where x is Chi2 distributed with df=v.
//  Returns -1 if in error.   0.000002 < prob < 0.999998.
//
// RATNEST FORTRAN by
//  Best, D. J. and D. E. Roberts.  1975.  The percentage points of the
//     Chi2 distribution.  Applied Statistics 24:385-388.  (AS91)
//
//  Converted into C by Ziheng Yang, Oct. 1993.
//
//-----------------------------------------------------------------------------
double PointChi2 (double prob, double v)        {
	double 		e = 0.5e-6, aa = 0.6931471805, p = prob, g,
					xx, c, ch, a = 0.0, q = 0.0, p1 = 0.0, p2 = 0.0, t = 0.0,
					x = 0.0, b = 0.0, s1, s2, s3, s4, s5, s6;

	if (p < 0.000002 || p > 0.999998 || v <= 0.0)
		return (-1.0);
	g = LnGamma (v/2.0);
	xx = v/2.0;
	c = xx - 1.0;
	if (v >= -1.24*log(p))
		goto l1;
	ch = pow((p*xx*exp(g+xx*aa)), 1.0/xx);
	if (ch-e<0)
		return (ch);
	goto l4;
	l1:
		if (v > 0.32)
			goto l3;
		ch = 0.4;
		a = log(1.0-p);
	l2:
		q = ch;
		p1 = 1.0+ch*(4.67+ch);
		p2 = ch*(6.73+ch*(6.66+ch));
		t = -0.5+(4.67+2.0*ch)/p1 - (6.73+ch*(13.32+3.0*ch))/p2;
		ch -= (1.0-exp(a+g+0.5*ch+c*aa)*p2/p1)/t;
		if (fabs(q/ch-1.0)-0.01 <= 0.0)
			goto l4;
		else
			goto l2;
	l3:
		x = PointNormal (p);
		p1 = 0.222222/v;
		ch = v*pow((x*sqrt(p1)+1.0-p1), 3.0);
		if (ch > 2.2*v+6.0)
			ch = -2.0*(log(1.0-p)-c*log(0.5*ch)+g);
	l4:
		q = ch;
		p1 = 0.5*ch;
		if ((t = IncompleteGamma (p1, xx, g)) < 0.0)
			{
			printf ("\nerr IncompleteGamma");
			return (-1.0);
			}
		p2 = p-t;
		t = p2*exp(xx*aa+g+p1-c*log(ch));
		b = t/ch;
		a = 0.5*t-b*c;
		s1 = (210.0+a*(140.0+a*(105.0+a*(84.0+a*(70.0+60.0*a))))) / 420.0;
		s2 = (420.0+a*(735.0+a*(966.0+a*(1141.0+1278.0*a))))/2520.0;
		s3 = (210.0+a*(462.0+a*(707.0+932.0*a)))/2520.0;
		s4 = (252.0+a*(672.0+1182.0*a)+c*(294.0+a*(889.0+1740.0*a)))/5040.0;
		s5 = (84.0+264.0*a+c*(175.0+606.0*a))/2520.0;
		s6 = (120.0+c*(346.0+127.0*c))/5040.0;
		ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
		if (fabs(q/ch-1.0) > e)
			goto l4;
		return (ch);

}

double  PointNormal (double prob)

{
	double 		a0 = -0.322232431088, a1 = -1.0, a2 = -0.342242088547, a3 = -0.0204231210245,
 					a4 = -0.453642210148e-4, b0 = 0.0993484626060, b1 = 0.588581570495,
 					b2 = 0.531103462366, b3 = 0.103537752850, b4 = 0.0038560700634,
 					y, z = 0, p = prob, p1;

	p1 = (p<0.5 ? p : 1-p);
	if (p1<1e-20)
	   return (-9999);
	y = sqrt (log(1/(p1*p1)));
	z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
	return (p<0.5 ? -z : z);

}

#define PointGamma(prob,alpha,beta) 		PointChi2(prob,2.0*(alpha))/(2.0*(beta)) //yan23.dec2004


