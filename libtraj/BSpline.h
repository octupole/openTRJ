/*
 * BSpline.h
 *
 *  Created on: Jun 26, 2015
 *      Author: marchi
 */

#ifndef SRC_BSPLINE_H_
#define SRC_BSPLINE_H_
#include <vector>

using std::vector;
struct spline{
	vector<double> x;
	vector<double> dx;
};
class BSpline {
	static int order;
	spline theta;
	void allocate();
	void Init(const double);
	void OnePass(const double ,const int);
	void Diff();
public:
	BSpline();
	BSpline(int);
	virtual ~BSpline();
	spline & operator()(const double);
	static int Order(){return order;}
	static void SetOrder(int MyOrder){order=MyOrder;}

};

#endif /* SRC_BSPLINE_H_ */
