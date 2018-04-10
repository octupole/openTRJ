/*
 * LagrangeInterpolation.cpp
 *
 *  Created on: Jul 7, 2015
 *      Author: marchi
 */

#include "LagrangeInterpolation.h"
auto x2=[](double x){return x*x;};
auto x3=[](double x){return x*x*x;};
auto x4=[](double x){return x*x*x*x;};
auto x5=[](double x){return x*x*x*x*x;};
auto x6=[](double x){return x*x*x*x*x;};
static bool ok(double x){return x-1.0 > 1.0e-2 || x < 0.0?false:true;}

map<int,map<int,function<double(double)> > >  LagrangeInterpolation::Interpolant=
{{
		{2,{{
				{-1,[](double x){return ok(x)?
						0.5*(x2(x)-x)
						:0.0;}},
				{ 0,[](double x){return ok(x)?
						(-x2(x)+1)
						:0.0;}},
				{ 1,[](double x){return ok(x)?
						0.5*(x2(x)+x)
						:0.0;}}
		}}},
		{4,{{
				{-2,[](double x){return ok(x)?
						(x4(x)-2.0*x3(x)-x2(x)+2.0*x)/24.0
						:0.0;}},
				{-1,[](double x){return ok(x)?
						(-4*x4(x)+4*x3(x)+16*x2(x)-16*x)/24.0
						:0.0;}},
				{ 0,[](double x){return ok(x)?
						(6*x4(x)-30*x2(x)+24)/24.0
						:0.0;}},
				{ 1,[](double x){return ok(x)?
						(-4*x4(x)-4*x3(x)+16*x2(x)+16*x)/24.0
						:0.0;}},
				{ 2,[](double x){return ok(x)?
						(x4(x)+2.0*x3(x)-x2(x)-2.0*x)/24.0
						:0.0;}}

		}}},
		{6,{{
				{-3,[](double x){return ok(x)?
						(x6(x)-3*x5(x)-5*x4(x)+15*x3(x)+4*x2(x)-12*x)/720.0
						:0.0;}},
				{-2,[](double x){return ok(x)?
						(-6*x6(x)+12*x5(x)+60*x4(x)-120*x3(x)-54*x2(x)+108*x)/720.0
						:0.0;}},
				{-1,[](double x){return ok(x)?
						(15*x6(x)-15*x5(x)-195*x4(x)+195*x3(x)+540*x2(x)-540*x)/720.0
						:0.0;}},
				{ 0,[](double x){return ok(x)?
						(-20*x6(x)+280*x4(x)-980*x2(x)+720)/720.0
						:0.0;}},
				{ 1,[](double x){return ok(x)?
						(15*x6(x)+15*x5(x)-195*x4(x)-195*x3(x)+540*x2(x)+540*x)/720.0
						:0.0;}},
				{ 2,[](double x){return ok(x)?
						(-6*x6(x)-12*x5(x)+60*x4(x)+120*x3(x)-54*x2(x)-108*x)/720.0
						:0.0;}},
				{ 3,[](double x){return ok(x)?
						(x6(x)+3*x5(x)-5*x4(x)-15*x3(x)+4*x2(x)+12*x)/720.0
						:0.0;}}

		}}}
}};


LagrangeInterpolation::LagrangeInterpolation() {
	// TODO Auto-generated constructor stub
	MyPoly=Interpolant[order];
}
LagrangeInterpolation::LagrangeInterpolation(int MyOrder) {
	try{
		stringstream ss;ss<<MaxOrder;
		if(MyOrder%2) throw string("Cannot run with an odd order!");
		if(MyOrder>MaxOrder) throw string("Cannot run with an order larger than ")+ss.str()+string(". \n");
	} catch(const string & s){
		cout << s << endl;
		exit(0);
	}
	order=MyOrder;
	MyPoly=Interpolant[order];
}

LagrangeInterpolation::~LagrangeInterpolation() {
	// TODO Auto-generated destructor stub
}

