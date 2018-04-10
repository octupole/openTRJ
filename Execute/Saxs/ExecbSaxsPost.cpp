/*
 * ExecbSaxsPost.cpp
 *
 *  Created on: Mar 6, 2016
 *      Author: marchi
 */

#include "ExecbSaxsPost.h"
ExecbSaxsPost::ExecbSaxsPost(trj::TrjRead & MyIn): ExecbSaxs(MyIn){
	__SetUp(MyIn);

};

void ExecbSaxsPost::__SetUp(trj::TrjRead & MyIn){
	ExecbSaxs::__SetUp(MyIn);
	fin1x=MyIn.gFin1();
	fin2x=MyIn.gFin2();
	fin_contrast=MyIn.gFin_contrast();
	fout_saxsx=MyIn.gFoutsaxs();

}

void ExecbSaxsPost::operator()(MAtoms * atm){
	__Differences();
	__GofR();
}
void ExecbSaxsPost::__GofR(){
	MySaxs->GofR();
}
void ExecbSaxsPost::__Differences(){
	Saxs MySaxs1(Myd,Mycut);
	Saxs MySaxs2(Myd,Mycut);
	MySaxs1.SetMass(MassSolute);
	if(bnoSplineOut){
		MySaxs1.SetSplineout();
		MySaxs2.SetSplineout();
	}
	MySaxs1.SetupQdf();
	MySaxs2.SetupQdf();
	MySaxs1.bReadIq(*fin1x);
	if(fin2x)
		MySaxs2.bReadIq(*fin2x);

	MySaxs1-=MySaxs2;
	MySaxs=new Saxs(MySaxs1);
	bnoSplineOut=true;
	MySaxs->SetSplineout();
};

ExecbSaxsPost::~ExecbSaxsPost() {
	// TODO Auto-generated destructor stub
}
