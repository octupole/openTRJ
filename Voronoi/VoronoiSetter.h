/*
 * VoronoiPrint.h
 *
 *  Created on: Oct 13, 2011
 *      Author: marchi
 */

#ifndef VORONOIPRINT_H_
#define VORONOIPRINT_H_

class VoronoiSetter {
public:
	static int pGroup;
	static bool bPrintVols;
	static bool bPrintAreas;
	static bool bPrintShell;
	static int maxLevel;
	VoronoiSetter();
	virtual ~VoronoiSetter();
};

#endif /* VORONOIPRINT_H_ */
