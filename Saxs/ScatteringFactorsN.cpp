/*
 * ScatteringFactorsN.cpp
 *
 *  Created on: Feb 1, 2017
 *      Author: marchi
 */

#include "ScatteringFactorsN.h"

ScatteringFactorsN::ScatteringFactorsN():ScatteringFactors::ScatteringFactors() {
	it1992.clear();
	it1992=	{
			{"DU",  {{1.0   }}},
			{"H",   {{-3.7406}}},
			{"D",   {{6.671  }}},
			{"C",   {{6.6460 }}},
			{"O",   {{5.803  }}},
			{"N",   {{9.36   }}},
			{"S",   {{2.847  }}},
			{"P",   {{5.13   }}},
			{"K",   {{3.67   }}},
			{"Iz",  {{-4.5758}}},
			{"Iy",  {{-0.8352}}},
			{"Ix",  {{2.9054 }}},
			{"Na",  {{3.63   }}},
			{"Cl",  {{9.5770 }}},
			{"Ca",  {{4.70   }}},
			{"Fe",  {{9.45   }}},
			{"Mg",  {{5.375  }}}
	};
}

ScatteringFactorsN::~ScatteringFactorsN() {
	// TODO Auto-generated destructor stub
}

