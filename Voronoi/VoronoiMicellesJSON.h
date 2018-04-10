/*
 * VoronoiMicellesJSON.h
 *
 *  Created on: Dec 17, 2017
 *      Author: marchi
 */

#ifndef VORONOI_VORONOIMICELLESJSON_H_
#define VORONOI_VORONOIMICELLESJSON_H_

#include "VoronoiMicelles.h"
#include "json.hpp"

namespace Voro {
using nlohmann::json;

class VoronoiMicellesJSON: public VoronoiMicelles {
protected:
	void WriteIt(std::ofstream &);
public:
	using VoronoiMicelles::VoronoiMicelles;
	void WriteLastJSON(std::ofstream &);

	virtual ~VoronoiMicellesJSON();
};

} /* namespace Voro */

#endif /* VORONOI_VORONOIMICELLESJSON_H_ */
