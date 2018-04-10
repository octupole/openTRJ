/*
 * Timer.h
 *
 *  Created on: Jun 17, 2015
 *      Author: marchi
 */

#ifndef SRC_TIMER_H_
#define SRC_TIMER_H_
#include <ostream>
#include <chrono>

namespace timer {
typedef std::chrono::high_resolution_clock high_resolution_clock;
typedef std::chrono::milliseconds milliseconds;
typedef std::chrono::seconds seconds;

class Timer {
private:
	high_resolution_clock::time_point _start;
public:
	explicit Timer()
	{
		Reset();
	}
	void Reset()
	{
		_start = high_resolution_clock::now();
	}

	std::chrono::duration<double> Elapsed() const
	{
		return std::chrono::duration_cast<std::chrono::duration<double> >(high_resolution_clock::now() - _start);
	}
	double operator/(double y){return this->Elapsed().count()/y;}
	double operator*(double y){return this->Elapsed().count()*y;}
	template <typename T, typename Traits>
	friend std::basic_ostream<T, Traits>& operator<<(std::basic_ostream<T, Traits>& out, const Timer& timer)
	{
		return out << timer.Elapsed().count();
	}

	double val()
	{
		return this->Elapsed().count();
	}
	virtual ~Timer();
};

} /* namespace timer */

#endif /* SRC_TIMER_H_ */
