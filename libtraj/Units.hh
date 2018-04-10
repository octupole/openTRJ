/*
 * Units.hh
 *
 *  Created on: Jan 5, 2016
 *      Author: marchi
 */

#ifndef TRJLIB_UNITS_HH_
#define TRJLIB_UNITS_HH_

#include <iostream>
#include <string>


namespace Units {
template<int M, int K, int S, int Q, int T>
struct Unit {
	enum {
		m = M, kg = K, s = S, q = Q, t = T
	};
};
using M=Unit<1,0,0,0,0>;
using Kg=Unit<0,1,0,0,0>;
using S=Unit<0,0,1,0,0>;
using C=Unit<0,0,0,1,0>;
using T=Unit<0,0,0,0,1>;

using A=Unit<0,0,-1,1,0>;
using MpS=Unit<1,0,-1,0,0>;
using MpS2=Unit<1,0,-2,0,0>;
using Newt=Unit<1,1,-2,0,0>;
using Joules=Unit<2,1,-2,0,0>;
using Volts=Unit<2,1,-2,-1,0>;
using Farads=Unit<-2,-1,2,2,0>;
using Adim=Unit<0,0,0,0,0>;

template<typename U1, typename U2>
struct Uplus {
	using type=Unit<U1::m+U2::m,U1::kg+U2::kg,U1::s+U2::s,U1::q+U2::q,U1::t+U2::t>;
};
template<typename U1, typename U2>
struct Uminus {
	using type=Unit<U1::m-U2::m,U1::kg-U2::kg,U1::s-U2::s,U1::q-U2::q,U1::t-U2::t>;
};

template<typename U1, typename U2>
using Unit_plus=typename Uplus<U1,U2>::type;

template<typename U1, typename U2>
using Unit_minus=typename Uminus<U1,U2>::type;

template<typename U>
struct Quantity {
	long double val;
	explicit constexpr Quantity(long double d) :
			val { d } {
	}
};
template<typename U>
Quantity<U> operator+(Quantity<U> x, Quantity<U> y) {
	return Quantity<U> { x.val + y.val };
}

template<typename U>
Quantity<U> operator-(Quantity<U> x, Quantity<U> y) {
	return Quantity<U> { x.val - y.val };
}

template<typename U1, typename U2>
Quantity<Unit_plus<U1, U2> > operator *(Quantity<U1> x, Quantity<U2> y) {
	return Quantity<Unit_plus<U1, U2> > { x.val * y.val };
}

template<typename U1, typename U2>
Quantity<Unit_minus<U1, U2> > operator /(Quantity<U1> x, Quantity<U2> y) {
	return Quantity<Unit_minus<U1, U2> > { x.val / y.val };
}

template<typename U>
Quantity<U> operator *(Quantity<U> x, double y) {
	return Quantity<U> { x.val * y };
}

template<typename U>
Quantity<U> operator *(double x, Quantity<U> y) {
	return Quantity<U> { x * y.val };
}
constexpr Quantity<M> operator"" _m(long double d) {
	return Quantity<M> { d };
}

constexpr Quantity<M> operator"" _nm(long double d) {
	return Quantity<M> { d / 1e9 };
}

constexpr Quantity<M> operator"" _cm(long double d) {
	return Quantity<M> { d / 1e2 };
}

constexpr Quantity<M> operator"" _A(long double d) {
	return Quantity<M> { d / 1e10 };
}

constexpr Quantity<Kg> operator"" _kg(long double d) {
	return Quantity<Kg> { d };
}

constexpr Quantity<Kg> operator"" _g(long double d) {
	return Quantity<Kg> { d / 1e3 };
}

constexpr Quantity<Kg> operator"" _mg(long double d) {
	return Quantity<Kg> { d / 1e6 };
}

constexpr Quantity<S> operator"" _s(long double d) {
	return Quantity<S> { d };
}

constexpr Quantity<S> operator"" _ms(long double d) {
	return Quantity<S> { d / 1e3 };
}

constexpr Quantity<S> operator"" _us(long double d) {
	return Quantity<S> { d / 1e6 };
}

constexpr Quantity<S> operator"" _ns(long double d) {
	return Quantity<S> { d / 1e9 };
}

constexpr Quantity<Joules> operator"" _J(long double d) {
	return Quantity<Joules> { d };
}

constexpr Quantity<Joules> operator"" _KJ(long double d) {
	return Quantity<Joules> { d * 1e3 };
}

constexpr Quantity<Farads> operator"" _Far(long double d) {
	return Quantity<Farads> { d };
}

constexpr Quantity<Adim> operator"" _adim(long double d) {
	return Quantity<Adim> { d };
}

constexpr Quantity<C> operator"" _C(long double d) {
	return Quantity<C> { d };
}

constexpr Quantity<A> operator"" _Amp(long double d) {
	return Quantity<A> { d };
}
constexpr Quantity<Volts> operator"" _Volt(long double d) {
	return Quantity<Volts> { d };
}
constexpr Quantity<C> operator"" _el(long double d) {
	return Quantity<C> { d * 1.6021766208e-19 };
}

constexpr Quantity<T> operator"" _K(long double d) {
	return Quantity<T> { d };
}


const auto eps0 = 8.854187817e-12_Far / 1.0_m;
const auto avogad = 6.0225e23_adim;
const auto R = 8.3144598e-3_KJ / 1.0_K;
const auto kb = R / avogad;
const auto kt300 = R * 300.0_K;

template<typename U>
constexpr Quantity<Unit_plus<U, U> > square(Quantity<U> x) {
	return Quantity<Unit_plus<U, U> > { x.val * x.val };
}

template<typename U>
constexpr bool operator==(Quantity<U> x, Quantity<U> y) {
	return x.val == y.val;
}

template<typename U>
constexpr bool operator!=(Quantity<U> x, Quantity<U> y) {
	return x.val != y.val;
}

inline std::string suffix(int u, const char *x) {
	std::string suf;
	if (u) {
		suf += x;
		if (u > 1) {
			suf += "^{";
			suf += '0' + u;
			suf += "}";
		}
		if (u < 0) {
			suf += "^{-";
			suf += '0' - u;
			suf += "}";
		}
	}
	return suf;
}
template<typename U>
std::ostream & operator<<(std::ostream & os, Quantity<U> v) {
	os << v.val << " " << suffix(U::kg, "kg") << suffix(U::m, "m")
			<< suffix(U::s, "s") << suffix(U::q, "C") << suffix(U::t, "K");
	return os;
}

}




#endif /* TRJLIB_UNITS_HH_ */
