#pragma once
#include "./any.hh"

#include "./sine.hh"
#include "../dilated.hh"
#include "../square.hh"
#include "../root.hh"

XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_pow=1, int M_car=0>
XTAL_REQ inclusive_q<M_ism, 1, 2, -1, -2> and inclusive_q<M_pow, 1, -1> and inclusive_q<M_car, -0, -1>
XTAL_TYP polylogarithm;

template <int M_ism=1, int M_pow=1, int M_car=0>
XTAL_USE polylogarithm_t = process::confined_t<polylogarithm<M_ism, M_pow, M_car>>;

template <int M_ism=1, int M_pow=1, int M_car=0, int ...Ns>
XTAL_DEF_(return,inline)
XTAL_FN1 polylogarithm_f(auto &&o)
XTAL_0EX
{
	return polylogarithm_t<M_ism, M_pow, M_car>::template function<Ns...>(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` as the polylogarithm `-Log[1 - #]`, \
approximated by `#/Sqrt[1 - #]`. \

template <int M_ism, int M_pow> requires inclusive_q<M_ism, 1, 2>
struct polylogarithm<M_ism, M_pow, -0>
{
	using subprocess = process::link_t<void
	,	bond::compose<dilated<1>, taylor::sine<-2>>
	,	bond::compose<discarded<M_pow, +1>, polylogarithm<M_ism, M_pow, -1>>
	>;

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline)
		XTAL_FN1 function(auto &&o)
		XTAL_0EX
		{
			using _std::log;

			using _op = bond::operate<decltype(o)>;
			auto constexpr _1 = _op::alpha_1;
			auto constexpr _i = _op::alpha_1*-sign_n<M_ism&1, -1>;

			if constexpr (N_lim < 0) {
				return root_f<M_pow>(_i*log(_1 + _i*XTAL_REF_(o)));
			}
			else {
				return subprocess::template function<N_lim>(XTAL_REF_(o));
			}
		}

	};
};
///\
Defines `function` as the antipolylogarithm `1 - Exp[-#]`, \
approximated by `(Sqrt[1 + (#/2)^2] - (#/2))*(#)`. \

template <int M_ism, int M_pow> requires inclusive_q<M_ism,-1,-2>
struct polylogarithm<M_ism, M_pow, -0>
{
	using subprocess = process::link_t<void
	,	bond::compose<discarded<M_pow, +1>, polylogarithm<M_ism, M_pow, -1>>
	,	bond::compose<dilated<1>, taylor::sine<+2>>
	>;

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline)
		XTAL_FN1 function(auto &&o)
		XTAL_0EX
		{
			using _std::exp;

			using _op = bond::operate<decltype(o)>;
			auto constexpr _1 = _op::alpha_1;
			auto constexpr _i = _op::alpha_1*-sign_n<M_ism&1, -1>;

			if constexpr (N_lim < 0) {
				return root_f<M_pow>(exp(XTAL_REF_(o)*_i)*_i - _i);
			}
			else {
				return subprocess::template function<N_lim>(XTAL_REF_(o));
			}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` as the cardinal polylogarithm `-Log[1 - #]/#`, \
approximated by `1/Sqrt[1 - #]`. \

template <int M_ism, int M_pow> requires inclusive_q<M_ism, 1, 2>
struct polylogarithm<M_ism, M_pow, -1>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;
		using S0 = bond::compose_s<S, polylogarithm<M_ism, M_pow, -0>>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline)
		XTAL_FN1 function(auto &&u)
		XTAL_0EX
		{
			using _op = bond::operate<decltype(u)>;
			auto constexpr _1 = _op::alpha_1;
			auto constexpr _i = _op::alpha_1*-sign_n<M_ism&1, -1>;

			if constexpr (N_lim < 0) {
				XTAL_IF0
				XTAL_0IF (M_pow ==  1) {return S0::template function<-0, N_lim>(u)/XTAL_REF_(u);}
				XTAL_0IF (M_pow == -1) {return XTAL_REF_(u)/S0::template function<-0, N_lim>(u);}
			}
			else {
				return root_f<M_pow*-2>(_1 + _i*XTAL_REF_(u));
			}
		}

	};
};
///\
Defines `function` as the cardinal antipolylogarithm `(1 - Exp[-#])/#`, \
approximated by `Sqrt[1 + (#/2)^2] + (#/2)`. \

template <int M_ism, int M_pow> requires inclusive_q<M_ism,-1,-2>
struct polylogarithm<M_ism, M_pow, -1>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;
		using S0 = bond::compose_s<S, polylogarithm<M_ism, M_pow, -0>>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline)
		XTAL_FN1 function(auto &&u)
		XTAL_0EX
		{
			using _op = bond::operate<decltype(u)>;
			auto constexpr _1 = _op::alpha_1;
			auto constexpr _i = _op::alpha_1*-sign_n<M_ism&1, -1>;

			if constexpr (N_lim < 0) {
				XTAL_IF0
				XTAL_0IF (M_pow ==  1) {return S0::template function<-0, N_lim>(u)/XTAL_REF_(u);}
				XTAL_0IF (M_pow == -1) {return XTAL_REF_(u)/S0::template function<-0, N_lim>(u);}
			}
			else {
				auto v = XTAL_REF_(u)*_op::haplo_1;
				return root_f<M_pow>(root_f<2>(horner::term_f<1>(_1, v, v)) + _i*v);
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
