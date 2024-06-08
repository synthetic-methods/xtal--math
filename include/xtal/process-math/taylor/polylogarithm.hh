#pragma once
#include "./any.hh"

#include "./sine.hh"
#include "../dilating.hh"
#include "../square.hh"
#include "../root.hh"

XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_car=0>
XTAL_REQ inclusive_q<M_ism, 1, 2,-1,-2> and inclusive_q<M_car,-0,-1>
XTAL_TYP polylogarithm;

template <int M_ism=1, int M_car=0>
XTAL_USE polylogarithm_t = process::confined_t<polylogarithm<M_ism, M_car>>;

template <int M_ism=1, int M_car=0>
XTAL_FN2 polylogarithm_f(auto &&o)
XTAL_0EX
{
	return polylogarithm_t<M_ism, M_car>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` as the polylogarithm `-Log[1 - #]`, \
approximated by `#/Sqrt[1 - #]`. \

template <int M_ism> requires inclusive_q<M_ism, 1, 2>
struct polylogarithm<M_ism, -0>
{
	using subprocess = process::link_t<void
	,	bond::compose<dilating<+1>, taylor::sine<-2>>
	,	bond::compose<discarding<1>, polylogarithm<M_ism, -1>>
	>;

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_FN2 function(auto &&o)
		XTAL_0EX
		{
			using _std::log;

			using _op = bond::operate<decltype(o)>;
			auto constexpr _1 = _op::alpha_1;
			auto constexpr _i = _op::alpha_1*sign_n<(M_ism&1)^1, -1>;

			if constexpr (N_lim < 0) {
				return _i*log(_1 + _i*XTAL_REF_(o));
			}
			else {
				return subprocess::template function<N_lim>(typename _op::alpha_t(XTAL_REF_(o)));
			}
		}

	};
};
///\
Defines `function` as the antipolylogarithm `1 - Exp[-#]`, \
approximated by `(Sqrt[1 + (#/2)^2] - (#/2))*(#)`. \

template <int M_ism> requires inclusive_q<M_ism,-1,-2>
struct polylogarithm<M_ism, -0>
{
	using subprocess = process::link_t<void
	,	bond::compose<discarding<1>, polylogarithm<M_ism, -1>>
	,	bond::compose<dilating<+1>, taylor::sine<+2>>
	>;

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_FN2 function(auto &&o)
		XTAL_0EX
		{
			using _std::exp;

			using _op = bond::operate<decltype(o)>;
			auto constexpr _1 = _op::alpha_1;
			auto constexpr _i = _op::alpha_1*sign_n<(M_ism&1)^1, -1>;

			if constexpr (N_lim < 0) {
				return exp(XTAL_REF_(o)*_i)*_i - _i;
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

template <int M_ism> requires inclusive_q<M_ism, 1, 2>
struct polylogarithm<M_ism, -1>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_FN2 function(auto &&u)
		XTAL_0EX
		{
			using _op = bond::operate<decltype(u)>;
			auto constexpr _1 = _op::alpha_1;
			auto constexpr _i = _op::alpha_1*sign_n<(M_ism&1)^1, -1>;

			if constexpr (N_lim < 0) {
				return polylogarithm_t<M_ism, -0>::template function<-1>(u)/XTAL_REF_(u);
			}
			else {
				return root_f<-2>(_1 + _i*XTAL_REF_(u));
			}
		}

	};
};
///\
Defines `function` as the cardinal antipolylogarithm `(1 - Exp[-#])/#`, \
approximated by `Sqrt[1 + (#/2)^2] + (#/2)`. \

template <int M_ism> requires inclusive_q<M_ism,-1,-2>
struct polylogarithm<M_ism, -1>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_FN2 function(auto &&o)
		XTAL_0EX
		{
			using _op = bond::operate<decltype(o)>;
			auto constexpr _1 = _op::alpha_1;
			auto constexpr _i = _op::alpha_1*sign_n<(M_ism&1)^1, -1>;

			if constexpr (N_lim < 0) {
				return polylogarithm_t<M_ism, -0>::template function<-1>(o)/XTAL_REF_(o);
			}
			else {
				auto u = XTAL_REF_(o)*_op::haplo_f(1);
				return root_f< 2>(horner::term_f<1>(_1, u, u)) + _i*u;
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
