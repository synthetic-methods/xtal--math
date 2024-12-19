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

template <int M_ism=1, int M_pow=1, int M_car=0> requires in_q<M_ism, 1, 2, -1, -2> and in_q<M_pow, 1, -1> and in_q<M_car, -0, -1>
struct   monologarithm;

template <auto ...Ms>
using    monologarithm_t = process::confined_t<monologarithm<Ms...>>;

template <auto ...Ms>
XTAL_DEF_(short)
XTAL_LET monologarithm_f(auto &&o, constant_q auto ...oo)
noexcept -> decltype(auto)
{
	return monologarithm_t<Ms...>::template function<oo...>(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` as the monologarithm `-Log[1 - #]`, \
approximated by `#/Sqrt[1 - #]`. \

template <int M_ism, int M_pow> requires in_q<M_ism, 1, 2>
struct monologarithm<M_ism, M_pow, -0>
{
	using superprocess = process::lift_t<void
	,	bond::compose<dilated<2>, taylor::sine<-2>>
	,	bond::compose<discarded<1, M_pow>, monologarithm<M_ism, M_pow, -1>>
	>;
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(o)>;
			auto constexpr _1 = _op::alpha_1;
			auto constexpr _i = _op::alpha_1*signum_n<M_ism&1, -1>;

			if constexpr (N_lim < 0) {
			//	XTAL_0IF (M_pow == 1) {return               _i*log(_1 + _i*XTAL_REF_(o)) ;}
			//	XTAL_0IF (M_pow == 1) {return root_f<M_pow>(_i*log(_1 + _i*XTAL_REF_(o)));}
				return root_f<M_pow>(_i*log(_1 + _i*XTAL_REF_(o)));
			}
			else {
				return superprocess::template function<N_lim>(XTAL_REF_(o));
			}
		}

	};
};
///\
Defines `function` as the antimonologarithm `1 - Exp[-#]`, \
approximated by `(Sqrt[1 + (#/2)^2] - (#/2))*(#)`. \

template <int M_ism, int M_pow> requires in_q<M_ism,-1,-2>
struct monologarithm<M_ism, M_pow, -0>
{
	using superprocess = process::lift_t<void
	,	bond::compose<discarded<1, M_pow>, monologarithm<M_ism, M_pow, -1>>
	,	bond::compose<dilated<2>, taylor::sine<+2>>
	>;
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(o)>;
			auto constexpr _1 = _op::alpha_1;
			auto constexpr _i = _op::alpha_1*signum_n<M_ism&1, -1>;

			if constexpr (N_lim < 0) {
				return root_f<M_pow>(term_f(-_i, _i, exp(_i*XTAL_REF_(o))));
			}
			else {
				return superprocess::template function<N_lim>(XTAL_REF_(o));
			}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` as the cardinal monologarithm `-Log[1 - #]/#`, \
approximated by `1/Sqrt[1 - #]`. \

template <int M_ism, int M_pow> requires in_q<M_ism, 1, 2>
struct monologarithm<M_ism, M_pow, -1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;
		using S0 = bond::compose_s<S, monologarithm<M_ism, M_pow, -0>>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&u)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(u)>;
			auto constexpr _1 = _op::alpha_1;
			auto constexpr _i = _op::alpha_1*signum_n<M_ism&1, -1>;

			if constexpr (N_lim < 0) {
//				echo("M_car<-1>");
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
Defines `function` as the cardinal antimonologarithm `(1 - Exp[-#])/#`, \
approximated by `Sqrt[1 + (#/2)^2] + (#/2)`. \

template <int M_ism, int M_pow> requires in_q<M_ism,-1,-2>
struct monologarithm<M_ism, M_pow, -1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;
		using S0 = bond::compose_s<S, monologarithm<M_ism, M_pow, -0>>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&u)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(u)>;
			auto constexpr _1 = _op::alpha_1;
			auto constexpr _i = _op::alpha_1*signum_n<M_ism&1, -1>;

			if constexpr (N_lim < 0) {
//				echo("M_car<-2>");
				XTAL_IF0
				XTAL_0IF (M_pow ==  1) {return S0::template function<-0, N_lim>(u)/XTAL_REF_(u);}
				XTAL_0IF (M_pow == -1) {return XTAL_REF_(u)/S0::template function<-0, N_lim>(u);}
			}
			else {
				auto v = XTAL_REF_(u)*_op::haplo_1;
				return root_f<M_pow>(root_f<2>(term_f(_1, v, v)) + _i*v);
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
