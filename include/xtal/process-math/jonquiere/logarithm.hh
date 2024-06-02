#pragma once
#include "./any.hh"

#include "../dilate.hh"
#include "../square.hh"
#include "../taylor/sine.hh"


XTAL_ENV_(push)
namespace xtal::process::math::jonquiere
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_car=0> struct logarithm {static_assert(M_ism);};
template <int M_ism=1, int M_car=0> using  logarithm_t = process::confined_t<logarithm<M_ism, M_car>>;
template <int M_ism=1, int M_car=0>
XTAL_FN2 logarithm_f(auto &&o)
XTAL_0EX
{
	return logarithm_t<M_ism, M_car>::function(XTAL_REF_(o));
}

////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` as the polylogarithm `-Log[1 - #]`, \
approximated by `#/Sqrt[1 - #]`. \

template <int M_ism> requires (0 < M_ism)
struct logarithm<M_ism, -0>
{
	using subkind = process::link<void
	,	bond::compose<dilate<+1>, taylor::sine<-2>>
	,	bond::compose<discard<1>, logarithm<M_ism, -1>>
	>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_FN2 function(auto &&o)
		XTAL_0EX
		{
			using _std::log;

			using op = bond::operate<decltype(o)>;
			auto constexpr _1 = op::alpha_1;
			auto constexpr _i = op::alpha_1*sign_n<(M_ism&1)^1, -1>;

			if constexpr (N_lim < 0) {
				return _i*log(_1 + _i*XTAL_REF_(o));
			}
			else {
				return S_::template function<N_lim>(typename op::alpha_t(XTAL_REF_(o)));
			}
		}

	};
};
///\
Defines `function` as the antipolylogarithm `1 - Exp[-#]`, \
approximated by `(Sqrt[1 + (#/2)^2] - (#/2))*(#)`. \

template <int M_ism> requires (M_ism < 0)
struct logarithm<M_ism, -0>
{
	using subkind = process::link<void
	,	bond::compose<discard<1>, logarithm<M_ism, -1>>
	,	bond::compose<dilate<+1>, taylor::sine<+2>>
	>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_FN2 function(auto &&o)
		XTAL_0EX
		{
			using _std::exp;

			using op = bond::operate<decltype(o)>;
			auto constexpr _1 = op::alpha_1;
			auto constexpr _i = op::alpha_1*sign_n<(M_ism&1)^1, -1>;

			if constexpr (N_lim < 0) {
				return exp(XTAL_REF_(o)*_i)*_i - _i;
			}
			else {
				return S_::template function<N_lim>(XTAL_REF_(o));
			}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` as the cardinal polylogarithm `-Log[1 - #]/#`, \
approximated by `1/Sqrt[1 - #]`. \

template <int M_ism> requires (0 < M_ism)
struct logarithm<M_ism, -1>
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
			using op = bond::operate<decltype(u)>;
			auto constexpr _1 = op::alpha_1;
			auto constexpr _i = op::alpha_1*sign_n<(M_ism&1)^1, -1>;

			if constexpr (N_lim < 0) {
				return logarithm_t<M_ism, -0>::template function<-1>(u)/XTAL_REF_(u);
			}
			else {
				return square_f<-1,-1>(_1 + _i*XTAL_REF_(u));
			}
		}

	};
};
///\
Defines `function` as the cardinal antipolylogarithm `(1 - Exp[-#])/#`, \
approximated by `Sqrt[1 + (#/2)^2] + (#/2)`. \

template <int M_ism> requires (M_ism < 0)
struct logarithm<M_ism, -1>
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
			using op = bond::operate<decltype(o)>;
			auto constexpr _1 = op::alpha_1;
			auto constexpr _i = op::alpha_1*sign_n<(M_ism&1)^1, -1>;

			if constexpr (N_lim < 0) {
				return logarithm_t<M_ism, -0>::template function<-1>(o)/XTAL_REF_(o);
			}
			else {
				auto u = XTAL_REF_(o)*op::haplo_f(1);
				return square_f<-1>(horner::term_f<1>(_1, u, u)) + _i*u;
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
