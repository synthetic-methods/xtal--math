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

template <int M_ism=1, int M_car=0> requires in_v<M_ism, 1, 2, -1, -2> and in_v<M_car, -0, -1>
struct  monologarithm;

template <auto ...Ms>
using   monologarithm_t = process::confined_t<monologarithm<Ms...>>;


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines `function` as the monologarithm `-Log[1 - #]`.

Approximated by `#/Sqrt[1 - #]`.
*/

template <int M_ism> requires in_v<M_ism, 1, 2>
struct monologarithm<M_ism, -0>
{
	static int constexpr I_ism = M_ism&1;

	using superprocess = process::lift_t<void
	,	bond::compose<dilated<2>, taylor::sine<-2>>
	,	bond::compose<discarded<1>, monologarithm<M_ism, -1>>
	>;
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			auto constexpr k = _fit::alpha_f(sign_v<I_ism, -1>);

			XTAL_IF0
			XTAL_0IF (0 <= N_lim) {return  superprocess::template method_f<N_lim>(XTAL_REF_(o));}
			XTAL_0IF (0 == I_ism) {return -log(one - XTAL_REF_(o));}
			XTAL_0IF (1 == I_ism) {return  log(one + XTAL_REF_(o));}
		}

	};
};
/*!
\brief   Defines `function` as the antimonologarithm `1 - Exp[-#]`,

Approximated by `(Sqrt[1 + (#/2)^2] - (#/2))*(#)`.
*/

template <int M_ism> requires in_v<M_ism,-1,-2>
struct monologarithm<M_ism, -0>
{
	static int constexpr I_ism = M_ism&1;

	using superprocess = process::lift_t<void
	,	bond::compose<discarded<1>, monologarithm<M_ism, -1>>
	,	bond::compose<dilated<2>, taylor::sine<+2>>
	>;
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> auto
		{
			using U = XTAL_ALL_(o);
			XTAL_IF0
			XTAL_0IF (0 <= N_lim) {return superprocess::template method_f<N_lim>(XTAL_REF_(o));}
#if XTAL_SYS_(builtin)
			XTAL_0IF (0 == I_ism and real_q<U>) {return -__builtin_expm1(-XTAL_REF_(o));}
			XTAL_0IF (1 == I_ism and real_q<U>) {return  __builtin_expm1( XTAL_REF_(o));}
#endif
			XTAL_0IF (0 == I_ism) {return one - exp(-XTAL_REF_(o));}
			XTAL_0IF (1 == I_ism) {return exp( XTAL_REF_(o)) - one;}
			XTAL_0IF_(void)
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines `function` as the cardinal monologarithm `-Log[1 - #]/#`,

Approximated by `1/Sqrt[1 - #]`.
*/
template <int M_ism> requires in_v<M_ism, 1, 2>
struct monologarithm<M_ism, -1>
{
	static int constexpr I_ism = M_ism&1;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;
		using S0 = bond::compose_s<S, monologarithm<M_ism, -0>>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (N_lim <  0) {return S0::template method_f<-0, N_lim>(o)/o;}
			XTAL_0IF (0 == I_ism) {return root_f<-2>(one - XTAL_REF_(o));}
			XTAL_0IF (1 == I_ism) {return root_f<-2>(one + XTAL_REF_(o));}
		}

	};
};
/*!
\brief   Defines `function` as the cardinal antimonologarithm.
\tparam  `M_ism` Selects between `(Exp@# - 1)/#&` and `(1 - Exp@-#)/#&` when `M_ism=-1` and `M_ism=-2` respectively.

When `M_ism=-1`, approximated by `Sqrt[1 + (#/2)^2] + (#/2)`.
*/
template <int M_ism> requires in_v<M_ism,-1,-2>
struct monologarithm<M_ism, -1>
{
	static int constexpr I_ism = M_ism&1;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;
		using S0 = bond::compose_s<S, monologarithm<M_ism, -0>>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> auto
		{
			auto const u = half*XTAL_REF_(o);
			XTAL_IF0
			XTAL_0IF (N_lim <  0) {return S0::template method_f<-0, N_lim>(o)/o;}
			XTAL_0IF (0 == I_ism) {return root_f<2>(term_f(one, u, u)) - u;}
			XTAL_0IF (1 == I_ism) {return root_f<2>(term_f(one, u, u)) + u;}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <auto N>
XTAL_DEF_(return,inline,let)
monologarithm_f(auto &&...oo)
noexcept -> decltype(auto)
{
	auto constexpr N_sgn =   sign_v<N>;
	auto constexpr N_abs = N*sign_v<N>;
	return monologarithm_t<1*N_sgn>::template method_f<N_abs>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
