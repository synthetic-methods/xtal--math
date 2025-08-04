#pragma once
#include "./any.hh"

#include "./sine.hh"




XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines a class of tangent approximations indexed by `M_ism`.

\tparam  M_ism Indexes circular/hyperbolic functions when `M_ism={1,2}`,
         and their respective inverses when `M_ism={-1,-2}`.

\note    Derives the tangent from `sine` using `(#/Sqrt[1 + #^2] &)`.
*/
template <int M_ism=0, int M_car=0>
XTAL_TYP_(new) tangent;

template <int M_ism=0, int M_car=0>
XTAL_TYP_(let) tangent_t = process::confined_t<tangent<M_ism, M_car>>;

template <int M_ism=0, int M_car=0, int N_lim=4>
XTAL_DEF_(let) tangent_f = [] XTAL_1FN_(call) (tangent_t<M_ism, M_car>::template method_f<N_lim>);


////////////////////////////////////////////////////////////////////////////////

template <int M_ism>
struct tangent<M_ism, -0>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0> requires (N_lim <  0)
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (M_ism ==  2) {return  tanh(XTAL_REF_(o));}
			XTAL_0IF (M_ism ==  1) {return  tan (XTAL_REF_(o));}
			XTAL_0IF (M_ism == -1) {return atan (XTAL_REF_(o));}
			XTAL_0IF (M_ism == -2) {return atanh(XTAL_REF_(o));}
		}
		template <int N_lim=0> requires (0 <= N_lim)
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
		{
			auto const o_sgn = decompose_f<signed>(o);
			auto x = objective_f(XTAL_REF_(o));
			XTAL_IF0
			XTAL_0IF (0 < M_ism) {
				x *= sine_t<M_ism, -1>::template method_f<N_lim>(x);
				x *= root_f<-2>(term_f(one,  cosign_v<M_ism>, square_f(x)));
			}
			XTAL_0IF (M_ism < 0) {
				x *= root_f<-2>(term_f(one, -cosign_v<M_ism>, square_f(x)));
				x *= sine_t<M_ism, -1>::template method_f<N_lim>(x);
			}
			return x;
		}

	};
};
template <int M_ism>
struct tangent<M_ism, -1>
{
	using superkind = tangent<M_ism, -0>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto const &o)
		noexcept -> decltype(auto)
		{
			return S_::template method_f<Ns...>(o)*root_f<-1, 1>(o);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
