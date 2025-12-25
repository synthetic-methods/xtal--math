#pragma once
#include "./any.hh"
#include "./arc.hh"
#include "./unify.hh"
#include "./tangy.hh"
#include "../taylor/logarithm.hh"
#include "../taylor/octarithm.hh"

XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines `function` by `(-1)^(2 #) &`; spiritually equivalent to `1^# &`.

\tparam M_ism
Specifies the underlying morphism \f$\in {1, 2}\f$,
generating either circular or hyperbolic `{cosine, sine}` pairs.
*/
template <int M_ism=1, int M_car=0>
struct unity;


////////////////////////////////////////////////////////////////////////////////

template <int M_ism, int M_car> requires in_n<M_ism, 0>
struct unity<M_ism, M_car>
:	identity<M_ism, M_car>
{
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism, int M_car> requires in_n<M_ism, 1, 2> and in_n<M_car, 0, 1>
struct unity<M_ism, M_car>
{
	using superprocess = process::lift_t<unify<M_ism>, _detail::impunity<M_ism>>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1, class U>
		XTAL_DEF_(return,inline,set)
		method_f(_std::initializer_list<U> o)
		noexcept -> decltype(auto)
		{
			using _fit = bond::fit<decltype(o)>;
			_std::complex<U> w; auto &m = destruct_f(w);
			_std::copy_n(point_f(o), 2, m);
			return method_f<N_lim>(w);
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(complex_field_q auto const &u)
		noexcept -> decltype(auto)
		{
			return method_f<N_lim>(u.real(), u.imag());
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&t_re, simplex_field_q auto &&t_im)
		noexcept -> decltype(auto)
		{
			auto constexpr _exp_2pi = [] XTAL_1FN_(call) (taylor::octarithm_t<-1>::template method_f<2>);
			return method_f<N_lim>(XTAL_REF_(t_re))*_exp_2pi(XTAL_REF_(t_im));
		}

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(simplex_field_q auto const &u)
		noexcept -> decltype(auto)
		{
			using U       = XTAL_ALL_(u);
			using U_fit   = bond::fit<U>;
			using U_alpha = typename U_fit::alpha_type;

			XTAL_IF0
			XTAL_0IF (M_car == 1) {
				return confined_t<unity<M_ism, 0>>::
					template method_f<N_lim>(modulo_f<-0>(u));
			}
			XTAL_0IF (0 <= N_lim) {
				auto const o_1 = objective_f(u);
				auto const o_2 = objective_f(modulo_f<-1>(o_1));
				return superprocess::template method_f<N_lim>(o_2)*
					term_f(one, -U_fit::diplo_f(2), decompose_f<unsigned>(o_1 - o_2));
			}
			XTAL_0IF (N_lim <  0) {
				XTAL_IF1_(consteval) {
					return method_f<below_v<(1<<2), (unsigned) N_lim>>(u);
				}
				XTAL_0IF_(else) {
					auto const w = objective_f(u*U_fit::patio_2);
					XTAL_IF0
					XTAL_0IF (1 == M_ism) {return complexion_f(cos (w), sin (w));}
					XTAL_0IF (2 == M_ism) {return complexion_f(cosh(w), sinh(w));}
				}
			}
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,set)
		method_f(atom::math::phason_q auto const &t_)
		noexcept -> decltype(auto)
		{
			using T_      = XTAL_ALL_(t_);
			using T_fit   = bond::fit<decltype(t_[0])>;// Underlying...
			using U_fit   = bond::fit<decltype(t_(0))>;
			using U_sigma = typename U_fit::sigma_type;
			using U_alpha = typename U_fit::alpha_type;

			XTAL_IF0
			XTAL_0IF (N_lim <  0) {
				return method_f<N_lim>(t_(0));
			}
			XTAL_0IF_(else) {
				auto f_  = modulo_f<-1>(t_);
				auto n1  = operative_f<[] XTAL_1FN_(call) (U_sigma)>(f_[0]^t_[0]);
				n1     <<= U_fit::full.depth - T_fit::full.depth;
				n1      &= U_fit::sign.mask;
				n1      |= U_fit::unit.mask;
				auto const u1 = operative_f<[] XTAL_1FN_(call) (_xtd::bit_cast<U_alpha>)>(n1);
				return superprocess::template method_f<N_lim>(f_(0))*u1;
			}
		}

	};
};
template <int M_ism, int M_car> requires in_n<M_ism,-1,-2> and in_n<M_car, 0, 1>
struct unity<M_ism, M_car>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(complex_field_q auto &&u)
		noexcept -> decltype(auto)
		{
		//	TODO: Handle `M_car == 1`?
			auto const   x_re = u.real();
			auto const   x_im = u.imag();
			auto const   y_re = half*arc_t<-0, M_car>::template method_f<N_lim>(x_im, x_re);
			auto const   y_im = half*arc_t<-1, M_car>::template method_f<N_lim>(square_f(x_re, x_im));
			return complexion_f(y_re, y_im);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_car=1>
XTAL_TYP_(let) unity_t = process::confined_t<unity<M_ism, M_car>>;

template <int M_ism=1, int M_car=1, int N_lim=2>
XTAL_DEF_(let) unity_f = [] XTAL_1FN_(call) (unity_t<M_ism, M_car>::template method_f<N_lim>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////


namespace xtal::occur
{////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <auto ..._s>
struct context<process::math::pade::unity<_s...>>
{
	using superkind = context<>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;
	
	public:
		using S_::S_;

		using order_attribute = occur::inferred_t<union ORDER, bond::seek_s<(1<<3)>>;

		template <extent_type N_mask=1>
		struct dispatch
		{
			template <class R>
			using subtype = bond::compose_s<R, typename S_::template dispatch<N_mask>
			,	provision::voiced<void
				,	typename T_::order_attribute::template dispatch<N_mask>
				>
			>;

		};

	};
};


///////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
