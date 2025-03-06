#pragma once
#include "./any.hh"

#include "../gudermannian/tang.hh"




XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines `Tan[Pi #] &` and `Tanh[Pi #] &`.

\tparam  M_ism
Specifies the underlying morphism \f$\in {1, 2}\f$,
generating either the circular or hyperbolic tangent.
*/
template <int M_ism=0, int M_car=0>
struct tangy;


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_n<M_ism, 1, 2>
struct tangy<M_ism,-0>
{
	static constexpr int I_sgn = cosign_v<M_ism>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(simplex_field_q auto &&o)
		noexcept -> auto
		{
			using _fit = bond::fit<decltype(o)>;
			
			XTAL_IF0
			XTAL_0IF (N_lim <  0) {
				XTAL_IF0
				XTAL_0IF (1 == M_ism) {return tan (XTAL_REF_(o)*_fit::patio_1);}
				XTAL_0IF (2 == M_ism) {return tanh(XTAL_REF_(o)*_fit::patio_1);}
			}
			XTAL_0IF (N_lim == 0) {
				return _fit::patio_1*gudermannian::tang_t<M_ism>::template method_f<N_lim>(XTAL_REF_(o));
			}
			XTAL_0IF (0 == (N_lim&1)) {
				auto const [x1, y1] = destruct_f(_detail::impunity_t<M_ism,-0>::template method_f<N_lim>(o*_fit::haplo_1));
				auto const x2 =  square_f(x1) + I_sgn*square_f(y1);
				auto const y2 = _fit::diplo_1*x1*y1;
				return y2*root_f<-1, 1>(x2);
			}
			XTAL_0IF (1 == (N_lim&1)) {
				auto const [x1, y1] = destruct_f(_detail::impunity_t<M_ism,-0>::template method_f<N_lim>(o));
				return y1*root_f<-1, 1>(x1);
			}
		}

	};
};
template <int M_ism> requires in_n<M_ism,-1,-2>
struct tangy<M_ism, 1>
{
	using superkind = tangy<M_ism,-0>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(simplex_field_q auto &&t)
		noexcept -> decltype(auto)
		{
			return method_f<N_lim>(XTAL_REF_(t), unstruct_u<decltype(t)>{one});
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(simplex_field_q auto &&v, simplex_field_q auto &&u)
		noexcept -> decltype(auto)
		{
			using _fit = bond::fit<decltype(v), decltype(u)>;
			using U_aphex = typename _fit::aphex_type;
			using U_alpha = typename _fit::alpha_type;
			using W_alpha = atom::couple_t<U_alpha[2]>;

			auto u_abs = u, u_sgn = aspect_t<signed>::edit_f(u_abs);
			auto v_abs = v, v_sgn = aspect_t<signed>::edit_f(v_abs);// v_sgn *= *_fit::haplo_1;

			W_alpha co(v_abs < u_abs);
			W_alpha up{v, u_abs}; up *= co;
			W_alpha dn{u_abs,-v}; dn *= co;

			auto const &[co_0, co_1]  = co;
			auto const u_flp = _fit::haplo_1 - _fit::haplo_1*co_0*u_sgn;

			return term_f(u_flp*v_sgn, u_sgn, S_::template method_f<N_lim>(up.sum()/dn.sum()));
		}

	};
};
template <int M_ism>
struct tangy<M_ism,-0>
:	bond::compose<discarded<1>, tangy<M_ism,-1>>
{
};
template <int M_ism>
struct tangy<M_ism,-1>
:	bond::compose<discarded<2>, tangy<M_ism,-2>>
{
};
template <int M_ism>// requires in_n<M_ism,-1>
struct tangy<M_ism,-2>
{
	static constexpr int I_sgn = cosign_v<M_ism>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(simplex_field_q auto &&o)
		noexcept -> decltype(auto)
		{
			auto const [t_re, t_im] = destruct_f(_detail::impunity_t<M_ism,-2>::template method_f<N_lim>(XTAL_REF_(o)));
			return t_im*root_f<-1, 1>(t_re);
		}

	};
};

///////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_car=0>
using tangy_t = process::confined_t<tangy<M_ism, M_car>>;

template <int M_ism=1, int M_car=0, int ...Ns>
XTAL_DEF_(let)
tangy_f = [] XTAL_1FN_(call) (tangy_t<M_ism, M_car>::template method_f<Ns...>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
