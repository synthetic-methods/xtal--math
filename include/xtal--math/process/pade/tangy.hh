#pragma once
#include "./any.hh"

#include "./impunity.hh"
#include "../gudermannian/tang.hh"



XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Defines `Tan[Pi #] &` and `Tanh[Pi #] &`. \

///\param M_ism \f$\in {1, 2}\f$ specifies the underlying morphism, \
generating either the circular or hyperbolic tangent. \

template <int M_ism=0, int M_car=0, typename ...As>// requires in_n<M_ism, 0, 1,-1, 2,-2> and in_n<M_car, 1,-0,-1,-2>
struct tangy
:	process::lift<tangy<M_ism, M_car>, bond::compose<As...>>
{
};
template <>
struct   tangy<>
{
	using limit_type = occur::math::limit_t<(1<<3)>;

	template <class S>
	using subtype = bond::compose_s<S, provision::context<void
	,	typename limit_type::template dispatch<>
	>>;

};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_n<M_ism, 1, 2>
struct tangy<M_ism,-0>
{
	static constexpr int I_sgn = sign_n<(M_ism&1)^1, -1>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		static_method(simplex_field_q auto &&o)
		noexcept -> auto
		{
			using _fix = bond::fixture<decltype(o)>;
			
			XTAL_IF0
			XTAL_0IF (N_lim <  0) {
				XTAL_IF0
				XTAL_0IF (1 == M_ism) {return tan (XTAL_REF_(o)*_fix::patio_1);}
				XTAL_0IF (2 == M_ism) {return tanh(XTAL_REF_(o)*_fix::patio_1);}
			}
			XTAL_0IF (N_lim == 0) {
				return _fix::patio_1*gudermannian::tang_t<M_ism>::template static_method<N_lim>(XTAL_REF_(o));
			}
			XTAL_0IF (0 == (N_lim&1)) {
				auto const [x1, y1] = destruct_f(impunity_t<M_ism,-0>::template static_method<N_lim>(o*_fix::haplo_1));
				auto const x2 =  square_f(x1) + I_sgn*square_f(y1);
				auto const y2 = _fix::diplo_1*x1*y1;
				return y2*root_f<-1, 1>(x2);
			}
			XTAL_0IF (1 == (N_lim&1)) {
				auto const [x1, y1] = destruct_f(impunity_t<M_ism,-0>::template static_method<N_lim>(o));
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
		static_method(simplex_field_q auto &&t)
		noexcept -> decltype(auto)
		{
			return static_method<N_lim>(XTAL_REF_(t), absolve_u<decltype(t)>{one});
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		static_method(simplex_field_q auto &&v, simplex_field_q auto &&u)
		noexcept -> decltype(auto)
		{
			using _fix = bond::fixture<decltype(v), decltype(u)>;
			using U_aphex = typename _fix::aphex_type;
			using U_alpha = typename _fix::alpha_type;
			using W_alpha = atom::couple_t<U_alpha[2]>;

			auto u_abs = u, u_sgn = signum_t<>::static_edit(u_abs);
			auto v_abs = v, v_sgn = signum_t<>::static_edit(v_abs);// v_sgn *= *_fix::haplo_1;

			W_alpha co(v_abs < u_abs);
			W_alpha up{v, u_abs}; up *= co;
			W_alpha dn{u_abs,-v}; dn *= co;

			auto const &[co_0, co_1]  = co;
			auto const u_flp = _fix::haplo_1 - _fix::haplo_1*co_0*u_sgn;

			return term_f(u_flp*v_sgn, u_sgn, S_::template static_method<N_lim>(up.sum()/dn.sum()));
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
	static constexpr int I_sgn = sign_n<(M_ism&1)^1, -1>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		static_method(simplex_field_q auto &&o)
		noexcept -> decltype(auto)
		{
			return rate_f<-1>(impunity_t<M_ism,-2>::template static_method<N_lim>(XTAL_REF_(o)));
		}

	};
};

template <int M_ism=1, int M_car=0, typename ...As>
using tangy_t = process::confined_t<tangy<M_ism, M_car, As...>>;//, tangy<>>;


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
