#pragma once
#include "./any.hh"

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

template <int M_ism> requires in_n<M_ism, 0, 1, 2>
struct unity<M_ism, 0>
{
	using superprocess = process::lift_t<unify<M_ism>, _detail::impunity<M_ism>>;

	using superkind = any<unity>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

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
		method_f(complex_field_q auto const &t)
		noexcept -> decltype(auto)
		{
			return method_f<N_lim>(t.real(), t.imag());
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
		method_f(simplex_field_q auto o)
		noexcept -> decltype(auto)
		{
			using U       = XTAL_ALL_(o);
			using U_fit   = bond::fit<U>;
			using U_alpha = typename U_fit::alpha_type;

			if constexpr (N_lim < 0) {
				auto w = o*U_fit::patio_2;
				XTAL_IF0
				XTAL_0IF (1 == M_ism) {return complexion_f(cos (w), sin (w));}
				XTAL_0IF (2 == M_ism) {return complexion_f(cosh(w), sinh(w));}
			}
			else {
				auto const  w = wrap_f(o);
				auto const  m = objective_f(U_fit::haplo_1*wrap_f(U_fit::diplo_1*w));
				auto const _1 = U_fit::alpha_1*(((m == w) << 1) - 1);
				return superprocess::template method_f<N_lim>(m)*_1;
			}
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,set)
		method_f(atom::math::phason_q auto t_)
		noexcept -> decltype(auto)
		{
			using T_ = XTAL_ALL_(t_);
			using T_fit = bond::fit<decltype(t_[0])>;// Underlying...
			using U_fit = bond::fit<decltype(t_(0))>;
			using T_sigma = typename T_fit::sigma_type;
			using U_sigma = typename U_fit::sigma_type;
			using U_alpha = typename U_fit::alpha_type;

			if constexpr (N_lim < 0) {
				return method_f<N_lim>(t_(0));
			}
			else {
				T_sigma &u0 = t_[0];
				T_sigma  un = u0;
				un ^= un << 1;
				un &= T_fit::sign.mask;
				u0 ^= un;
				U_sigma vn{un};
				vn <<= U_fit::full.depth - T_fit::full.depth;
				vn  |= U_fit::unit.mask;
				return superprocess::template method_f<N_lim>(t_(0))*
					operative_f<[] XTAL_1FN_(call) (_xtd::bit_cast<U_alpha>)>(vn);
			}
		}

	};
};
template <int M_ism> requires in_n<M_ism, -1, -2>
struct unity<M_ism, 0>
{
	using superkind = any<unity>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method_f(complex_field_q auto &&o)
		noexcept -> decltype(auto)
		{
			using U       = XTAL_ALL_(o);
			using U_fit    = bond::fit<U>;

			auto constexpr N_lim_tan = N_lim;
			auto constexpr N_lim_log = 2;

			auto constexpr _1 = U_fit::haplo_0, _2pi = _1/U_fit::patio_2;
			auto constexpr _2 = U_fit::haplo_1, _4pi = _2/U_fit::patio_2;

			auto const &[x_re, x_im] = destruct_f(o);
			auto const   y_re = tangy_t<-1, 1>::template method_f<N_lim_tan>(x_im, x_re);
			auto const   w_im = square_f(x_re, x_im);
			//\
			auto const   y_im = log(w_im);
			auto const   y_im = taylor::logarithm_t< 1, 1>::template method_f<N_lim_log>(w_im);
			return complexion_f(y_re*_2, y_im*-_4pi);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_car=0>
using unity_t = process::confined_t<unity<M_ism, M_car>>;

template <int M_ism=1, int M_car=0, int ...Ns>
XTAL_DEF_(let)
unity_f = [] XTAL_1FN_(call) (unity_t<M_ism, M_car>::template method_f<Ns...>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
