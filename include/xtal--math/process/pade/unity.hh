#pragma once
#include "./any.hh"

#include "./impunity.hh"
#include "./disunity.hh"
#include "./tangy.hh"
#include "../taylor/logarithm.hh"


XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` by `(-1)^(2 #) &`; spiritually equivalent to `1^# &`. \

///\param M_ism \f$\in {1, 2}\f$ specifies the underlying morphism, \
generating either circular or hyperbolic `{cosine, sine}` pairs. \

template <int M_ism=0, typename ...As> requires in_n<M_ism, 0, 1,-1, 2,-2>
struct unity
:	process::lift<unity<M_ism>, bond::compose<As...>>
{
};
template <>
struct   unity<>
{
	using limit_type = occur::math::limit_t<(1<<3)>;

	template <class S>
	using subtype = bond::compose_s<S, provision::context<void
	,	typename limit_type::template dispatch<>
	>>;

};
template <int M_ism=1, typename ...As>
using    unity_t = process::confined_t<unity<M_ism, As...>, unity<>>;


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_n<M_ism, 0, 1, 2>
struct unity<M_ism> : unity<>
{
	using superprocess = process::lift_t<disunity<M_ism>, impunity<M_ism>>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1, class U>
		XTAL_DEF_(return,inline,set)
		static_method(_std::initializer_list<U> o)
		noexcept -> decltype(auto)
		{
			using _fix = bond::fixture<decltype(o)>;
			_std::complex<U> w; auto &m = destruct_f(w);
			_std::copy_n(point_f(o), 2, m);
			return static_method<N_lim>(w);
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		static_method(complex_field_q auto const &t)
		noexcept -> decltype(auto)
		{
			return static_method<N_lim>(t.real(), t.imag());
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		static_method(auto &&t_re, simplex_field_q auto &&t_im)
		noexcept -> decltype(auto)
		{
			auto constexpr exp = [] XTAL_0FN_(alias) (taylor::logarithm_t<-1, 1>::template static_method<2>);
			using U_fix = bond::fixture<decltype(t_re), decltype(t_im)>;
			return static_method<N_lim>(XTAL_REF_(t_re))*exp(XTAL_REF_(t_im)*U_fix::patio_f(-2));
		}

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		static_method(simplex_field_q auto o)
		noexcept -> decltype(auto)
		{
			using U       = XTAL_ALL_(o);
			using U_fix    = bond::fixture<U>;
			using U_alpha = typename U_fix::alpha_type;

			if constexpr (N_lim < 0) {
				auto w = o*U_fix::patio_2;
				XTAL_IF0
				XTAL_0IF (1 == M_ism) {return complexion_f(cos (w), sin (w));}
				XTAL_0IF (2 == M_ism) {return complexion_f(cosh(w), sinh(w));}
			}
			else {
				auto w = wrap_f(o);
				auto m = objective_f(wrap_f(w*U_fix::diplo_1)*U_fix::haplo_1);
				return superprocess::template static_method<N_lim>(m)*
					operative_f<[] XTAL_0FN_(alias) (U_alpha)>(((m == w) << 1) - 1);
			}
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,set)
		static_method(arrange::math::phason_q auto t_)
		noexcept -> decltype(auto)
		{
			using T_ = XTAL_ALL_(t_);
			using T_fix = bond::fixture<decltype(t_[0])>;// Underlying...
			using U_fix = bond::fixture<decltype(t_(0))>;
			using T_sigma = typename T_fix::sigma_type;
			using U_sigma = typename U_fix::sigma_type;
			using U_alpha = typename U_fix::alpha_type;

			if constexpr (N_lim < 0) {
				return static_method<N_lim>(t_(0));
			}
			else {
				T_sigma &u0 = t_[0];
				T_sigma  un = u0;
				un ^= un << 1;
				un &= T_fix::sign.mask;
				u0 ^= un;
				U_sigma vn{un};
				vn <<= U_fix::full.depth - T_fix::full.depth;
				vn  |= U_fix::unit.mask;
				return superprocess::template static_method<N_lim>(t_(0))*
					operative_f<[] XTAL_0FN_(alias) (_xtd::bit_cast<U_alpha>)>(vn);
			}
		}

	};
};
template <int M_ism> requires in_n<M_ism,-1,-2>
struct unity<M_ism> : unity<>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		static_method(complex_field_q auto &&o)
		noexcept -> decltype(auto)
		{
			using U       = XTAL_ALL_(o);
			using U_fix    = bond::fixture<U>;

			auto constexpr N_lim_tan = N_lim;
			auto constexpr N_lim_log = 2;

			auto constexpr _1 = U_fix::haplo_0, _2pi = _1/U_fix::patio_2;
			auto constexpr _2 = U_fix::haplo_1, _4pi = _2/U_fix::patio_2;

			auto const &[x_re, x_im] = destruct_f(o);
			auto const   y_re = tangy_t<-1, 1>::template static_method<N_lim_tan>(x_im, x_re);
			auto const   w_im = square_f(x_re, x_im);
			//\
			auto const   y_im = log(w_im);
			auto const   y_im = taylor::logarithm_t< 1, 1>::template static_method<N_lim_log>(w_im);
			return complexion_f(y_re*_2, y_im*-_4pi);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
