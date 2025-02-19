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
struct  unity<>
{
	using limit_type = occur::math::limit_t<(1<<3)>;

	template <class S>
	using subtype = bond::compose_s<S, provision::voiced<void
	,	typename limit_type::template dispatch<>
	>>;

};
template <int M_ism=1, typename ...As>
using   unity_t = process::confined_t<unity<M_ism, As...>, unity<>>;


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
			auto constexpr exp = [] XTAL_1FN_(call) (taylor::logarithm_t<-1, 1>::template method_f<2>);
			using U_fit = bond::fit<decltype(t_re), decltype(t_im)>;
			return method_f<N_lim>(XTAL_REF_(t_re))*exp(XTAL_REF_(t_im)*U_fit::patio_f(-2));
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


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
