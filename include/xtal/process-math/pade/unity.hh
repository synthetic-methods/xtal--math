#pragma once
#include "./any.hh"

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
XTAL_TYP unity<>
{
	using limit_type = occur::math::limit_t<(1<<3)>;

	template <class S>
	using subtype = bond::compose_s<S, resource::invoice<void
	,	typename limit_type::template dispatch<>
	>>;

};
template <int M_ism=1, typename ...As>
XTAL_USE unity_t = process::confined_t<unity<M_ism, As...>, unity<>>;


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires in_n<M_ism, 0, 1, 2>
struct unity<M_ism> : unity<>
{
	using superprocess = process::lift_t<square<M_ism, 0>, _detail::subunity<M_ism,-0>>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1, class U>
		XTAL_DEF_(return,inline)
		XTAL_SET function(_std::initializer_list<U> o)
		XTAL_0EX -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			_std::complex<U> w; auto &m = involved_f(w);
			_std::copy_n(point_f(o), 2, m);
			return function<N_lim>(w);
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline)
		XTAL_SET function(complex_field_q auto const &t)
		XTAL_0EX -> decltype(auto)
		{
			return function<N_lim>(t.real(), t.imag());
		}
		template <int N_lim=-1>
		XTAL_DEF_(return,inline)
		XTAL_SET function(auto &&t_re, simplex_field_q auto &&t_im)
		XTAL_0EX -> decltype(auto)
		{
			using _op = bond::operate<decltype(t_re), decltype(t_im)>;
			return function<N_lim>(XTAL_REF_(t_re))*exp(XTAL_REF_(t_im)*_op::patio_f(-2));
		}

		template <int N_lim=-1>
		XTAL_DEF_(return)
		XTAL_SET function(simplex_field_q auto o)
		XTAL_0EX -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			using Op_alpha = typename _op::alpha_type;

			if constexpr (N_lim < 0) {
				using namespace _std;

				auto w = o*_op::patio_2;
				XTAL_IF0
				XTAL_0IF (1 == M_ism) {return complexion_f(cos (w), sin (w));}
				XTAL_0IF (2 == M_ism) {return complexion_f(cosh(w), sinh(w));}
			}
			else {
				auto constexpr assigned_f = [] (int i) XTAL_0FN -> Op_alpha {return (i << 1) - 1;};
				auto w = wrap_f(o);
				auto m = objective_f(wrap_f(w*_op::diplo_1)*_op::haplo_1);
				return superprocess::template function<N_lim>(m)*operative_f(assigned_f, m == w);
			}
		}
		template <int N_lim=-1>
		XTAL_DEF_(return)
		XTAL_SET function(algebra::bicycle_q auto t_)
		XTAL_0EX -> decltype(auto)
		{
			using T_ = XTAL_ALL_(t_);
			using T_op = bond::operate<decltype(t_[0])>;// Underlying...
			using U_op = bond::operate<decltype(t_(0))>;
			using T_sigma = typename T_op::sigma_type;
			using U_sigma = typename U_op::sigma_type;
			using U_alpha = typename U_op::alpha_type;
			XTAL_USE alpha_f = decltype([] XTAL_1FN_(_xtd::bit_cast<U_alpha>));

			if constexpr (N_lim < 0) {
				return function<N_lim>(t_(0));
			}
			else {
				T_sigma &u0 = t_[0];
				T_sigma  un = u0;
				un ^= un << 1;
				un &= T_op::sign.mask;
				u0 ^= un;
				U_sigma vn{un};
				vn <<= U_op::full.depth - T_op::full.depth;
				vn  |= U_op::unit.mask;
				return superprocess::template function<N_lim>(t_(0))*inoperative_f<alpha_f>(vn);
			}
		}

	};
};
template <int M_ism> requires in_n<M_ism,-1,-2>
struct unity<M_ism> : unity<>
{
	using superprocess = process::lift_t<square<M_ism, 0>, _detail::subunity<M_ism,-0>>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline)
		XTAL_SET function(complex_field_q auto &&o)
		XTAL_0EX -> decltype(auto)
		{
			using horner::term_f;
			using horner::square_f;
			using _op = bond::operate<decltype(o)>;

			XTAL_LET N_lim_tan = N_lim;
			XTAL_LET N_lim_log = 2;

			XTAL_LET _1 = _op::haplo_0, _2pi = _1/_op::patio_2;
			XTAL_LET _2 = _op::haplo_1, _4pi = _2/_op::patio_2;

			auto const &[x_re, x_im] = involved_f(o);
			auto const   y_re = tangy_t<-1, 1>::template function<N_lim_tan>(x_im, x_re);
			auto const   w_im = square_f(x_re, x_im);
			//\
			auto const   y_im = log(w_im);
			auto const   y_im = taylor::logarithm_t< 1, 1, 1>::template function<N_lim_log>(w_im);
			return complexion_f(y_re*_2, y_im*-_4pi);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
