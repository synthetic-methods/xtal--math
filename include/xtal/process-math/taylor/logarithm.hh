#pragma once
#include "./any.hh"

#include "./sine.hh"
#include "./monologarithm.hh"
#include "../gudermannian/sigmoid.hh"


XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_pow=1, int M_car=0>
	requires in_n<M_ism, 1,-1> and in_n<M_pow, 1,-1> and in_n<M_car, 0, 1>
XTAL_TYP logarithm;

template <int ...Ms>
XTAL_USE logarithm_t = process::confined_t<logarithm<Ms...>>;

template <int ...Ms>
XTAL_DEF_(return,inline)
XTAL_LET logarithm_f(auto &&o, nominal_q auto ...oo)
XTAL_0EX -> decltype(auto)
{
	return logarithm_t<Ms...>::template function<oo...>(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
/*/
template <int M_ism, int M_car>
struct logarithm<M_ism,-1, M_car>
:	process::lift<root<-1>, logarithm<M_ism, 1, M_car>>
{
};
/*/
template <int M_car>
struct logarithm< 1,-1, M_car>
:	process::lift<root<-1>, logarithm< 1, 1, M_car>>
{
};
template <int M_car>
struct logarithm<-1,-1, M_car>
{
	using superprocess = logarithm_t<-1, 1, M_car>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	//	TODO: Define `complex` variant!

		template <int N_lim=0>
		XTAL_DEF_(return)
		XTAL_SET function(auto &&o)
		XTAL_0EX -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (N_lim <  0) {
				return exp(-XTAL_REF_(o));
			}
			XTAL_0IF (0 <= N_lim) {
				return root_f<-1>(superprocess::template function<N_lim>(XTAL_REF_(o)));
			}

		}

	};
};

/***/
////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` as the logarithm `Log[#]`, \
approximated by `(# - 1)/Sqrt[#]`. \

template <>
struct logarithm< 1, 1, 0>
{
	using superprocess = process::confined_t<dilated<1>, taylor::sine<-2>>;
	
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	//	TODO: Define `complex` variant?

		template <int N_lim=0>
		XTAL_DEF_(return,inline)
		XTAL_SET function(auto &&o)
		XTAL_0EX -> decltype(auto)
		{
			using _std::log;

			using _op = bond::operate<decltype(o)>;
			auto constexpr _1 = _op::alpha_1;

			XTAL_IF0
			XTAL_0IF (N_lim <  0) {
				return log(XTAL_REF_(o));
			}
			XTAL_0IF (0 <= N_lim) {
				return superprocess::template function<N_lim>(
					(o - _1)*root_f<-2>(XTAL_REF_(o))
				);
			}
		}

	};
};
///\
Defines `function` as the antilogarithm `Exp[#]`, \
approximated by `(Sqrt[(#/2)^2 + 1] + (#/2))*# + 1`. \

template <>
struct logarithm<-1, 1, 0>
{
	using superprocess = process::confined_t<dilated<1>, taylor::sine<+2>>;
	
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline)
		XTAL_SET function(auto &&o)
		XTAL_0EX -> decltype(auto)
		{
			using _std::exp;

			using _op = bond::operate<decltype(o)>;
			auto constexpr _1 =   _op::alpha_1;
			auto const      u = o*_op::haplo_1;

			XTAL_IF0
			XTAL_0IF (N_lim <  0) {
				return exp(XTAL_REF_(o));
			}
			XTAL_0IF (0 == N_lim) {
				return horner::term_f(_1, XTAL_REF_(o), u + root_f<2>(horner::term_f(_1, u, u)));
			}
			XTAL_0IF (1 <= N_lim) {
				return monologarithm_t<-1>::template function<N_lim>(XTAL_REF_(o)) + _1;
			}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
///\
Defines argument-reduced approximation of the logarithm `Log[#]`. \

template <>
struct logarithm< 1, 1, 1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	//	TODO: Define `complex` variant!

		template <int N_lim=0>
		XTAL_DEF_(return,inline)
		XTAL_SET function(complex_number_q auto u)
		XTAL_0EX -> decltype(auto)
		{
			using _op = bond::operate<decltype(u)>;

			auto constexpr up = _op::alpha_1/_op::patio_1;
			auto constexpr dn =              _op::patio_1;

			auto const [u_re, u_im] = invalued_f(XTAL_REF_(u));
			auto const w_re = square_f(u_re);
			auto const w_im = square_f(u_im);

			auto const y_re = _op::haplo_1*function<N_lim>(w_re + w_im);
			auto const y_im = gudermannian::sigmoid_t<-1>::template function<-!!N_lim>(up*u_im/u_re)*dn;
			return complexion_f(y_re, y_im);
		}
		template <int N_lim=0>
		XTAL_DEF_(return,inline)
		XTAL_SET function(real_number_q auto o)
		XTAL_0EX -> decltype(auto)
		{
			using _std::log;

			using _op = bond::operate<decltype(o)>;
			using U_alpha = typename _op::alpha_type;
			using U_sigma = typename _op::sigma_type;
			using U_delta = typename _op::delta_type;

			XTAL_IF0
			XTAL_0IF (N_lim <  0) {
				return log(o);
			}
			XTAL_0IF (0 <= N_lim) {
			//	Log[m 2^x]
			//	Log[m] + Log[2^x]
			//	Log[m] + Log[2]*x
			//	Log[m/Sqrt[2]] + Log[2]*(x + 1/2)
			//	Log[m/Sqrt[2]] + Log[2]*(x*2 + 1)/2

				U_alpha constexpr          _1 = 1;
				U_alpha constexpr N_sqrt_half = 0.7071067811865475244008443621048490393e+0L;
				U_alpha constexpr N_half_log2 = 0.3465735902799726547086160607290882840e+0L;

				U_sigma m = _xtd::bit_cast<U_sigma>(o);
				U_delta n = (m << 1) - (_op::unit.mask << 1);
				m     &= _op::fraction.mask;
				m     |= _op::unit.mask;
				n    >>= _op::unit.shift;
				n     |= 1;
				auto w = N_sqrt_half*_xtd::bit_cast<U_alpha>(m);
				auto u = N_half_log2*   static_cast<U_alpha>(n);
				return logarithm_t<1, 1, 0>::template function<N_lim>(w) + u;
			}
		}

	};
};
///\
Defines argument-reduced approximation of the antilogarithm `Exp[#]`. \

template <>
struct logarithm<-1, 1, 1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	//	TODO: Define `complex` variant!

		template <int N_lim=0>
		XTAL_DEF_(return,inline)
		XTAL_SET function(real_number_q auto o)
		XTAL_0EX -> decltype(auto)
		{
			using _std::log;

			using _op = bond::operate<decltype(o)>;
			using U_alpha = typename _op::alpha_type;

			XTAL_IF0
			XTAL_0IF (N_lim <  0) {
				return exp(o);
			}
			XTAL_0IF (0 <= N_lim) {
				U_alpha constexpr    _1 = 1;
				U_alpha constexpr M_ln2 = 1.4426950408889634073599246810018921374e+0L;// 1/Log[2]
				U_alpha constexpr N_ln2 = 0.6931471805599453094172321214581765681e+0L;// 1*Log[2]
				o *= M_ln2;
				auto const n = _std::round(o);
				o -= n;
				o *= N_ln2;
				return _std::ldexp(logarithm_t<-1>::template function<N_lim>(o), n);
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
