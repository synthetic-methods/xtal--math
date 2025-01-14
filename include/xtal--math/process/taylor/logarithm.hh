#pragma once
#include "./any.hh"

#include "./sine.hh"
#include "./monologarithm.hh"
#include "../pade/tangy.hh"


XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_car=0> requires in_n<M_ism, 1,-1> and in_n<M_car, 0, 1>
struct   logarithm;

template <auto ...Ms>
using    logarithm_t = process::confined_t<logarithm<Ms...>>;

template <auto ...Ms>
XTAL_DEF_(short)
XTAL_LET logarithm_f(auto &&o, constant_q auto ...oo)
noexcept -> decltype(auto)
{
	return logarithm_t<Ms...>::template function<oo...>(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
///\
Defines `function` as the logarithm `Log[#]`, approximated by `(# - 1)/Sqrt[#]`. \

template <>
struct logarithm< 1, 0>
{
	using superprocess = process::confined_t<dilated<2>, taylor::sine<-2>>;
	
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)                 {return dysfunction<N_lim>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)                  {return dysfunction<   ~0>(XTAL_REF_(o));}
#if XTAL_SYS_(builtin)
			XTAL_0IF (real_variable_q<decltype(o)>) {return      __builtin_log(XTAL_REF_(o));}
#endif
			XTAL_0IF_(else)                       {return                log(XTAL_REF_(o));}
		}

	protected:
		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET dysfunction(auto o)
		noexcept -> decltype(auto)
		{
			return superprocess::template function<N_lim>(roots_f<2>(XTAL_MOV_(o)).template sum<-1>());
		}

	};
};
///\
Defines `function` as the antilogarithm `Exp[#]`, \
approximated by `(Sqrt[(#/2)^2 + 1] + (#/2))*# + 1`. \

template <>
struct logarithm<-1, 0>
{
	using superprocess = process::confined_t<dilated<2>, taylor::sine<+2>>;
	
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)                 {return dysfunction<N_lim>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)                  {return dysfunction<   ~0>(XTAL_REF_(o));}
#if XTAL_SYS_(builtin)
			XTAL_0IF (real_variable_q<decltype(o)>) {return      __builtin_exp(XTAL_REF_(o));}
#endif
			XTAL_0IF_(else)                       {return                exp(XTAL_REF_(o));}
		}

	protected:
		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET dysfunction(auto &&o)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;

			if constexpr (0 == N_lim) {
				auto u = o*_op::haplo_1; u += root_f<2>(term_f(one, u, u));
				return term_f(one, XTAL_MOV_(u), XTAL_REF_(o));
			}
			else {
				/**/
				auto constexpr N = below_m<0x10, (unsigned) N_lim> << 2;
				return square_f<N>(function<0>(XTAL_REF_(o)*_op::haplo_f(N)));
				/*/
				return monologarithm_t<-1>::template function<N_lim>(XTAL_REF_(o)) + one;
				/***/
			}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
///\
Defines the argument-restricted approximation of the logarithm `Log[#]`. \

template <>
struct logarithm< 1, 1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim) {return dysfunction<N_lim>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)  {return dysfunction<   ~0>(XTAL_REF_(o));}
			XTAL_0IF_(else)       {return                log(XTAL_REF_(o));}
		}

	protected:
		template <int N_lim=0>
		XTAL_DEF_(long,static)
		XTAL_LET dysfunction(real_variable_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			using _op = bond::operate<decltype(o)>;
			using U_alpha = typename _op::alpha_type;
			using U_sigma = typename _op::sigma_type;
			using U_delta = typename _op::delta_type;

		//	Log[m 2^x]
		//	Log[m] + Log[2^x]
		//	Log[m] + Log[2]*x
		//	Log[m/2^(1/2)] + Log[2]*x + Log[2^(1/2)]
		//	Log[m/2^(1/2)] + (x + 1/2)*Log[2]
		//	Log[m/2^(1/2)] + (x*2 + 1)*Log[2]/2

			U_sigma m = _xtd::bit_cast<U_sigma>(o);
			U_delta n = m - _op::unit.mask;
			m  &= _op::fraction.mask;
			m  |= _op::unit.mask;
			n >>= _op::unit.shift - one;
			n  |= one;

			U_alpha constexpr w1 =                       root_f<-2>(2.) ;
			U_alpha constexpr u1 =           logarithm_f(root_f< 2>(2.));
			auto const w    = w1 *  _xtd::bit_cast<U_alpha>(XTAL_MOV_(m));
			auto const u    = u1 *     static_cast<U_alpha>(XTAL_MOV_(n));

			return logarithm_t<1>::template function<N_lim>(XTAL_MOV_(w)) + XTAL_MOV_(u);
		}
		template <int N_lim=0>
		XTAL_DEF_(long,static)
		XTAL_LET dysfunction(complex_variable_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			using _op = bond::operate<decltype(o)>;

			auto constexpr up = one/_op::patio_1;
			auto constexpr dn =     _op::patio_1;

			auto const [u_re, u_im] = destruct_f(XTAL_REF_(o));
			auto const w_re = square_f(u_re);
			auto const w_im = square_f(u_im);

			auto const y_re = _op::haplo_1*dysfunction<N_lim>(w_re + w_im);
			auto const y_im = pade::tangy_t<-1, 1>::template function<N_lim>(u_im, u_re)*_op::patio_1;
			return {y_re, y_im};
		}

	};
};
///\
Defines argument-restricted approximation of the antilogarithm `Exp[#]`. \

template <>
struct logarithm<-1, 1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	//	TODO: Define `complex` variant!

		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim) {return dysfunction<N_lim>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)  {return dysfunction<   ~0>(XTAL_REF_(o));}
			XTAL_0IF_(else)       {return                exp(XTAL_REF_(o));}
		}

	protected:
		template <int N_lim=0>
		XTAL_DEF_(short,static)
		XTAL_LET dysfunction(real_variable_q auto o)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;
			using U_alpha = typename _op::alpha_type;

			U_alpha constexpr N_log2 = _std::numbers::ln2_v<U_alpha>;
			o *= one/N_log2; auto const n = round(o); o -= n;
			o *=     N_log2;
			return ldexp(logarithm_t<-1>::template function<N_lim>(XTAL_MOV_(o)), XTAL_MOV_(n));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
