#pragma once
#include "./any.hh"

#include "./sine.hh"
#include "./tangent.hh"
#include "./monologarithm.hh"
#include "../pade/tangy.hh"

XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_car=0>
struct  logarithm;

template <int M_ism=1, int M_car=0>
using   logarithm_t = process::confined_t<logarithm<M_ism, M_car>>;

template <int M_ism=1, int M_car=0, int ...Ns>
XTAL_DEF_(let)
logarithm_f = [] XTAL_1FN_(call) (logarithm_t<M_ism, M_car>::template method_f<Ns...>);


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines the logarithm `Log`, approximated by `((# - 1)/Sqrt[#] &)`.
*/
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
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)                   {return approximate<N_lim>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)                    {return approximate<   ~0>(XTAL_REF_(o));}
#if XTAL_SYS_(builtin)
			XTAL_0IF (real_variable_q<decltype(o)>) {return      __builtin_log(XTAL_REF_(o));}
#endif
			XTAL_0IF_(else)                         {return                log(XTAL_REF_(o));}
		}

	protected:
		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		approximate(auto o)
		noexcept -> decltype(auto)
		{
			return superprocess::template method_f<N_lim>(roots_f<2>(XTAL_MOV_(o)).template sum<-1>());
		}

	};
};
/*!
\brief   Defines the antilogarithm `Exp[#]`,

Approximated by `((Sqrt[(#/2)^2 + 1] + (#/2))*# + 1 &)`.
*/
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
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)          {return approximate<N_lim>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)           {return approximate<   ~0>(XTAL_REF_(o));}
#if XTAL_SYS_(builtin)
			XTAL_0IF (real_q<decltype(o)>) {return      __builtin_exp(XTAL_REF_(o));}
#endif
			XTAL_0IF_(else)                {return                exp(XTAL_REF_(o));}
		}

	protected:
		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		approximate(auto &&o)
		noexcept -> decltype(auto)
		{
			using U_fit = bond::fit<decltype(o)>;
			/*/
			auto constexpr zoom_up2 = U_fit::haplo_f(N_lim*2 + 1);

			auto v = aspect_f<signed>(o), w = square_f(XTAL_REF_(o));
			w = term_f(one, zoom_up2, XTAL_MOV_(w));

			#pragma unroll
			for (int i{}; i < N_lim; ++i) {
				w = term_f(-one, two, square_f(XTAL_MOV_(w)));
			}
			return root_f<2>(term_f(-one, two, w, term_f(w, v, root_f<2>(term_f<1, 2>(-one, w)))));
			/*/
			if constexpr (0 == N_lim) {
				auto u = o*U_fit::haplo_1; u += root_f<2>(term_f(one, u, u));
				return term_f(one, XTAL_MOV_(u), XTAL_REF_(o));
			}
			else {
				auto constexpr N = below_v<0x10, (unsigned) N_lim> << 2;
				return square_f<N>(method_f<0>(XTAL_REF_(o)*U_fit::haplo_f(N)));
			//	return monologarithm_t<-1>::template method_f<N_lim>(XTAL_REF_(o)) + one;
			}
			/***/
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines the argument-restricted approximation of the logarithm `Log[#]`.
*/
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
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim) {return approximate<N_lim>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)  {return approximate<   ~0>(XTAL_REF_(o));}
			XTAL_0IF_(else)       {return                log(XTAL_REF_(o));}
		}

	protected:
		template <int N_lim=0>
		XTAL_DEF_(return,set)
		approximate(real_variable_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			using _fit = bond::fit<decltype(o)>;
			using U_alpha = typename _fit::alpha_type;
			using U_sigma = typename _fit::sigma_type;
			using U_delta = typename _fit::delta_type;

		//	Log[m 2^x]
		//	Log[m] + Log[2^x]
		//	Log[m] + Log[2]*x
		//	Log[m/2^(1/2)] + Log[2]*x + Log[2^(1/2)]
		//	Log[m/2^(1/2)] + (x + 1/2)*Log[2]
		//	Log[m/2^(1/2)] + (x*2 + 1)*Log[2]/2

			U_sigma m = _xtd::bit_cast<U_sigma>(o);
			U_delta n = m - _fit::unit.mask;
			m  &= _fit::fraction.mask;
			m  |= _fit::unit.mask;
			n >>= _fit::unit.shift - one;
			n  |= one;

			U_alpha constexpr w1 =                       root_f<-2>(2.) ;
			U_alpha constexpr u1 =        logarithm_f<1>(root_f< 2>(2.));
			auto const w    = w1 *  _xtd::bit_cast<U_alpha>(XTAL_MOV_(m));
			auto const u    = u1 *     static_cast<U_alpha>(XTAL_MOV_(n));

			return logarithm_t<1>::template method_f<N_lim>(XTAL_MOV_(w)) + XTAL_MOV_(u);
		}
		template <int N_lim=0>
		XTAL_DEF_(return,set)
		approximate(complex_variable_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			using _fit = bond::fit<decltype(o)>;

			auto constexpr up = one/_fit::patio_1;
			auto constexpr dn =     _fit::patio_1;

			auto const [u_re, u_im] = destruct_f(XTAL_REF_(o));
			auto const w_re = square_f(u_re);
			auto const w_im = square_f(u_im);

			auto const y_re = _fit::haplo_1*approximate<N_lim>(w_re + w_im);
			auto const y_im = pade::tangy_t<-1, 1>::template method_f<N_lim>(u_im, u_re)*_fit::patio_1;
			return {y_re, y_im};
		}

	};
};
/*!
\brief   Defines argument-restricted approximation of the antilogarithm `Exp[#]`.
*/
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
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim) {return approximate<N_lim>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)  {return approximate<   ~0>(XTAL_REF_(o));}
			XTAL_0IF_(else)       {return                exp(XTAL_REF_(o));}
		}

	protected:
		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		approximate(real_variable_q auto o)
		noexcept -> decltype(auto)
		{
			using _fit = bond::fit<decltype(o)>;
			using U_alpha = typename _fit::alpha_type;

			U_alpha constexpr N_log2 = _std::numbers::ln2_v<U_alpha>;
			o *= one/N_log2; auto const n = round(o); o -= n;
			o *=     N_log2;
			return ldexp(logarithm_t<-1>::template method_f<N_lim>(XTAL_MOV_(o)), XTAL_MOV_(n));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
