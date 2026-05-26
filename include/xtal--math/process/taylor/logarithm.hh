#pragma once
#include "./any.hh"

#include "./sine.hh"
#include "./tangent.hh"
#include "../pade/tangy.hh"
#include "../pade/arc.hh"

XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_car=0>
XTAL_TYP_(new) logarithm;

template <int M_ism=1, int M_car=0>
XTAL_TYP_(let) logarithm_t = process::confined_t<
	logarithm  <M_ism, M_car>
>;
template <int M_ism=1, int M_car=0, int ...Ns>
XTAL_DEF_(let) logarithm_f = [] XTAL_1FN_(call) (
	logarithm_t<M_ism, M_car>{}.template method<Ns...>
);


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines the logarithm `Log`, approximated by `((# - 1)/Sqrt[#] &)`.
*/
template <>
struct logarithm< 1, 0>
{
	XTAL_DEF_(set) sinh_aexp_ = [] (auto &&u)
		XTAL_0FN_(to) (roots_f<2>(XTAL_REF_(u)).template sum<-1>());
	
	using superprocess = process::confined_t<typename dilate_t<2>::template infix<>, taylor::sine<-2>>;
	
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		method(auto &&o)
		noexcept -> decltype(auto)
		requires un_v<atom::quantify_q<decltype(o)>>
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)                   {return methox<N_lim>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)                    {return methox<   ~0>(XTAL_REF_(o));}
	#if   XTAL_SYS_(builtin)
			XTAL_0IF (real_variable_q<decltype(o)>) {return __builtin_log(XTAL_REF_(o));}
	#endif
			XTAL_0IF_(else)                         {return           log(XTAL_REF_(o));}
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(auto &&o)
		noexcept -> decltype(auto)
		requires in_v<atom::quantify_q<decltype(o)>>
		{
			return XTAL_ALL_(o)::template zip_from<[]
				XTAL_1FN_(call) (method<Ns...>)>(XTAL_REF_(o));
		}

	protected:
		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		methox(auto &&o)
		noexcept -> decltype(auto)
		{
			return superprocess{}.template method<N_lim>(sinh_aexp_(XTAL_REF_(o)));
		}

	};
};
/*!
\brief   Defines the antilogarithm `Exp[#]`,

Approximated by `((Sqrt[(#/2)^2 + 1] + (#/2))^2 &)`.
*/
template <>
struct logarithm<-1, 0>
{
	XTAL_DEF_(set) exp_asinh_ = [] (auto &&u)
		XTAL_0FN_(to) (u + root_f<2>(term_f<1, 2>(one, u)));
	
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		method(auto &&o)
		noexcept -> decltype(auto)
		requires un_v<atom::quantify_q<decltype(o)>>
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)          {return methox<N_lim>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)           {return methox<   ~0>(XTAL_REF_(o));}
	#if	XTAL_SYS_(builtin)
			XTAL_0IF (real_q<decltype(o)>) {return __builtin_exp(XTAL_REF_(o));}
	#endif
			XTAL_0IF_(else)                {return           exp(XTAL_REF_(o));}
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(auto &&o)
		noexcept -> decltype(auto)
		requires in_v<atom::quantify_q<decltype(o)>>
		{
			return XTAL_ALL_(o)::template zip_from<[]
				XTAL_1FN_(call) (method<Ns...>)>(XTAL_REF_(o));
		}

//	protected:
		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		methox(auto &&o)
		noexcept -> decltype(auto)
		{
			return square_f(exp_asinh_(taylor::sine_t<2, 0>{}.
				template method<N_lim>(half*XTAL_REF_(o))));// Actual inverse...
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
		method(auto &&o)
		noexcept -> XTAL_ALL_(o)
		requires un_v<atom::quantify_q<decltype(o)>>
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim) {return methox<N_lim>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)  {return methox<   ~0>(XTAL_REF_(o));}
			XTAL_0IF_(else)       {return           log(XTAL_REF_(o));}
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(auto &&o)
		noexcept -> decltype(auto)
		requires in_v<atom::quantify_q<decltype(o)>>
		{
			return XTAL_ALL_(o)::template zip_from<[]
				XTAL_1FN_(call) (method<Ns...>)>(XTAL_REF_(o));
		}

	protected:
		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		methox(real_variable_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			auto constexpr K_pi  =  bond::fit<decltype(o)>::patio_2;
			return -K_pi*pade::arc_t<-1, 1>{}.template method<N_lim>(XTAL_MOV_(o));
		}
		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		methox(complex_variable_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			auto constexpr K_pi  =  bond::fit<decltype(o)>::patio_1;
			auto const    [u_re,
			               u_im] =  destruct_f(XTAL_REF_(o));
			auto const     y_re  = -K_pi*pade::arc_t<~0, 1>{}.template method<N_lim>(u_im, u_re);
			auto const     y_im  =  K_pi*pade::arc_t< 0, 1>{}.template method<N_lim>(u_im, u_re);
			return  {y_re, y_im};
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
		method(auto &&o)
		noexcept -> XTAL_ALL_(o)
		requires un_v<atom::quantify_q<decltype(o)>>
		{
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)          {return methox<N_lim>(XTAL_REF_(o));}
			XTAL_0IF_(consteval)           {return methox<   ~0>(XTAL_REF_(o));}
	#if	XTAL_SYS_(builtin)
			XTAL_0IF (real_q<decltype(o)>) {return __builtin_exp(XTAL_REF_(o));}
	#endif
			XTAL_0IF_(else)                {return           exp(XTAL_REF_(o));}
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(auto &&o)
		noexcept -> decltype(auto)
		requires in_v<atom::quantify_q<decltype(o)>>
		{
			return XTAL_ALL_(o)::template zip_from<[]
				XTAL_1FN_(call) (method<Ns...>)>(XTAL_REF_(o));
		}

	protected:
		template <int N_lim=0> requires (0 == N_lim)
		XTAL_DEF_(return,inline,set)
		methox(real_variable_q auto x)
		noexcept -> XTAL_ALL_(x)
		{
			using U_fit   = bond::fit<decltype(x)>;
			using U_alpha = typename U_fit::alpha_type;
			using U_delta = typename U_fit::delta_type;
			using U_sigma = typename U_fit::sigma_type;

			U_alpha constexpr M  =  one << U_fit::fraction.depth;
			U_alpha constexpr X  =  std::numbers::log2e_v<U_alpha>;       // 1 / Log@2
			U_alpha constexpr X_ =  X*M;
		//	U_alpha constexpr K  =  0.021838724451800924121345965705447L; //     Log2[1/E/Log@2 + 1/2]/2
			U_alpha constexpr K_ =  0.978161275548199075878654034294553L; // 1 - Log2[1/E/Log@2 + 1/2]/2
			U_delta constexpr I  =  std::bit_cast<U_delta>(K_);

			auto i = I + static_cast<U_delta>(X_*x); x = std::bit_cast<U_alpha>(i);
			i &= U_fit::fraction.mask;
			i |= U_fit::    unit.mask;
			auto w = std::bit_cast<U_alpha>(i);
			U_alpha constexpr a0 = -1.15254531398409900;
			U_alpha constexpr a1 =  1.25646677591560070;
			U_alpha constexpr a2 = -0.41619451269060210;
			U_alpha constexpr a3 =  0.06863680575030143;
			U_alpha constexpr a_ =  0.97104234634093000;
			x *= a_ + process::math::square_f(process::math::termial_f(w, a0, a1, a2, a3));
			return x;
		}
		template <int N_lim=0> requires (0 != N_lim)
		XTAL_DEF_(return,set)
		methox(real_variable_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			using U_fit = bond::fit<decltype(o)>;
			using U_alpha = typename U_fit::alpha_type;
			using U_sigma = typename U_fit::sigma_type;
			using U_delta = typename U_fit::delta_type;

			U_alpha constexpr  N_log2 = one*std::numbers::ln2_v<U_alpha>;
			U_alpha constexpr _N_log2 = one/std::numbers::ln2_v<U_alpha>;
			U_alpha constexpr _N_dns2 = U_fit::dnsilon_f(1)*half;
			U_alpha n;
			U_delta N;
			o *= _N_log2;
			XTAL_IF1_(consteval) {
				N = static_cast<U_delta>(o + part_f<signed>(o)*_N_dns2);
				n = static_cast<U_alpha>(N);
			}
			XTAL_0IF_(else) {
				n = round(o);
				N = static_cast<U_delta>(n);
			}
			o -= n;
			o *= N_log2;
			o  = logarithm_t<-1, 0>::template method<N_lim>(XTAL_MOV_(o));
			XTAL_IF1_(consteval) {
				auto  m = xtd::bit_cast<U_sigma>(o);
				N <<= U_fit::exponent.shift;
				m  += N;
				return xtd::bit_cast<U_alpha>(m);
			}
			XTAL_0IF_(else) {
				return xtd::ldexp(o, N);

			}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
