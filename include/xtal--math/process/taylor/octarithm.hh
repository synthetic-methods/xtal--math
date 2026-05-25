#pragma once
#include "./any.hh"

#include "./logarithm.hh"




XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Produces the base-two logarithm and its inverse.
\todo    Either define `subdivision` to apply `1/N_ind` prescaling, or supply an `attach`ed parameter to do so.
\todo    Allow for recentering w.r.t. `A4`.
*/
template <int M_ism=2, int M_div=1>
struct  octarithm;


////////////////////////////////////////////////////////////////////////////////

template <int M_ism, int M_div> requires in_v<M_ism,  1,  2>
struct octarithm<M_ism, M_div>
{
	using superkind = logarithm<1, 1>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

	//	TODO: Define `complex` variant!

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(auto &&...oo)
		noexcept -> auto
		requires (2 <= sizeof...(oo))
		{
			return method<Ns...>((XTAL_REF_(oo) *...* one));
		}

		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		method(auto o)
		noexcept -> XTAL_ALL_(o)
		requires un_v<atom::quantify_q<decltype(o)>>
		{
			using U_alpha = typename bond::fit<decltype(o)>::alpha_type;
			U_alpha constexpr N_ =  one*U_alpha{M_div};
			U_alpha constexpr N2 =  one*std::numbers::ln2_v<U_alpha>;
			U_alpha constexpr N1 = -two*std::numbers:: pi_v<U_alpha>;

			o = S_::template method<N_lim>(XTAL_MOV_(o));
			if constexpr (M_div  >  1) {o *=     N_;}
			if constexpr (M_ism ==  1) {o *= one/N1;}
			if constexpr (M_ism ==  2) {o *= one/N2;}
			return o;
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

	};
};
template <int M_ism, int M_div> requires in_v<M_ism, -1, -2>
struct octarithm<M_ism, M_div>
{
//	using superkind = logarithm<-1, 1>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	//	TODO: Define `complex` variant!

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method(auto &&...oo)
		noexcept -> auto
		requires (2 <= sizeof...(oo))
		{
			return method<Ns...>((XTAL_REF_(oo) *...* one));
		}

		template <int N_lim=0>
		XTAL_DEF_(return,set)
		method(integral_q auto &&x)
		noexcept -> auto
		{
			using U_fit   = bond::fit<decltype(x)>;
			using U_alpha = typename U_fit::alpha_type;

			static_assert(M_ism == -2);// For now...

			auto constexpr K0 = U_fit::alpha_f(1);
			auto constexpr KN = U_fit::alpha_f(2);
			auto constexpr K1 = root_f<M_div>(KN);
			auto constexpr K_ = [&]<auto ...I> (bond::seek_in_t<I...>)
				XTAL_0FN_(to) (std::array{monomial_f<I>(K1)...}) (bond::seek_to_t<M_div>{});

			XTAL_IF1_(consteval) {
				return U_fit::diplo_f(U_fit::ratio_f(XTAL_REF_(x), M_div));
			}
			XTAL_0IF (1 == M_div) {
				return xtd::ldexp(KN, XTAL_REF_(x));
			}
			XTAL_0IF (2 <= M_div) {
				auto const r = modulo_f<(unsigned) M_div>(x);
				auto const n = x - r;
				return xtd::ldexp(K_[r], n/M_div);
			}
		}
		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		method(real_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			using U_alpha = typename bond::fit<decltype(o)>::alpha_type;
			U_alpha constexpr _N =  one/U_alpha{M_div};
			U_alpha constexpr N2 =  one*std::numbers::ln2_v<U_alpha>;
			U_alpha constexpr N1 = -two*std::numbers:: pi_v<U_alpha>;

			if constexpr (M_div  >  1) {o *=    _N;}
			if constexpr (M_ism == -1) {o *= N1/N2;}
		//	if constexpr (M_ism == -2) {o *= N2/N2;}
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)          {return  methox<N_lim>(XTAL_MOV_(o));}
			XTAL_0IF_(consteval)           {return  methox<   ~0>(XTAL_MOV_(o));}
#if XTAL_SYS_(builtin)
			XTAL_0IF (real_q<decltype(o)>) {return __builtin_exp2(XTAL_MOV_(o));}
#endif
			XTAL_0IF_(else)                {return         exp(N2*XTAL_MOV_(o));}
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
		XTAL_DEF_(return,set)
		methox(real_variable_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			using U_fit   = bond::fit<decltype(o)>;
			using U_alpha = typename U_fit::alpha_type;
			using U_delta = typename U_fit::delta_type;
			static_assert(-1 <= N_lim);
			auto constexpr L4 = term_f(4, 4, N_lim&0b11);
			auto constexpr U1 = U_fit::haplo_f(L4)*std::numbers::ln2_v<U_alpha>;
			auto constexpr W2 = U_fit::haplo_f(1)*U1*U1;
			auto const n = nearest_f<>(o); o -= n;
			o = term_f(one, o, term_f(U1, W2, o));
			#pragma unroll
			for (int i{}; i < L4; ++i) {o *= o;}
			return xtd::ldexp(o, static_cast<U_delta>(n));
		}

	};
};

////////////////////////////////////////////////////////////////////////////////

template <int M_ism=2, int M_div=1>
XTAL_TYP_(let) octarithm_t = process::confined_t<
	octarithm  <M_ism, M_div>
>;
template <int M_ism=2, int M_div=1, int N_lim=2>
XTAL_DEF_(let) octarithm_f = [] XTAL_1FN_(call) (
	octarithm_t<M_ism, M_div>::template method<N_lim>
);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
