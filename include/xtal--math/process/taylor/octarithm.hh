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
		requires un_v<atom::groupoid_q<decltype(o)>>
		{
			using U_alpha = typename bond::fit<decltype(o)>::alpha_type;
			U_alpha constexpr N_ =  one*U_alpha{M_div};
			U_alpha constexpr N2 =  one*_std::numbers::ln2_v<U_alpha>;
			U_alpha constexpr N1 = -two*_std::numbers:: pi_v<U_alpha>;

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
		requires in_v<atom::groupoid_q<decltype(o)>>
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
				XTAL_0FN_(to) (_std::array{monomial_f<I>(K1)...}) (bond::seek_to_t<M_div>{});

			XTAL_IF1_(consteval) {
				return U_fit::diplo_f(U_fit::ratio_f(XTAL_REF_(x), M_div));
			}
			XTAL_0IF (1 == M_div) {
				return _std::ldexp(KN, XTAL_REF_(x));
			}
			XTAL_0IF (2 <= M_div) {
				auto const r = modulo_f<(unsigned) M_div>(x);
				auto const n = x - r;
				return _std::ldexp(K_[r], n/M_div);
			}
		}
		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		method(real_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			using U_alpha = typename bond::fit<decltype(o)>::alpha_type;
			U_alpha constexpr _N =  one/U_alpha{M_div};
			U_alpha constexpr N2 =  one*_std::numbers::ln2_v<U_alpha>;
			U_alpha constexpr N1 = -two*_std::numbers:: pi_v<U_alpha>;

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
		requires in_v<atom::groupoid_q<decltype(o)>>
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
			using U_sigma = typename U_fit::sigma_type;
			using U_delta = typename U_fit::delta_type;
			using V_delta = int;

			XTAL_IF1_(consteval) {
				auto const d = U_fit::dnsilon_f(1)*half*part_f<signed>(o);
				auto const N = static_cast<U_delta>(o + d);
				auto const n = static_cast<U_alpha>(N);
				methoxy<N_lim>(o, n);
				auto m = _xtd::bit_cast<U_sigma>(o);
				m += N << U_fit::exponent.shift;
				return _xtd::bit_cast<U_alpha>(m);
			}
			XTAL_0IF_(else) {
				auto const n = round(o);
				auto const N = static_cast<V_delta>(n);
				methoxy<N_lim>(o, n);
				return _std::ldexp(XTAL_MOV_(o), XTAL_MOV_(N));
			}
		}
		template <int N_lim=0>
		XTAL_DEF_(inline,set)
		methoxy(real_variable_q auto &o, real_variable_q auto n)
		noexcept -> void
		{
			static_assert(-1 <= N_lim);
			using U_fit = bond::fit<decltype(o)>;
			using U_alpha = typename U_fit::alpha_type;

			auto constexpr L4 = term_f(4, 4, N_lim&0b11);
			auto constexpr U1 = U_fit::haplo_f(L4)*_std::numbers::ln2_v<U_alpha>;
			auto constexpr W2 = U_fit::haplo_f(1)*U1*U1;

			o -= n;
			o  = term_f(one, o, term_f(U1, W2, o));
			#pragma unroll
			for (int i{}; i < L4; ++i) {
				o *= o;
			}
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
