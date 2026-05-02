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
\todo    Allow for recentering w.r.t. `A4`, e.g. `440*Exp2[-69/12]*Exp2[#/12]&`,
         or `440*0.1858136117191751604527105712350021e-1L*octarithm_f<-2>(n/m_div)`.
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
		XTAL_DEF_(return,inline,let)
		method(auto &&...oo)
		const noexcept -> auto
		requires (2 <= sizeof...(oo))
		{
			return method<Ns...>((XTAL_REF_(oo) *...* one));
		}

		template <int N_lim=0>
		XTAL_DEF_(return,inline,let)
		method(auto o)
		const noexcept -> XTAL_ALL_(o)
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
		XTAL_DEF_(return,inline,let)
		method(auto &&o)
		const noexcept -> decltype(auto)
		requires in_v<atom::groupoid_q<decltype(o)>>
		{
			return XTAL_ALL_(o)::template zip_from<[]
				XTAL_1FN_(call) (subtype{}.template method<Ns...>)>(XTAL_REF_(o));
		}

	};
};
template <int M_ism, int M_div> requires in_v<M_ism, -1, -2>
struct octarithm<M_ism, M_div>
{
	using superkind = logarithm<-1, 1>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

	//	TODO: Define `complex` variant!

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto &&...oo)
		const noexcept -> auto
		requires (2 <= sizeof...(oo))
		{
			return method<Ns...>((XTAL_REF_(oo) *...* one));
		}

		template <int N_lim=0>
		XTAL_DEF_(return,inline,let)
		method(integral_q auto &&x)
		const noexcept -> auto
		{
			using X = _xtd::make_signed_t<unstruct_t<decltype(x)>>;
			using U_fit = bond::fit<X>;

			static_assert(M_ism == -2);// For now...

			XTAL_IF1_(consteval) {
				return U_fit::diplo_f(U_fit::ratio_f(1, M_div)*XTAL_REF_(x));
			}
			XTAL_0IF (1 == M_div) {
				return _std::ldexp(U_fit::alpha_f(2), XTAL_REF_(x));
			}
			XTAL_0IF (2 <= M_div) {
				auto constexpr N     = root_f<M_div>(U_fit::diplo_1);
				auto           n     = U_fit::alpha_1;
				auto const     o     = _std::div((X) XTAL_REF_(x), (X) M_div);
				auto const     r     =        o.rem;
				auto const    _r     = 0x10 - o.rem;
				if constexpr (M_div <= 0x10) {
					switch (_r) {
					case 0x0: n *= N;
					case 0x1: n *= N;
					case 0x2: n *= N;
					case 0x3: n *= N;
					case 0x4: n *= N;
					case 0x5: n *= N;
					case 0x6: n *= N;
					case 0x7: n *= N;
					case 0x8: n *= N;
					case 0x9: n *= N;
					case 0xA: n *= N;
					case 0xB: n *= N;
					case 0xC: n *= N;
					case 0xD: n *= N;
					case 0xE: n *= N;
					case 0xF: n *= N;
					}
				}
				else {
					for (X i{}; i < r; ++i) {n *= N;}
				}
				return _std::ldexp(XTAL_MOV_(n), XTAL_MOV_(o).quot);
			}
		}
		template <int N_lim=0>
		XTAL_DEF_(return,inline,let)
		method(real_q auto o)
		const noexcept -> XTAL_ALL_(o)
		{
			using U_alpha = typename bond::fit<decltype(o)>::alpha_type;
			U_alpha constexpr _N =  one/U_alpha{M_div};
			U_alpha constexpr N2 =  one*_std::numbers::ln2_v<U_alpha>;
			U_alpha constexpr N1 = -two*_std::numbers:: pi_v<U_alpha>;

			if constexpr (M_div  >  1) {o *=    _N;}
			if constexpr (M_ism == -1) {o *= N1/N2;}
		//	if constexpr (M_ism == -2) {o *= N2/N2;}
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)          {return  method_approximant<N_lim>(XTAL_MOV_(o));}
			XTAL_0IF_(consteval)           {return  method_approximant<   ~0>(XTAL_MOV_(o));}
		//	XTAL_0IF (0 <= N_lim)          {return S_::template method<N_lim>(XTAL_MOV_(o)*N2);}
		//	XTAL_0IF_(consteval)           {return S_::template method<   ~0>(XTAL_MOV_(o)*N2);}
#if XTAL_SYS_(builtin)
			XTAL_0IF (real_q<decltype(o)>) {return             __builtin_exp2(XTAL_MOV_(o));}
#endif
			XTAL_0IF_(else)                {return                     exp(N2*XTAL_MOV_(o));}
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto &&o)
		const noexcept -> decltype(auto)
		requires in_v<atom::groupoid_q<decltype(o)>>
		{
			return XTAL_ALL_(o)::template zip_from<[]
				XTAL_1FN_(call) (subtype{}.template method<Ns...>)>(XTAL_REF_(o));
		}


	protected:
		template <int N_lim=0>
		XTAL_DEF_(return,inline,let)
		method_approximant(real_variable_q auto o)
		const noexcept -> XTAL_ALL_(o)
		{
			using U_fit = bond::fit<decltype(o)>;
			using U_alpha = typename U_fit::alpha_type;
			using U_sigma = typename U_fit::sigma_type;
			using U_delta = typename U_fit::delta_type;
			using V_delta = int;

			XTAL_IF1_(consteval) {
				auto const d = U_fit::dnsilon_f(1)*half*part_f<signed>(o);
				auto const N = static_cast<U_delta>(o + d);
				auto const n = static_cast<U_alpha>(N);
				method_approximate<N_lim>(o, n);
				auto m = _xtd::bit_cast<U_sigma>(o);
				m += N << U_fit::exponent.shift;
				return _xtd::bit_cast<U_alpha>(m);
			}
			XTAL_0IF_(else) {
				auto const n = round(o);
				auto const N = static_cast<V_delta>(n);
				method_approximate<N_lim>(o, n);
				return ldexp(XTAL_MOV_(o), XTAL_MOV_(N));
			}
		}
		template <int N_lim=0>
		XTAL_DEF_(inline,let)
		method_approximate(real_variable_q auto &o, real_variable_q auto n)
		const noexcept -> void
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
XTAL_TYP_(let) octarithm_t = process::confined_t<octarithm<M_ism, M_div>>;

template <int M_ism=2, int M_div=1, int N_lim=2>
XTAL_DEF_(let) octarithm_f = [] XTAL_1FN_(call) (octarithm_t<M_ism, M_div>{}.template method<N_lim>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
