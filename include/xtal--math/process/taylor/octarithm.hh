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

template <int M_ism, int M_div> requires in_n<M_ism,  1,  2>
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

		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		method_f(auto o)
		noexcept -> XTAL_ALL_(o)
		{
			using U_alpha = typename bond::fit<decltype(o)>::alpha_type;
			U_alpha constexpr N_ =  one*U_alpha{M_div};
			U_alpha constexpr N2 =  one*_std::numbers::ln2_v<U_alpha>;
			U_alpha constexpr N1 = -two*_std::numbers:: pi_v<U_alpha>;

			o = S_::template method_f<N_lim>(XTAL_MOV_(o));
			if constexpr (M_div  >  1) {o *=     N_;}
			if constexpr (M_ism ==  1) {o *= one/N1;}
			if constexpr (M_ism ==  2) {o *= one/N2;}
			return o;
		}

	};
};
template <int M_ism, int M_div> requires in_n<M_ism, -1, -2>
struct octarithm<M_ism, M_div>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	//	TODO: Define `complex` variant!

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(atom::groupoid_q auto &&o)
		noexcept -> XTAL_ALL_(o)
		{
			using  U = XTAL_ALL_(o);
			return U::template zip_from<[] XTAL_1FN_(call) (method_f<Ns...>)>(XTAL_REF_(o));
		}
		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		method_f(integral_q auto &&x)
		noexcept -> auto
		{
			using X = _xtd::make_signed_t<unstruct_u<decltype(x)>>;
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
		XTAL_DEF_(return,inline,set)
		method_f(real_q auto o)
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
			XTAL_0IF (0 <= N_lim)          {return approximant_f<N_lim>(XTAL_MOV_(o));}
			XTAL_0IF_(consteval)           {return approximant_f<   ~0>(XTAL_MOV_(o));}
#if XTAL_SYS_(builtin)
			XTAL_0IF (real_q<decltype(o)>) {return       __builtin_exp2(XTAL_MOV_(o));}
#endif
			XTAL_0IF_(else)                {return               exp(N2*XTAL_MOV_(o));}
		}

	protected:
		template <int N_lim=0>
		XTAL_DEF_(return,inline,set)
		approximant_f(real_variable_q auto o)
		noexcept -> XTAL_ALL_(o)
		{
			using U_fit = bond::fit<decltype(o)>;
			using U_alpha = typename U_fit::alpha_type;
			using U_sigma = typename U_fit::sigma_type;
			using U_delta = typename U_fit::delta_type;
			using V_delta = int;

			XTAL_IF1_(consteval) {
				auto const d = U_fit::dnsilon_f(1)*half*decompose_f<signed>(o);
				auto const N = static_cast<U_delta>(o + d);
				auto const n = static_cast<U_alpha>(N);
				approximate_f<N_lim>(o, n);
				auto m = _xtd::bit_cast<U_sigma>(o);
				m += N << U_fit::exponent.shift;
				return _xtd::bit_cast<U_alpha>(m);
			}
			XTAL_0IF_(else) {
				auto const n = round(o);
				auto const N = static_cast<V_delta>(n);
				approximate_f<N_lim>(o, n);
				return ldexp(XTAL_MOV_(o), XTAL_MOV_(N));
			}
		}
		template <int N_lim=0>
		XTAL_DEF_(inline,set)
		approximate_f(real_variable_q auto &o, real_variable_q auto n)
		noexcept -> void
		{
			using U_fit = bond::fit<decltype(o)>;
			using U_alpha = typename U_fit::alpha_type;

			auto constexpr I_lim = 4 + 4*below_v<4, (unsigned) N_lim>;
			auto constexpr U1 = U_fit::haplo_f(I_lim)*_std::numbers::ln2_v<U_alpha>;
			auto constexpr W2 = U_fit::haplo_f(1)*U1*U1;

			o -= n;
			o  = term_f(one, o, term_f(U1, W2, o));
			#pragma unroll
			for (int i{}; i < I_lim; ++i) {
				o *= o;
			}
		}

	};
};

////////////////////////////////////////////////////////////////////////////////

template <int M_ism=2, int M_div=1>
XTAL_TYP_(let) octarithm_t = process::confined_t<octarithm<M_ism, M_div>>;

template <int M_ism=2, int M_div=1, int N_lim=2>
XTAL_DEF_(let) octarithm_f = [] XTAL_1FN_(call) (octarithm_t<M_ism, M_div>::template method_f<N_lim>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
