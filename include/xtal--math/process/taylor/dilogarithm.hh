#pragma once
#include "./any.hh"

#include "../pade/unity.hh"




XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int M_car=0>
XTAL_TYP_(new) dilogarithm;

template <int M_ism=1, int M_car=0>
XTAL_TYP_(let) dilogarithm_t = process::confined_t<dilogarithm<M_ism, M_car>>;

template <int M_ism=1, int M_car=0, int ...Ns>
XTAL_DEF_(let) dilogarithm_f = [] XTAL_1FN_(call)
	(dilogarithm_t<M_ism, M_car>::template method<Ns...>);


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Approximates the dilogarithm `PolyLog[2, #]&`,
*/
template <int M_car>
struct dilogarithm<1, M_car>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=3>
		XTAL_DEF_(return,inline,set)
		method(complex_field_q auto const &q, complex_field_q auto const &o)
		noexcept -> decltype(auto)
		{
			using U_fit = bond::fit<decltype(o), decltype(q)>;
			if constexpr (N_lim == 0) {
				int constexpr K_lim = 8;
				auto o_i = o;
				auto q_i = q;
				auto z   = q*o.real();

				bond::seek_to_e<-K_lim>([&] (auto const I) XTAL_0FN {
					auto constexpr _K = U_fit::ratio_f(cosign_v<1 + I>, square_f(2 + I));
					z = term_f(XTAL_MOV_(z), (q_i *= q), (o_i *= o).real(), _K);
				});
				return -two*z;
			}
			else {
				return method<N_lim>(imagine_f<2, 0>(o)*(q)) + method<N_lim>(imagine_f<2, 1>(o)*(q));
			}
		}
		template <int N_lim=3>
		XTAL_DEF_(return,inline,set)
		method(complex_field_q auto const &q)
		noexcept -> decltype(auto)
		{
		//	Could file under `spence`, which is `-PolyLog[2,-#]&`
			using U = XTAL_ALL_(q);
			using V = unstruct_t<U>;
			using U_fit = bond::fit<decltype(q)>;
			
			V constexpr _03 = U_fit::ratio_f(1,  3), _Q3 = root_f<2>(_03);
			V constexpr _12 = U_fit::ratio_f(1, 12);
			V constexpr _27 = U_fit::ratio_f(1, 27);

			int constexpr K_lim = 1 << (0b11&N_lim);
		// Precision scales with the square, but...
		//	{0,  1,  2,  3,  4}
		//	{0,  1,  4,  9, 16} (* #^2& *)
		//	{1,  2,  4,  8, 16} (* 2^#& *)

			objective_t<U> dn{one}, up{one};
			dn *= term_f(two*_Q3 - one, _Q3, q);
			dn *= term_f(two*_Q3 + one, _Q3, q);
			bond::seek_to_e<-(K_lim)>([&] (auto const I) XTAL_0FN {
				up = term_f(one, q, XTAL_MOV_(up)*U_fit::template ratio_f<2>(I + 1, I + 4));
			});
			up = term_f(two - _12, _27*q, XTAL_MOV_(up));
			up = term_f(one + _03,     q, XTAL_MOV_(up));
			return root_f<-1>(dn)*term_f(q*up,
				one - square_f(q), imagine_f<1>(pade::unity_f<-1>(one - q))*U_fit::patio_2);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
