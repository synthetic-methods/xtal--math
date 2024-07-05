#pragma once
#include "./any.hh"

#include "./unity.hh"
#include "../gudermannian/sigmoid.hh"



XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Defines `(-1)^(2 #) &`; spiritually equivalent to `1^# &`. \

///\param M_ism \f$\in {1, 2}\f$ specifies the underlying morphism, \
generating either circular or hyperbolic `{cosine, sine}` pairs. \

template <int M_ism=1, typename ...As> XTAL_TYP tangy {static_assert(M_ism);};
template <int M_ism=1, typename ...As> XTAL_USE tangy_t = process::confined_t<tangy<M_ism, As...>>;


////////////////////////////////////////////////////////////////////////////////

template <int M_ism, bond::compose_q ...As> requires some_q<As...>
struct tangy<M_ism, As...>
:	process::lift<tangy<M_ism>, bond::compose<As...>>
{
};
template <int M_ism> requires (0 < M_ism)
struct tangy<M_ism>
{
	static constexpr int I_sgn = sign_n<(M_ism&1)^1, -1>;

	using limit_type = occur::math::limit_t<2>;

	using subkind = bond::compose<void
	,	limit_type::dispatch<>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return)
		XTAL_SET function(simplex_field_q auto &&o)
		XTAL_0EX
		{
			using _op = bond::operate<decltype(o)>;
			using namespace _std;
			
			XTAL_IF0
			XTAL_0IF (N_lim <  0) {
				XTAL_IF0
				XTAL_0IF (1 == M_ism) {return tan (XTAL_REF_(o)*_op::patio_1);}
				XTAL_0IF (2 == M_ism) {return tanh(XTAL_REF_(o)*_op::patio_1);}
			}
			XTAL_0IF (N_lim == 0) {
				return _op::patio_1*gudermannian::sigmoid_t<M_ism>::template function<N_lim>(XTAL_REF_(o));
			}
			XTAL_0IF (0 == (N_lim&1)) {
				auto const [x1, y1] = involved_f(_detail::subunity_t<M_ism,-0>::template function<N_lim>(o*_op::haplo_1));
				auto const x2 =  square_f(x1) + I_sgn*square_f(y1);
				auto const y2 = _op::diplo_1*x1*y1;
				return y2*root_f<-1, 1>(x2);
			}
			XTAL_0IF (1 == (N_lim&1)) {
				auto const [x1, y1] = involved_f(_detail::subunity_t<M_ism,-0>::template function<N_lim>(o));
				return y1*root_f<-1, 1>(x1);
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
