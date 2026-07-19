#pragma once
#include "./any.hh"

#include "../root.hh"




XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Half-period approximation of `{Cos[Pi #], Cosh[Pi #]}`, indexed by `M_ism`.

Evaluated using the adjusted Taylor series:
```m
(1 - Cosh[Sqrt[#1 - 1]*(Pi/2)]/#)/(1 - #) &;
CoefficientList[Normal@Series[%@x, {x, 0, 8}], x]
```
*/
template <int M_ism=1, int N_car=0>	struct  cosy {static_assert(M_ism);};
template <int M_ism=1, int N_car=0>	using   cosy_t = process::confined_t<cosy<M_ism, N_car>>;

template <int M_ism=1, int M_car=0, int ...Ns>
XTAL_DEF_(let)
cosy_f = [] XTAL_1FN_(call) (cosy_t<M_ism, M_car>::template method<Ns...>);


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (0 < M_ism)
struct cosy<M_ism, +2>
{
	using superkind = cosy<M_ism, -0>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,let)
		method(auto &&o
		)	const
		noexcept -> decltype(auto)
		{
			return square_f(S_::template method<N_lim>(XTAL_REF_(o)));
		}

	};
};
template <int M_ism> requires (0 < M_ism)
struct cosy<M_ism, -0>
{
	using superkind = bond::compose<typename discard_t<2>::template infix<>, cosy<M_ism, -2>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,let)
		method(auto &&o
		)	const
		noexcept -> decltype(auto)
		{
			using U_fit = bond::fit<decltype(o)>;
			XTAL_IF0
			XTAL_0IF_(consteval)  {return S_::template method<   ~0>(XTAL_REF_(o));}
			XTAL_0IF (0 <= N_lim) {return S_::template method<N_lim>(XTAL_REF_(o));}
			XTAL_0IF (M_ism == 2) {return cosh (U_fit::patio_1*XTAL_REF_(o));}
			XTAL_0IF (M_ism == 1) {return cos  (U_fit::patio_1*XTAL_REF_(o));}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (0 < M_ism)
struct cosy<M_ism, -2>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	public:
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		method(auto &&w)
		noexcept -> auto
		{
			using U_fit = bond::fit<decltype(w)>;
			XTAL_IF0
			XTAL_0IF (0 <= N_lim)  {return methox<N_lim>(XTAL_REF_(w));}
			XTAL_0IF_(consteval)   {return methox<   ~0>(XTAL_REF_(w));}
			XTAL_0IF (M_ism ==  2) {return cosh (U_fit::patio_1*root_f<2>(XTAL_REF_(w)));}
			XTAL_0IF (M_ism ==  1) {return cos  (U_fit::patio_1*root_f<2>(XTAL_REF_(w)));}
		}

	protected:
		template <int N_lim>
		XTAL_DEF_(return,inline,set)
		methox(auto &&w)
		noexcept -> auto
		{
			static_assert(-1 <= N_lim);
			using W = XTAL_ALL_(w);
			using U = unstruct_t<W>;
			W const w0 = U(4*cosign_v<M_ism>)*XTAL_REF_(w);
			W const w1 = U(1) + w0;
			W m(one);
			switch (N_lim&0b11) {
			case 3:
				m = term_f(one, U{0.00811288571346782540543271186180452L}*w1, XTAL_MOV_(m));
				m = term_f(one, U{0.01035856836301015154103279216061790L}*w1, XTAL_MOV_(m));
			case 2:
				m = term_f(one, U{0.01368644468307770278761249222497005L}*w1, XTAL_MOV_(m));
				m = term_f(one, U{0.01892392299447479546768315655526088L}*w1, XTAL_MOV_(m));
			case 1:
				m = term_f(one, U{0.02787253570624807830822963509172920L}*w1, XTAL_MOV_(m));
				m = term_f(one, U{0.04509227375660560857091132405389036L}*w1, XTAL_MOV_(m));
			case 0:
				m = term_f(one, U{0.08505190841862807638032318268573259L}*w1, XTAL_MOV_(m));
				m = term_f(one, U{0.21460183660255169038433915418012428L}*w0, XTAL_MOV_(m));
			}
			return m*w1;
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
