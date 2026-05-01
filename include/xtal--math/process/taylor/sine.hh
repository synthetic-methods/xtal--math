#pragma once
#include "./any.hh"

#include "../root.hh"




XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int N_car=0>	struct  sine {static_assert(M_ism);};
template <int M_ism=1, int N_car=0>	using   sine_t = process::confined_t<sine<M_ism, N_car>>;

template <int M_ism>
struct sine<M_ism, +1>
:	bond::compose<typename dilate_t<[] XTAL_1FN_(to) (-bond::fit<>::patio_2)>::template infix<>, sine<M_ism, -0>>
{
};
template <int M_ism>
struct sine<M_ism, -0>
{
	using superkind = bond::compose<typename discard_t<1>::template infix<>, sine<M_ism, -1>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,let)
		method(auto &&o)
		const noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF_(consteval)   {return S_::template method<   ~0>(XTAL_REF_(o));}
			XTAL_0IF (0 <= N_lim)  {return S_::template method<N_lim>(XTAL_REF_(o));}
			XTAL_0IF (M_ism ==  2) {return  sinh(XTAL_REF_(o));}
			XTAL_0IF (M_ism ==  1) {return  sin (XTAL_REF_(o));}
			XTAL_0IF (M_ism == -1) {return asin (XTAL_REF_(o));}
			XTAL_0IF (M_ism == -2) {return asinh(XTAL_REF_(o));}
		}

	};
};
template <int M_ism>
struct sine<M_ism, -1>
{
	using superkind = bond::compose<typename discard_t<2, M_ism>::template infix<>, sine<M_ism, -2>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,let)
		method(auto &&o)
		const noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF_(consteval)   {return S_::template method<   ~0>(XTAL_REF_(o));}
			XTAL_0IF (0 <= N_lim)  {return S_::template method<N_lim>(XTAL_REF_(o));}
			XTAL_0IF (M_ism ==  2) {return  sinh(o)/o;}
			XTAL_0IF (M_ism ==  1) {return  sin (o)/o;}
			XTAL_0IF (M_ism == -1) {return asin (o)/o;}
			XTAL_0IF (M_ism == -2) {return asinh(o)/o;}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism>
struct sine<M_ism, -2>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	public:
		template <int N_lim=-1>
		XTAL_DEF_(return,inline,let)
		method(auto &&w)
		const noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (0 == N_lim)  {return XTAL_ALL_(w) {one};}
			XTAL_0IF (1 <= N_lim)  {return approximate<N_lim>(XTAL_REF_(w));}
			XTAL_0IF_(consteval)   {return approximate<   ~0>(XTAL_REF_(w));}
			XTAL_0IF_(else) {
				auto const u = root_f<2>(w);
				XTAL_IF0
				XTAL_0IF (M_ism ==  2) {return sinh (u)/u;}
				XTAL_0IF (M_ism ==  1) {return sin  (u)/u;}
				XTAL_0IF (M_ism == -1) {return asin (u)/u;}
				XTAL_0IF (M_ism == -2) {return asinh(u)/u;}
			}
		}

	protected:
		template <int N_lim>
		XTAL_DEF_(return,inline,let)
		approximate(auto &&w)
		const noexcept -> auto
		{
			static_assert(-1 <= N_lim);
			int constexpr N = term_f(1, 2, N_lim&0b1111);

			using W = XTAL_ALL_(w); using _fit = bond::fit<W>;
			XTAL_IF0
			XTAL_0IF (0 < M_ism) {
				W x{one};

				bond::seek_to_e<-N>([&] (auto I) XTAL_0FN_(to) (
					x = term_f(one
					,	+_fit::ratio_f(1, (2 + 2*I)*(3 + 2*I))
					,	w
					,	x
					)
				));
				return x;
			}
			XTAL_0IF (M_ism < 0) {
				W x = _fit::ratio_f(1, 1 + 2*N);

				bond::seek_to_e<-N>([&] (auto I) XTAL_0FN_(to) (
					x = term_f(_fit::ratio_f(1, 1 + 2*I)
					,	-_fit::ratio_f(1 + 2*I, 2 + 2*I)
					,	w
					,	x
					)
				));
				return x;
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
