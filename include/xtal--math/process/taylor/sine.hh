#pragma once
#include "./any.hh"

#include "../root.hh"




XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int N_car=0> struct   sine {static_assert(M_ism);};
template <int M_ism=1, int N_car=0> using    sine_t = process::confined_t<sine<M_ism, N_car>>;

template <int M_ism>
struct sine<M_ism, +1>
:	bond::compose<dilated<[] XTAL_0FN_(value) (-bond::fixture<>::patio_2)>, sine<M_ism, -0>>
{
};
template <int M_ism>
struct sine<M_ism, -0>
{
	using superkind = bond::compose<discarded<1>, sine<M_ism, -1>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		static_method(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF_(consteval)   {return S_::template static_method<   ~0>(XTAL_REF_(o));}
			XTAL_0IF (0 <= N_lim)  {return S_::template static_method<N_lim>(XTAL_REF_(o));}
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
	using superkind = bond::compose<discarded<2, M_ism>, sine<M_ism, -2>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,set)
		static_method(auto &&o)
		noexcept -> decltype(auto)
		{
			XTAL_IF0
			XTAL_0IF_(consteval)   {return S_::template static_method<   ~0>(XTAL_REF_(o));}
			XTAL_0IF (0 <= N_lim)  {return S_::template static_method<N_lim>(XTAL_REF_(o));}
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
		XTAL_DEF_(return,inline,set)
		static_method(auto &&w)
		noexcept -> auto
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
		XTAL_DEF_(return,set)
		approximate(auto &&w)
		noexcept -> auto
		{
			auto constexpr N = 1 + 2*below_m<0x10, (unsigned) N_lim>;

			using W = XTAL_ALL_(w); using _fix = bond::fixture<W>;
			XTAL_IF0
			XTAL_0IF (0 < M_ism) {
				W x{one};

				bond::seek_backward_f<N>([&] (auto i) XTAL_0FN_(return) (
					x = term_f(one
					,	+_fix::ratio_f(1, (2 + 2*i)*(3 + 2*i))
					,	w
					,	x
					)
				));
				return x;
			}
			XTAL_0IF (M_ism < 0) {
				W x = _fix::ratio_f(1, 1 + 2*N);

				bond::seek_backward_f<N>([&] (auto i) XTAL_0FN_(return) (
					x = term_f(_fix::ratio_f(1, 1 + 2*i)
					,	-_fix::ratio_f(1 + 2*i, 2 + 2*i)
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
