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
:	bond::compose<dilated<XTAL_VAL_(-bond::operating::patio_2)>, sine<M_ism, -0>>
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
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&w)
		noexcept -> decltype(auto)
		{
			if constexpr (N_lim < 0) {
				using namespace _std;
				XTAL_IF0
				XTAL_0IF (M_ism ==  2) {return  sinh(XTAL_REF_(w));}
				XTAL_0IF (M_ism ==  1) {return  sin (XTAL_REF_(w));}
				XTAL_0IF (M_ism == -1) {return asin (XTAL_REF_(w));}
				XTAL_0IF (M_ism == -2) {return asinh(XTAL_REF_(w));}
			}
			else {
				return S_::template function<N_lim>(XTAL_REF_(w));
			}
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
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
		{
			using _op = bond::operate<decltype(o)>;

			auto const w = XTAL_REF_(o)*_op::haplo_f(N_lim);
			auto const u = w*_op::haplo_1;

			if constexpr (N_lim < 0) {
				using namespace _std;
				XTAL_IF0
				XTAL_0IF (M_ism ==  2) {return  sinh(o)/XTAL_REF_(o);}
				XTAL_0IF (M_ism ==  1) {return  sin (o)/XTAL_REF_(o);}
				XTAL_0IF (M_ism == -1) {return asin (o)/XTAL_REF_(o);}
				XTAL_0IF (M_ism == -2) {return asinh(o)/XTAL_REF_(o);}
			}
			else {
				return S_::template function<N_lim>(XTAL_REF_(o));
			}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (0 < M_ism)
struct sine<M_ism, -2>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(long,static)
		XTAL_LET function(auto &&w)
		noexcept -> auto
		{
			if constexpr (N_lim < 0) {
				using namespace _std;
				XTAL_IF0
				XTAL_0IF (M_ism ==  2) {return sinh(root_f<2>(w))/XTAL_REF_(w);}
				XTAL_0IF (M_ism ==  1) {return sin (root_f<2>(w))/XTAL_REF_(w);}
			}
			else {
				int constexpr N = (N_lim << 1) - (0 < N_lim);

				using W = XTAL_ALL_(w); using _op = bond::operate<W>;
				W x = _op::alpha_1;

				bond::seek_backward_f<N>([&] (auto i)
					XTAL_0FN_(x = term_f(_op::alpha_1
					,	+_op::ratio_f(1, (2 + 2*i)*(3 + 2*i))
					,	w
					,	x
					)
				));
				return x;
			}
		}

	};
};
template <int M_ism> requires (M_ism < 0)
struct sine<M_ism, -2>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(long,static)
		XTAL_LET function(auto &&w)
		noexcept -> auto
		{
			if constexpr (N_lim < 0) {
				using namespace _std;
				XTAL_IF0
				XTAL_0IF (M_ism == 1) {return asin (root_f<2>(XTAL_REF_(w)));}
				XTAL_0IF (M_ism == 2) {return asinh(root_f<2>(XTAL_REF_(w)));}
			}
			else {
				int constexpr N = (N_lim << 1) - (0 < N_lim);

				using W = XTAL_ALL_(w); using _op = bond::operate<W>;
				W x = _op::ratio_f(1, 1 + 2*N);

				bond::seek_backward_f<N>([&] (auto i)
					XTAL_0FN_(x = term_f(_op::ratio_f(1, 1 + 2*i)
					,	-_op::template ratio_f(1 + 2*i, 2 + 2*i)
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
