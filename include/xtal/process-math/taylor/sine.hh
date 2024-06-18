#pragma once
#include "./any.hh"

#include "../root.hh"




XTAL_ENV_(push)
namespace xtal::process::math::taylor
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=1, int N_car=0> XTAL_TYP sine {static_assert(M_ism);};
template <int M_ism=1, int N_car=0> XTAL_USE sine_t = process::confined_t<sine<M_ism, N_car>>;

template <int M_ism>
struct sine<M_ism, +1>
:	bond::compose<dilating<0, -1>, sine<M_ism, -0>>
{
};
template <int M_ism>
struct sine<M_ism, -0>
{
	using subkind = bond::compose<discarding<1, +1>, sine<M_ism, -1>>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,static)
		XTAL_RET function(auto &&u)
		XTAL_0EX
		{
			if constexpr (N_lim < 0) {
				using namespace _std;
				XTAL_IF0
				XTAL_0IF (M_ism ==  2) {return  sinh(XTAL_REF_(u));}
				XTAL_0IF (M_ism ==  1) {return  sin (XTAL_REF_(u));}
				XTAL_0IF (M_ism == -1) {return asin (XTAL_REF_(u));}
				XTAL_0IF (M_ism == -2) {return asinh(XTAL_REF_(u));}
			}
			else {
				return S_::template function<N_lim>(XTAL_REF_(u));
			}
		}

	};
};
template <int M_ism>
struct sine<M_ism, -1>
{
	using subkind = bond::compose<discarding<1, +2>, sine<M_ism, -2>>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,inline,static)
		XTAL_RET function(auto &&u)
		XTAL_0EX
		{
			if constexpr (N_lim < 0) {
				using namespace _std;
				XTAL_IF0
				XTAL_0IF (M_ism ==  2) {return  sinh(u)/XTAL_REF_(u);}
				XTAL_0IF (M_ism ==  1) {return  sin (u)/XTAL_REF_(u);}
				XTAL_0IF (M_ism == -1) {return asin (u)/XTAL_REF_(u);}
				XTAL_0IF (M_ism == -2) {return asinh(u)/XTAL_REF_(u);}
			}
			else {
				return S_::template function<N_lim>(XTAL_REF_(u));
			}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism> requires (0 < M_ism)
struct sine<M_ism, -2>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,static)
		XTAL_RET function(auto &&w)
		XTAL_0EX
		{
			if constexpr (N_lim < 0) {
				using namespace _std;
				XTAL_IF0
				XTAL_0IF (M_ism ==  2) {return sinh(root_f<2>(w))/XTAL_REF_(w);}
				XTAL_0IF (M_ism ==  1) {return sin (root_f<2>(w))/XTAL_REF_(w);}
			}
			else {
				int constexpr I_lim = (N_lim << 1) - (0 < N_lim);
				int constexpr I_sgn = sign_n<(M_ism&1)^1, -1>;

				using W = XTAL_ALL_(w); using _op = bond::operate<W>;
				W x = _op::alpha_1;

				bond::seek_backward_f<I_lim>([&] (auto i)
					XTAL_0FN_(x = horner::term_f<I_sgn>(_op::alpha_1
					,	_op::ratio_f(1, (2 + 2*i)*(3 + 2*i))
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
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>
		XTAL_DEF_(return,static)
		XTAL_RET function(auto &&w)
		XTAL_0EX
		{
			if constexpr (N_lim < 0) {
				using namespace _std;
				XTAL_IF0
				XTAL_0IF (M_ism == 1) {return asin (root_f<2>(XTAL_REF_(w)));}
				XTAL_0IF (M_ism == 2) {return asinh(root_f<2>(XTAL_REF_(w)));}
			}
			else {
				int constexpr I_lim = (N_lim << 1) - (0 < N_lim);
				int constexpr I_sgn = sign_n<(M_ism&1)^0, -1>;

				using W = XTAL_ALL_(w); using _op = bond::operate<W>;
				W x = _op::ratio_f(1, 1 + 2*I_lim);

				bond::seek_backward_f<I_lim>([&] (auto i)
					XTAL_0FN_(x = horner::term_f<I_sgn>(_op::ratio_f(1, 1 + 2*i)
					,	_op::template ratio_f(1 + 2*i, 2 + 2*i)
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
