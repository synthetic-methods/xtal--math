#pragma once
#include "./any.hh"

#include "./meta.hh"
#include "../pade/tangy.hh"
#include "../taylor/tangy.hh"


XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Differentiates the incoming signal.

If the modulator and/or phasor are supplied,
the result is renormalized w.r.t. the chain-rule of differentiation.

\note    The `method` uses local arena-like allocation to manage state.
\note    he maximum `sizeof(input) + sizeof(modulator) <= 64`.
\todo    Enable logarithmic differentiation by including the original signal with the normalization factor?
*/
template <typename ...As>	struct  differ;
template <typename ...As>	using   differ_t = process::confined_t<differ<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <class ..._s>
struct differ
{
	using cotype = occur::meta_t<differ>;

	using data_type = typename cotype::data_type;

	using superkind = bond::compose<bond::tag<differ_t>
	,	provision::memorized<data_type[2]>
	,	meta<_s...>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;
		
	public:// FUNC*

		template <int N_ord=1, int N_sub=0> requires un_v<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &x, auto &&...oo)
		const noexcept -> auto
		{
			return x;
		}
		template <int N_ord=1, int N_sub=0> requires in_v<N_ord, 1> and (1 == data_type::size())
		XTAL_DEF_(return,inline,let)
		method(auto x)
		const noexcept -> auto
		{
			using X = XTAL_ALL_(x);
			using U_fit = bond::fit<X>;
			XTAL_IF0
			XTAL_0IF (N_sub == 0) {
				auto  [x_] = S_::memory(x);
				return x - x_;
			}
			XTAL_0IF (N_sub != 0) {
				auto constexpr m = warped_f<U_fit, N_ord, N_sub>();
				auto  [x_] = S_::template memory<X>();
				x  -=  x_;
				x  *=  m ;
				x_ +=  x ;
				return x ;
			}
		}
		template <int N_ord=1, int N_sub=0> requires in_v<N_ord, 1> and (2 == data_type::size())
		XTAL_DEF_(return,inline,let)
		method(auto x)
		const noexcept -> XTAL_ALL_(x)
		{
			using X = XTAL_ALL_(x);
			using U_fit = bond::fit<X>;
			auto constexpr u = warp_f<U_fit, N_ord, N_sub>();
			auto constexpr m = two/term_f(one, u, u + root_f<2>(U_fit::alpha_f(2)));

			auto [s0, s1] = S_::template memory<X, X>();
			auto const y1 = m*term_f(s1, x - s0, u); s1 = y1 - XTAL_MOV_(s1);
			auto const y0 =   term_f(s0*two, y1, u); s0 = y0 - XTAL_MOV_(s0);
			return u*y1;
		}
		template <int N_ord=1, int N_sub=0> requires in_v<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto x, atom::math::phason_q auto const &t_)
		const noexcept -> auto
		{
			return method<N_ord, N_sub>(XTAL_MOV_(x))*root_f<-1>(t_(1));
		}
		template <int N_ord=1, int N_sub=0> requires in_v<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto x, auto y, atom::math::phason_q auto const &t_)
		const noexcept -> auto
		{
			using X = XTAL_ALL_(x);
			using Y = XTAL_ALL_(y);
			using U_fit = bond::fit<X, Y>;
			XTAL_IF0
			XTAL_0IF (N_sub == 0) {
				auto   [x_, y_] = S_::memory(x, y);
				return (x - x_)*root_f<-1>(t_(1) + (y - y_));
			}
			XTAL_0IF (N_sub != 0) {
				auto constexpr m = warped_f<U_fit, N_ord, N_sub>();
				auto [x_, y_] =  S_::template memory<X, Y>();
				x  -= x_; y  -=  y_;
				x  *= m ; y  *=  m ;
				x_ += x ; y_ +=  y ;
				return x*root_f<-1>(y + t_(1));
			}
		}

	protected:
		template <class U_fit, int N_ord=1, int N_sub=0> requires in_v<N_ord, 1>
		XTAL_DEF_(return,inline,set)
		warped_f()
		noexcept -> auto
		{
			auto constexpr u = warp_f<U_fit, N_ord, N_sub>();
			return (two*u)/(one + u);
		}
		template <class U_fit, int N_ord=1, int N_sub=0>
		XTAL_DEF_(return,inline,set)
		warp_f()
		noexcept -> auto
		{
			static_assert(0 <= N_sub);
			return pade::tangy_f<1>(U_fit::haplo_f(2 + N_sub));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////


namespace xtal::occur
{////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <class ..._s>
struct meta<process::math::zavalishin::differ<_s...>>
{
	using superkind = meta<process::math::zavalishin::meta<_s...>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;
	
	public:
		using S_::S_;

		template <extent_type N_mask=1>
		struct dispatch : bond::compose<void
		,	provision::voiced<void
			,	typename T_::   order_attribute::template dispatch<N_mask>
			>
		,	typename S_::template dispatch<N_mask>
		>
		{};

	};
};


///////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
