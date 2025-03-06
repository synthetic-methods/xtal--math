#pragma once
#include "./any.hh"






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
struct any<differ<_s...>>
{
	using superkind = any<class_template<_s...>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;
	
	public:
		using S_::S_;
		//\

		template <extent_type N_mask=-1>
		struct dispatch
		{
			template <class R>
			using subtype = bond::compose_s<R, provision::voiced<void
			,	typename T_::order_type::template dispatch<N_mask>
			>>;

		};

	};
};
template <scalar_q A>
struct any<differ<A>> : any<differ<A[1]>>
{
};
template <>
struct any<differ< >> : any<differ<typename bond::fit<>::alpha_type>>
{
};


////////////////////////////////////////////////////////////////////////////////

template <class ..._s>
struct differ
{
	using metakind = any  <differ>;
	using metatype = any_t<differ>;

	using state_type = typename metatype::state_type;
//	using slope_type = typename metatype::slope_type;

	using superkind = bond::compose<bond::tag<differ_t>
	,	metakind
	,	provision::memorized<state_type>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;
		
	public:// FUNC*

		template <int N_ord=1> requires un_n<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &u, auto &&...oo)
		noexcept -> auto
		{
			return u;
		}
		template <int N_ord=1> requires in_n<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &u)
		noexcept -> auto
		{
			auto [u_] = S_::memory(u);
			return (u - u_);
		}
		template <int N_ord=1> requires in_n<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &u, atom::math::phason_q auto const &t_)
		noexcept -> auto
		{
			auto [u_] = S_::memory(u);
			return (u - u_)*root_f<-1>(t_(1));
		}
		template <int N_ord=1> requires in_n<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &u, auto const &v, atom::math::phason_q auto const &t_)
		noexcept -> auto
		{
			auto [u_, v_] = S_::memory(u, v);
			return (u - u_)*root_f<-1>(t_(1) + v - v_);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
