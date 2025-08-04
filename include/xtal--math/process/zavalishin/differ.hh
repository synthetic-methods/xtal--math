#pragma once
#include "./any.hh"

#include "./scaffold.hh"




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
	using metatype = traits_t<differ>;

	using state_type = typename metatype::state_type;
//	using slope_type = typename metatype::slope_type;

	using superkind = bond::compose<bond::tag<differ_t>
	,	provision::memorized<state_type>
	,	scaffold<_s...>
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
		method(auto const &u, auto &&...oo) const
		noexcept -> auto
		{
			return u;
		}
		template <int N_ord=1> requires in_n<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &u) const
		noexcept -> auto
		{
			auto [u_] = S_::memory(u);
			return (u - u_);
		}
		template <int N_ord=1> requires in_n<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &u, atom::math::phason_q auto const &t_) const
		noexcept -> auto
		{
			auto [u_] = S_::memory(u);
			return (u - u_)*root_f<-1>(t_(1));
		}
		template <int N_ord=1> requires in_n<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &u, auto const &v, atom::math::phason_q auto const &t_) const
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

#include "./differ.hh_"
