#pragma once
#include "./any.hh"

#include "./meta.hh"




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
	using cotype = occur::codex_t<differ>;

	using data_type = typename cotype::data_type;

	using superkind = bond::compose<bond::tag<differ_t>
	,	provision::memorized<data_type>
	,	meta<_s...>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;
		
	public:// FUNC*

		template <int N_ord=1, int N_nyq=0> requires un_v<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto const &u, auto &&...oo)
		const noexcept -> auto
		{
			return u;
		}
		template <int N_ord=1, int N_nyq=0> requires in_v<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto u)
		const noexcept -> auto
		{
			if constexpr (N_nyq == 0) {
				auto  [u_] = S_::memory(u);
				return (u - u_);
			}
			else {
				using U = XTAL_ALL_(u);
				auto constexpr f = fact_f<N_ord, N_nyq>(u);

				auto [_u] = S_::template memory<U>();
				u  -= _u;
				u  *=  f;
				_u +=  u;
				return u;
			}
		}
		template <int N_ord=1, int N_nyq=0> requires in_v<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto u, atom::math::phason_q auto const &t_)
		const noexcept -> auto
		{
			return method<N_ord, N_nyq>(XTAL_MOV_(u))*root_f<-1>(t_(1));
		}
		template <int N_ord=1, int N_nyq=0> requires in_v<N_ord, 1>
		XTAL_DEF_(return,inline,let)
		method(auto u, auto v, atom::math::phason_q auto const &t_)
		const noexcept -> auto
		{
			if constexpr (N_nyq == 0) {
				auto   [u_, v_] = S_::memory(u, v);
				return (u - u_)*root_f<-1>(t_(1) + v - v_);
			}
			else {
				using U = XTAL_ALL_(u);
				using V = XTAL_ALL_(v);
				auto constexpr f = fact_f<N_ord, N_nyq>(u, v);

				auto [_u, _v] = S_::template memory<U, V>();
				u  -= _u; v  -= _v;
				u  *=  f; v  *=  f;
				_u +=  u; _v +=  v;
				return u*root_f<-1>(v + t_(1));
			}
		}

	protected:
		template <int N_ord=1, int N_nyq=0> requires in_v<N_ord, 1>
		XTAL_DEF_(return,inline,set)
		fact_f(auto const &...oo)
		noexcept -> auto
		{
			using U_fit = bond::fit<decltype(oo)...>;
			auto constexpr w = pade::tangy_f<1>(U_fit::haplo_f(1 + N_nyq));
			return (two*w)/(one + w);
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
struct codex<process::math::zavalishin::differ<_s...>>
{
	using superkind = codex<process::math::zavalishin::meta<_s...>>;

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
