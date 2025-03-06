#pragma once
#include "./any.hh"

#include "./filter.hh"
#include "../pade/tangy.hh"
#include "../taylor/logarithm.hh"


XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Scales the frequency using the input/state difference and the supplied `reshape`.
\note    Input is restricted to `U_pole` because the filter-state is managed out-of-band.
*/
template <auto ...As>	struct  vactrol;
template <auto ...As>	using   vactrol_t = process::confined_t<vactrol<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <auto ...As>
struct any<vactrol<As...>>
{
	template <class S>
	class subtype : public S
	{
	public:// CONSTRUCT
		using S::S;

		template <extent_type N_mask=-1>
		struct attach
		{
			template <class R>
			using subtype = bond::compose_s<R,
				typename R::reshape_type::template attach<N_mask>>;

		};

	};
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...As>
struct vactrol
{
	using superkind = typename any_t<vactrol>::template attach<>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::reshape_type;
		using typename S_::  state_type;

	public:// OPERATE

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto x, auto s_gain, auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto constexpr abs = [] XTAL_1FN_(call) (taylor::logarithm_t<-1>::template method_f<0>);

			auto const &[d0, d1] = S_::template head<reshape_type>().head();
			auto const  [s_]     = S_::template memory<state_type>();

			auto const v = pade::tangy_t<1>::template method_f< 1>(half*d1);
			auto const w = abs(x - s_.sum());
			s_gain /= d0*term_f(v, one - v, w);

			return S_::template method<Ns...>(x, s_gain, XTAL_REF_(oo)...);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
