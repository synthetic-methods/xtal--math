#pragma once
#include "./any.hh"

#include "./filter.hh"
#include "../taylor/logarithm.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Scales `damp` by the dot-product of the internal state and supplied `reshape`.
\note    Input is restricted to `U_pole` because the filter-state is managed out-of-band.
*/
template <auto ...As>	struct  vectrol;
template <auto ...As>	using   vectrol_t = process::confined_t<vectrol<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <auto ...As>
struct any<vectrol<As...>>
{
	template <class S>
	class subtype : public S
	{
	public:// CONSTRUCT
		using S::S;

		template <size_type N_mask=1>
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
struct vectrol
{
	using superkind = typename any_t<vectrol>::template attach<>;

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
		method(auto x, auto s_gain, auto s_damp, auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto constexpr abs = [] XTAL_1FN_(call) (taylor::logarithm_t<-1>::template method_f<0>);

			auto const [s_]       = S_::template memory<  state_type>();
			auto const &s_reshape = S_::template   head<reshape_type>().head();
			auto const  s_product = dot_f(s_, s_reshape);
			
			s_damp *= abs(s_product*half);

			return S_::template method<Ns...>(x, s_gain, s_damp, XTAL_REF_(oo)...);
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
