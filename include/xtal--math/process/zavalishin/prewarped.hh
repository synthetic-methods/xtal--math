#pragma once
#include "./any.hh"

#include "../pade/tangy.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
struct prewarped;


////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
struct prewarped
{
	using superkind = bond::compose<As..., bond::tag<prewarped>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(short)
		XTAL_LET method(auto &&u, real_variable_q auto &&f, auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto const g = pade::tangy_t<1>::template static_method<2>(S_::sample().period()*XTAL_REF_(f));
			return S_::template method<Ns...>(XTAL_REF_(u), g, XTAL_REF_(oo)...);
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
