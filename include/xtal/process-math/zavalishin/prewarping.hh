#pragma once
#include "./any.hh"

#include "../gudermannian/all.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
struct prewarping;


////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
struct prewarping
{
	using subkind = bond::compose<As..., resource::example<>, bond::tag<prewarping>>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_RET functor(auto &&u, auto &&f, auto &&...oo)
		XTAL_0EX
		{
			auto const t = gudermannian::tangent_f<1>(S_::sample().period()*XTAL_REF_(f));
			return S_::template functor<Is...>(XTAL_REF_(u), t, XTAL_REF_(oo)...);
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
