#pragma once
#include "./any.hh"

#include "../pade/tangy.hh"




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
		XTAL_RET method(auto &&u, auto &&f, auto &&...oo)
		XTAL_0EX
		{
			using _op = bond::operate<decltype(u)>;
			//\
			auto const t = _std::tan(_op::patio_1*S_::sample().period()*XTAL_REF_(f));
			auto const t = pade::tangy_t<1>::template function<2>(S_::sample().period()*XTAL_REF_(f));
			return S_::template method<Is...>(XTAL_REF_(u), t, XTAL_REF_(oo)...);
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
