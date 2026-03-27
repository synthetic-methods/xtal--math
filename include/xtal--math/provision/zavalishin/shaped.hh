#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::provision::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

namespace _detail
{///////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct shaper
{
	template <class S>
	class subtype : public S
	{
	public:
		using S::S;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&...oo)
		noexcept -> decltype(auto)
		{
			return S::template method_f<Ns...>(XTAL_REF_(oo)...);
		};

	};
};


}///////////////////////////////////////////////////////////////////////////////

template <template <auto ...> class A_=_detail::shaper>	struct  shaped;
template <template <auto ...> class A_=_detail::shaper>	using   shaped_t = confined_t<shaped<A_>>;
template <                                 class ..._s>	concept shaped_q = bond::tab_inner_p<shaped<>, _s...>;


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Provides the member-`shaper_t`.
*/
template <template <auto ...> class A_>
struct shaped
{
	using superkind = bond::tab<shaped<>>;
	
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		
	public:
		using S_::S_;
		
		template <auto ...Ms>
		using shaper_t = process::confined_t<A_<Ms...>>;

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
