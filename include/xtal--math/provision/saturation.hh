#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::provision::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

namespace _detail
{///////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct saturator
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

template <template <auto ...> class A_=_detail::saturator>	struct  saturation;
template <template <auto ...> class A_=_detail::saturator>	using   saturation_t = confined_t<saturation<A_>>;
template <                                    class ..._s>	concept saturation_q = bond::tab_in_p<saturation<>, _s...>;


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Provides the member-`saturate_t`.
*/
template <template <auto ...> class A_>
struct saturation
{
	using superkind = bond::tab<saturation<>>;
	
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		
	public:
		using S_::S_;
		
		template <auto ...Ms>
		using saturate_t = process::confined_t<A_<Ms...>>;

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
