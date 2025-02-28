#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::provision
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

template <template <auto ...> class A_=_detail::saturator>	struct  saturated;
template <template <auto ...> class A_=_detail::saturator>	using   saturated_t = confined_t<saturated<A_>>;
template <                                   class ..._s>	concept saturated_q = bond::tabbed_with_p<saturated<>, _s...>;


////////////////////////////////////////////////////////////////////////////////
///\
Provides a specialization of `atom::store`. \

template <template <auto ...> class A_>
struct saturated
{
	using superkind = bond::tab<saturated<>>;
	
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
