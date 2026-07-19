#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::scheme::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

namespace _detail
{///////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct distortion
{
	template <class S>
	class subtype : public S
	{
	public:
		using S::S;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto &&...oo)
		const noexcept -> decltype(auto)
		{
			return S::template method<Ns...>(XTAL_REF_(oo)...);
		};

	};
};


}///////////////////////////////////////////////////////////////////////////////

template <template <auto ...> class A_=_detail::distortion>	XTAL_TYP_(new) distorted;
template <template <auto ...> class A_=_detail::distortion>	XTAL_TYP_(let) distorted_t = confined_t<distorted<A_>>;
template <                                     class ..._s>	XTAL_TYP_(ask) distorted_q = bond::tab_inner_p<distorted<>, _s...>;


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Provides the member-`distortion_t`.
*/
template <template <auto ...> class A_>
struct distorted
{
	using superkind = bond::tab<distorted<>>;
	
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		
	public:
		using S_::S_;
		
		template <auto ...Ms>
		using distortion_t = process::confined_t<A_<Ms...>>;

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
