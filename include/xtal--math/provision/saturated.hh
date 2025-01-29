#pragma once
#include "./any.hh"

#include "../process/identity.hh"




XTAL_ENV_(push)
namespace xtal::provision
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <template <auto ...> class F=process::math::identity> struct   saturated;
template <template <auto ...> class F=process::math::identity> using    saturated_t = confined_t<saturated<F>>;
template <                                        class ..._s> concept  saturated_q = bond::tab_p<saturated<>, _s...>;


////////////////////////////////////////////////////////////////////////////////
///\
Provides a specialization of `atom::store`. \

template <template <auto ...> class F>
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
		using saturate_t = process::confined_t<F<Ms...>>;

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
