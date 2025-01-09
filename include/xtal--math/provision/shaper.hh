#pragma once
#include "./any.hh"

#include "../process/identity.hh"




XTAL_ENV_(push)
namespace xtal::provision
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <template <auto ...> class F=process::math::identity> struct   shaper;
template <template <auto ...> class F=process::math::identity> using    shaper_t = confined_t<shaper<F>>;
template <                                        class ..._s> concept  shaper_q = bond::tab_p<shaper<>, _s...>;


////////////////////////////////////////////////////////////////////////////////
///\
Provides a specialization of `arrange::store`. \

template <template <auto ...> class F>
struct shaper
{
	using superkind = bond::tab<shaper<>>;
	
	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;
		
	public:
		using S_::S_;
		
		template <auto ...Ms>
		using shape_t = process::confined_t<F<Ms...>>;

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
