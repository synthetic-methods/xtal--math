#pragma once
#include "./any.hh"

#include "./dot.hh"
#include "./root.hh"
#include "./roots.hh"


XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Resolves the element indicated by the supplied value.
*/
template <auto ...>
struct component;


////////////////////////////////////////////////////////////////////////////////

template <int M_sel>
struct component<M_sel>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(complex_field_q auto &&o)
		noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (M_sel == -1) {return XTAL_REF_(o)       ;}
			XTAL_0IF (M_sel ==  0) {return XTAL_REF_(o).real();}
			XTAL_0IF (M_sel ==  1) {return XTAL_REF_(o).imag();}
		}

	};
};
template <>
struct component<> : component<-1>
{
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
using component_t = process::confined_t<component<Ms...>>;

template <auto ...Ms>
XTAL_DEF_(let)
component_f = [] XTAL_1FN_(call) (component_t<Ms...>::method_f);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
