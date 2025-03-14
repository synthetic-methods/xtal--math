#pragma once
#include "./any.hh"

#include "./aspect.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Provides the opposing value w.r.t. unit-addition when `N_sel&1`.
\todo    Define via `crossfade` or something?
*/
template <auto ...Ms>
struct oppose;


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
struct oppose
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_sel=1>
		XTAL_DEF_(return,inline,set)
		method_f(auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			auto constexpr I_sel = N_sel&1;
			XTAL_IF0
			XTAL_0IF (I_sel == 0) {return       aspect_t<unsigned>::method_f(XTAL_REF_(o));}
			XTAL_0IF (I_sel == 1) {return one - aspect_t<unsigned>::method_f(XTAL_REF_(o));}
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <auto ...Ms>
using oppose_t = process::confined_t<oppose<Ms...>>;

template <int N_sel=1>
XTAL_DEF_(let)
oppose_f = [] XTAL_1FN_(call) (oppose_t<>::template method_f<N_sel>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
