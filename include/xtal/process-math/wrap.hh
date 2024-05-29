#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ...As> XTAL_TYP wrap;
template <typename ...As> XTAL_USE wrap_t = process::confined_t<wrap<As...>>;
template <typename ...As> XTAL_FN2 wrap_f(auto &&o)
	XTAL_0EX {return wrap_t<As...>::function(XTAL_REF_(o));};


////////////////////////////////////////////////////////////////////////////////

template <typename ...As> requires (1 <= sizeof...(As))
struct wrap<As...>
:	process::chain<wrap<>, As...>
{};
template <typename ...As> requires (0 == sizeof...(As))
struct wrap<As...>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_FN2 function(auto const &o)
		XTAL_0EX
		{
			using _std::round;

			return o - round(o);
		}
		template <auto ...>
		XTAL_FN2 function(algebra::differential::circular_q auto &&o)
		XTAL_0EX
		{
			return XTAL_REF_(o);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
