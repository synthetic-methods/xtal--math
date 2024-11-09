#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ...As> XTAL_TYP wrap;
template <typename ...As> XTAL_USE wrap_t = process::confined_t<wrap<As...>>;
template <typename ...As>
XTAL_DEF_(return,inline)
XTAL_LET wrap_f(auto &&o)
XTAL_0EX -> decltype(auto)
{
	return wrap_t<As...>::function(XTAL_REF_(o));
};


////////////////////////////////////////////////////////////////////////////////

template <typename ...As> requires (1 <= sizeof...(As))
struct wrap<As...>
:	process::lift<wrap<>, bond::compose<As...>>
{};
template <typename ...As> requires (0 == sizeof...(As))
struct wrap<As...>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline)
		XTAL_SET function(auto &&o)
		XTAL_0EX -> auto
		{
			return o - round(o);
		}
		template <auto ...>
		XTAL_DEF_(return,inline)
		XTAL_SET function(algebra::bicycle_q auto &&t_)
		XTAL_0EX -> decltype(auto)
		{
			return XTAL_REF_(t_);
		}
		template <auto ...>
		XTAL_DEF_(return,inline)
		XTAL_SET function(complex_field_q auto &&o)
		XTAL_0EX -> decltype(auto)
		{
			if constexpr (complex_number_q<decltype(o)>) {
				auto &xy = involved_f(o);
				return complexion_f(function(xy[0]), function(xy[1]));
			}
			else {
				return complexion_f(function(o.real()), function(o.imag()));
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
