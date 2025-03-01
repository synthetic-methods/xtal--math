#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ...As>	struct  wrap;
template <typename ...As>	using   wrap_t = process::confined_t<wrap<As...>>;
template <typename ...As>
XTAL_DEF_(return,inline,let)
wrap_f(auto &&o)
noexcept -> decltype(auto)
{
	return wrap_t<As...>::method_f(XTAL_REF_(o));
};


////////////////////////////////////////////////////////////////////////////////

template <typename ...As> requires some_q<As...>
struct wrap<As...>
:	process::lift<wrap<>, bond::compose<As...>>
{};
template <typename ...As> requires none_q<As...>
struct wrap<As...>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> auto
		{
			return bond::math::bit_fraction_f(o);
		}
		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(atom::math::phason_q auto &&t_)
		noexcept -> decltype(auto)
		{
			return XTAL_REF_(t_);
		}
		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(complex_field_q auto &&o)
		noexcept -> decltype(auto)
		{
			if constexpr (complex_variable_q<decltype(o)>) {
				auto &xy = destruct_f(o);
				return complexion_f(method_f(xy[0]), method_f(xy[1]));
			}
			else {
				return complexion_f(method_f(o.real()), method_f(o.imag()));
			}
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
