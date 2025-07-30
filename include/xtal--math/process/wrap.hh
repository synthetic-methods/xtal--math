#pragma once
#include "./any.hh"

#include "./nearest.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <typename ...As>	struct  wrap;
template <typename ...As>	using   wrap_t = process::confined_t<wrap<As...>>;
template <typename ...As>
XTAL_DEF_(return,inline,let)
wrap_f(auto &&z)
noexcept -> decltype(auto)
{
	return wrap_t<As...>::method_f(XTAL_REF_(z));
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

		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&z)
		noexcept -> auto
		{
			using Z = XTAL_ALL_(z);

			XTAL_IF1_(to) (bond::math::bit_fraction_f(XTAL_REF_(z)))
			XTAL_0IF_(to) (z - nearest_f(z))
			XTAL_0IF_(terminate)
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,set)
		method_f(atom::math::phason_q auto &&z_)
		noexcept -> decltype(auto)
		{
			return XTAL_REF_(z_);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
