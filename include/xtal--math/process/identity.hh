#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_ism=0, int M_car=0>
struct identity;

template <int M_ism=0, int M_car=0>
using identity_t = process::confined_t<identity<M_ism, M_car>>;

template <int M_ism=0, int M_car=0>
XTAL_VAL_(return,inline,let)
identity_f(auto &&x, auto &&...oo)
noexcept -> decltype(auto)
{
	return identity_t<M_ism, M_car>{}.method(XTAL_REF_(x), XTAL_REF_(oo)...);
};


////////////////////////////////////////////////////////////////////////////////

template <int M_ism, int M_car>
struct identity
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE
		template <auto ...>
		XTAL_VAL_(return,inline,let)
		method(auto &&x, auto &&...oo)
		const noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (M_car >= -0) {return XTAL_REF_(x);}
			XTAL_0IF (M_car == -1) {return XTAL_ALL_(x) {one};}
			XTAL_0IF (M_car == -2) {return XTAL_ALL_(x) {zero};}
		}

	public:// CLASS

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
