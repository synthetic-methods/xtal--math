#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
struct magnum
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(_std::unsigned_integral auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			return o;
		}
		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(_std::  signed_integral auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _op = bond::operate<decltype(o)>;
			auto const v = o >> _op::sign.shift;
			return (o^v) - v;
		}
		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(real_number_q auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _op = bond::operate<decltype(o)>;
			return _xtd::copysign(o, _op::alpha_1);
		}
		template <int ...Ns>
		XTAL_DEF_(short,static)
		XTAL_LET function(complex_number_q auto const &o)
		noexcept -> XTAL_ALL_(o)
		{
			using _op = bond::operate<decltype(o)>;
			return {function<Ns...>(o.real()), function<Ns...>(o.imag())};
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int ...Ms>
using    magnum_t = process::confined_t<magnum<Ms...>>;

template <auto ...Ns>
XTAL_DEF_(short)
XTAL_LET magnum_f(auto &&...oo)
noexcept -> decltype(auto)
{
	//\
	return magnum_t<Ms...>::function(XTAL_REF_(oo)...);
	return magnum_t<>::template function<Ns...>(XTAL_REF_(oo)...);
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
