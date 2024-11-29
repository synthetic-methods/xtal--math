#pragma once
#include "./any.hh"
#include "./term.hh"
#include "../square.hh"




XTAL_ENV_(push)
namespace xtal::process::math::horner
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Evaluates the term `(x^2 + xs^2...)` (using fused multiply-add, if supported by the compiler). \

template <int N_sign=1> XTAL_TYP square;
template <int N_sign=1> XTAL_USE square_t = process::confined_t<square<N_sign>>;
template <int N_sign=1>
XTAL_DEF_(return,inline)
XTAL_LET square_f(auto const &x, auto &&...xs)
XTAL_0EX -> auto
{
	static_assert(N_sign == 1);// For now...
	if constexpr (0 == sizeof...(xs)) {
		return math::square_f(x);
	}
	else {
		return term_f(square_f(XTAL_REF_(xs)...), x, x);
	}
};


////////////////////////////////////////////////////////////////////////////////

template <int N_sign>
struct square
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline)
		XTAL_SET function(auto &&...oo)
		XTAL_0EX -> decltype(auto)
		{
			return square_f<N_sign>(XTAL_REF_(oo)...);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
