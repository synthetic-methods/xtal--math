#pragma once
#include "./any.hh"
#include "./square.hh"





XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Wraps the super-`function`, applying domain/codomain transformations \
associated with partial implementations (cf. `Sin` vs `Sinc`). \

///\note\
Because it invokes the super-`function` directly, \
it must be applied via `{compose,confined}` (etc) rather than `process::{lift,link}`.

template <int M_pow=1, int M_car=0>
	requires inclusive_q<M_pow, 1, -1> and inclusive_q<M_car, 0, 1, 2>
XTAL_TYP discarding;


////////////////////////////////////////////////////////////////////////////////

template <int M_pow>
struct discarding<M_pow, +0>
{
	template <class S>
	using subtype = bond::compose_s<S, bond::tag<process::link>>;

};
template <int M_pow>
struct discarding<M_pow, +1>
{
	using subkind = bond::tag<process::link>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline,static)
		XTAL_RET function(auto &&u)
		XTAL_0EX
		{
			auto  v = S_::template function<Is...>(u);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);
			static_assert(is_q<U, V>);

			XTAL_IF0
			XTAL_0IF (M_pow ==  1) {return v*XTAL_REF_(u);}
			XTAL_0IF (M_pow == -1) {return v/XTAL_REF_(u);}
		};

	};
};
template <>
struct discarding<1, +1>
{
	using subkind = bond::tag<process::link>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline,static)
		XTAL_RET function(auto &&u)
		XTAL_0EX
		{
			auto  v = S_::template function<Is...>(u);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);

			XTAL_IF0
			XTAL_0IF (is_q<U, V>) {
				return v*XTAL_REF_(u);
			}
			XTAL_0IF (complex_number_q<V>) {
				devolved_f(v)[1] *= XTAL_REF_(u); return v;
			}
			XTAL_0IF (complex_field_q<V>) {
				return complexion_f(v.real(), v.imag()*XTAL_REF_(u));
			}
		};

	};
};
template <>
struct discarding<1, +2>
{
	using subkind = bond::tag<process::link>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline,static)
		XTAL_RET function(auto &&o)
		XTAL_0EX
		{
			return S_::template function<Is...>(square_f<1>(XTAL_REF_(o)));
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
