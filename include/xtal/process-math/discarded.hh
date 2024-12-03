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
	requires in_n<M_pow, 1, -1> and in_n<M_car, 0, 1, 2>
XTAL_TYP discarded;


////////////////////////////////////////////////////////////////////////////////

template <int M_pow>
struct discarded<M_pow, +0>
{
	template <class S>
	using subtype = bond::compose_s<S>;

};
template <int M_pow>
struct discarded<M_pow, +1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	//	TODO: Account for `const &` and `&`, or find a cleaner way to express...
		XTAL_DO2_(template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&u, auto &&...oo), -> decltype(auto)
		{
			auto  v = S_::template method<Is...>(u, XTAL_REF_(oo)...);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);
			static_assert(is_q<U, V>);

			XTAL_IF0
			XTAL_0IF (M_pow ==  1) {return v*XTAL_REF_(u);}
			XTAL_0IF (M_pow == -1) {return v/XTAL_REF_(u);}
		})
		template <auto ...Is>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&u, auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto  v = S_::template function<Is...>(u, XTAL_REF_(oo)...);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);
			static_assert(is_q<U, V>);

			XTAL_IF0
			XTAL_0IF (M_pow ==  1) {return v*XTAL_REF_(u);}
			XTAL_0IF (M_pow == -1) {return v/XTAL_REF_(u);}
		}

	};
};
template <>
struct discarded<1, +1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		XTAL_DO2_(template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&u, auto &&...oo), -> decltype(auto)
		{
			auto  v = S_::template method<Is...>(u, XTAL_REF_(oo)...);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);

			XTAL_IF0
			XTAL_0IF (is_q<U, V>) {
				return v*XTAL_REF_(u);
			}
			XTAL_0IF (complex_number_q<V>) {
				devalued_f(v)[1] *= XTAL_REF_(u); return v;
			}
			XTAL_0IF (complex_field_q<V>) {
				return complexion_f(v.real(), v.imag()*XTAL_REF_(u));
			}
		})
		template <auto ...Is>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&u, auto &&...oo)
		noexcept -> decltype(auto)
		{
			auto  v = S_::template function<Is...>(u, XTAL_REF_(oo)...);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);

			XTAL_IF0
			XTAL_0IF (is_q<U, V>) {
				return v*XTAL_REF_(u);
			}
			XTAL_0IF (complex_number_q<V>) {
				devalued_f(v)[1] *= XTAL_REF_(u); return v;
			}
			XTAL_0IF (complex_field_q<V>) {
				return complexion_f(v.real(), v.imag()*XTAL_REF_(u));
			}
		}

	};
};
template <>
struct discarded<1, +2>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		XTAL_DO2_(template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&u, auto &&...oo), -> decltype(auto)
		{
			return S_::template method<Is...>(square_f(XTAL_REF_(u)), XTAL_REF_(oo)...);
		})
		template <auto ...Is>
		XTAL_DEF_(return,inline,static)
		XTAL_LET function(auto &&u, auto &&...oo)
		noexcept -> decltype(auto)
		{
			return S_::template function<Is...>(square_f(XTAL_REF_(u)), XTAL_REF_(oo)...);
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
