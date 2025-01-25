#pragma once
#include "./any.hh"

#include "./root.hh"
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

template <int M_car=0, int M_aux=0> requires in_n<M_car, 0, 1, 2>
struct   discarded;


////////////////////////////////////////////////////////////////////////////////

template <int M_aux>
struct discarded<0, M_aux>
{
	template <class S>
	using subtype = bond::compose_s<S>;

};
template <int M_aux>
struct discarded<1, M_aux>
{
	XTAL_SET M_pow = sign_n<M_aux, 1>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	//	TODO: Account for `const &` and `&`, or find a cleaner way to express...
		XTAL_DO2_(template <auto ...Is>
		XTAL_DEF_(short)
		XTAL_LET        method(auto &&u, auto &&...oo),
		noexcept -> auto
		{
			auto  v = S_::template        method<Is...>(u, XTAL_REF_(oo)...);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);
			static_assert(same_q<U, V>);

			return v*root_f<M_pow, 1>(XTAL_REF_(u));
		})
		XTAL_DO0_(template <auto ...Is>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto &&u, auto &&...oo),
		noexcept -> auto
		{
			auto  v = S_::template static_method<Is...>(u, XTAL_REF_(oo)...);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);
			static_assert(same_q<U, V>);

			return v*root_f<M_pow, 1>(XTAL_REF_(u));
		})

	};
};
template <>
struct discarded<1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		XTAL_DO2_(template <auto ...Is>
		XTAL_DEF_(short)
		XTAL_LET        method(auto &&u, auto &&...oo),
		noexcept -> auto
		{
			auto  v = S_::template        method<Is...>(u, XTAL_REF_(oo)...);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);

			XTAL_IF0
			XTAL_0IF (same_q<U, V>) {
				return v*XTAL_REF_(u);
			}
			XTAL_0IF (complex_variable_q<V>) {
				destruct_f(v)[1] *= XTAL_REF_(u); return v;
			}
			XTAL_0IF (complex_field_q<V>) {
				return complexion_f(v.real(), v.imag()*XTAL_REF_(u));
			}
		})
		XTAL_DO0_(template <auto ...Is>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto &&u, auto &&...oo),
		noexcept -> auto
		{
			auto  v = S_::template static_method<Is...>(u, XTAL_REF_(oo)...);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);

			XTAL_IF0
			XTAL_0IF (same_q<U, V>) {
				return v*XTAL_REF_(u);
			}
			XTAL_0IF (complex_variable_q<V>) {
				destruct_f(v)[1] *= XTAL_REF_(u); return v;
			}
			XTAL_0IF (complex_field_q<V>) {
				return complexion_f(v.real(), v.imag()*XTAL_REF_(u));
			}
		})

	};
};
template <int M_aux>
struct discarded<2, M_aux>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		XTAL_DO2_(template <auto ...Is>
		XTAL_DEF_(short)
		XTAL_LET        method(auto &&u, auto &&...oo),
		noexcept -> decltype(auto)
		{
			using _fix = bond::fixture<decltype(u)>;
			XTAL_LET v = _fix::alpha_f(sign_n<(M_aux&1)^1, -1>);
			return S_::template        method<Is...>(v*square_f(XTAL_REF_(u)), XTAL_REF_(oo)...);
		})
		XTAL_DO0_(template <auto ...Is>
		XTAL_DEF_(short,static)
		XTAL_LET static_method(auto &&u, auto &&...oo),
		noexcept -> decltype(auto)
		{
			using _fix = bond::fixture<decltype(u)>;
			XTAL_LET v = _fix::alpha_f(sign_n<(M_aux&1)^1, -1>);
			return S_::template static_method<Is...>(v*square_f(XTAL_REF_(u)), XTAL_REF_(oo)...);
		})

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
