#pragma once
#include "./any.hh"

#include "./root.hh"
#include "./square.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Wraps the super-`function`, applying domain/codomain transformations
associated with partial implementations (cf. `Sin` vs `Sinc`).

\note    Because it invokes the super-`function` directly,
it must be applied via `{compose,confined}` (etc) rather than `process::{lift,link}`.
*/
template <int M_car=0, int M_aux=0> requires in_v<M_car, 0, 1, 2>
struct  discarded;


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
	static int constexpr M_pow = sign_v<M_aux, 1>;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

	//	TODO: Account for `const &` and `&`, or find a cleaner way to express...
		XTAL_FN2_(do) (template <auto ...Is>
		XTAL_DEF_(return,inline,let)
		method(auto &&u, auto &&...oo),
		noexcept -> auto
		{
			auto  v = S_::template method<Is...>(u, XTAL_REF_(oo)...);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);
			static_assert(same_q<U, V>);

			return v*root_f<M_pow, 1>(XTAL_REF_(u));
		})
		template <auto ...Is>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&u, auto &&...oo)
		noexcept -> auto
		{
			auto  v = S_::template method_f<Is...>(u, XTAL_REF_(oo)...);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);
			static_assert(same_q<U, V>);

			return v*root_f<M_pow, 1>(XTAL_REF_(u));
		}

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

		XTAL_FN2_(do) (template <auto ...Is>
		XTAL_DEF_(return,inline,let)
		method(auto &&u, auto &&...oo),
		noexcept -> auto
		{
			auto  v = S_::template method<Is...>(u, XTAL_REF_(oo)...);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);

			XTAL_IF0
			XTAL_0IF (same_q<U, V>) {
				return v*XTAL_REF_(u);
			}
			XTAL_0IF (complex_variable_q<V>) {
				XTAL_IF1_(consteval) {
					return V{v.real(), v.imag()*XTAL_REF_(u)};
				}
				XTAL_0IF_(else) {
					destruct_f(v)[1] *= XTAL_REF_(u); return v;
				}
			}
			XTAL_0IF (complex_field_q<V>) {
				return complexion_f(v.real(), v.imag()*XTAL_REF_(u));
			}
		})
		template <auto ...Is>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&u, auto &&...oo)
		noexcept -> auto
		{
			auto  v = S_::template method_f<Is...>(u, XTAL_REF_(oo)...);
			using V = XTAL_ALL_(v);
			using U = XTAL_ALL_(u);

			XTAL_IF0
			XTAL_0IF (same_q<U, V>) {
				return v*XTAL_REF_(u);
			}
			XTAL_0IF (complex_variable_q<V>) {
				XTAL_IF1_(consteval) {
					return V{v.real(), v.imag()*XTAL_REF_(u)};
				}
				XTAL_0IF_(else) {
					destruct_f(v)[1] *= XTAL_REF_(u); return v;
				}
			}
			XTAL_0IF (complex_field_q<V>) {
				return complexion_f(v.real(), v.imag()*XTAL_REF_(u));
			}
		}

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

		XTAL_FN2_(do) (template <auto ...Is>
		XTAL_DEF_(return,inline,let)
		method(auto &&u, auto &&...oo),
		noexcept -> decltype(auto)
		{
			using _fit = bond::fit<decltype(u)>;
			auto constexpr v = _fit::alpha_f(cosign_v<M_aux>);
			return S_::template        method<Is...>(v*square_f(XTAL_REF_(u)), XTAL_REF_(oo)...);
		})
		template <auto ...Is>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&u, auto &&...oo)
		noexcept -> decltype(auto)
		{
			using _fit = bond::fit<decltype(u)>;
			auto constexpr v = _fit::alpha_f(cosign_v<M_aux>);
			return S_::template method_f<Is...>(v*square_f(XTAL_REF_(u)), XTAL_REF_(oo)...);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
