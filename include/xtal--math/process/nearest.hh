#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Provides the nearest integer.
*/
template <int M_dir=0>
struct nearest;

template <int M_dir=0>
using nearest_t = process::confined_t<nearest<M_dir>>;

template <int M_dir=0>
XTAL_DEF_(return,inline,let)
nearest_f(auto &&...oo)
noexcept -> decltype(auto)
{
	return nearest_t<M_dir>{}.method(XTAL_REF_(oo)...);
};


////////////////////////////////////////////////////////////////////////////////

template <int M_dir>
struct nearest
{
	template <class S>
	class subtype : public S
	{
	protected:

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_base_f(auto &&z)
		noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (M_dir == -1) {return floor(XTAL_REF_(z));}
			XTAL_0IF (M_dir ==  0) {return round(XTAL_REF_(z));}
			XTAL_0IF (M_dir ==  1) {return ceil (XTAL_REF_(z));}
		}

	public:
		using S::S;

		template <auto ...>
		XTAL_DEF_(return,inline,let)
		method(auto &&z)
		const noexcept -> auto
		{
			using Z = XTAL_ALL_(z);
			using Z_fit = bond::fit<Z>;

			XTAL_IF1_(eval) {
				return method_base_f(XTAL_REF_(z));
			}
			XTAL_0IF (real_variable_q<Z>) {
				using Z_alpha = typename Z_fit::alpha_type;
				using Z_delta = typename Z_fit::delta_type;
				auto u = XTAL_REF_(z);
				XTAL_IF0
				XTAL_0IF (M_dir ==  0) {u += Z_fit::  haplo_f(1);}
				XTAL_0IF (M_dir ==  1) {u += Z_fit::dnsilon_f(1);}
				return Z_alpha(Z_delta(XTAL_REF_(z)));
			}
			XTAL_0IF_(else) {
				return method_base_f(XTAL_REF_(z));
			}
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(complex_field_q auto &&z)
		const noexcept -> decltype(auto)
		{
			return complexion_f(method<Ns...>(z.real()), method<Ns...>(z.imag()));
		}
		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(atom::groupoid_q auto &&z)
		const noexcept -> decltype(auto)
		{
			return XTAL_ALL_(z)::template zip_from<[] XTAL_1FN_(call) (subtype{}.template method<Ns...>)>(z);
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
