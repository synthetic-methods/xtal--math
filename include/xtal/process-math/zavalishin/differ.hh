#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Differentiates the incoming signal. \
When the modulating signal and/or phasor are supplied, \
and/or if the value is `std::complex`, \
the result is renormalized w.r.t. the chain-rule of differentiation. \

///\note\
The `method` uses local arena-like allocation to manage state, so `sizeof(input) <= 56`. \

///\todo\
Define a `brace`d version with static-state, \
possibly invoking the parent with the same `template method` parameters (e.g. `N_ord`er). \

template <typename ...As> XTAL_TYP differ;
template <typename ...As> XTAL_USE differ_t = process::confined_t<differ<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
struct differ
:	process::link<differ<>, As...>
{
};
template <>
struct differ<>
{
	using subkind = bond::compose<void
	,	resource::cached<>
	,	resource::example<>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:// CONSTRUCT
		using S_::S_;
		
	public:// FUNC*

		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&u)
		XTAL_0EX -> auto
		{
			using _op = bond::operate<decltype(u)>;
			
			XTAL_USE U = XTAL_ALL_(u);
			auto [u1] = S_::template cache<U>(); U u0 = XTAL_REF_(u);

			_std::swap(u0, u1);
			return u1 - u0;
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&u, algebra::d_::circular_q auto &&z_)
		XTAL_0EX -> auto
		{
			using _op = bond::operate<decltype(u)>;
			
			XTAL_USE U = XTAL_ALL_(u);
			XTAL_USE Z = debraced_t<XTAL_ALL_(z_)>;
			
			auto [u1] = S_::template cache<U>();
			U u0 = XTAL_REF_(u);
			
		//	Divides the difference by the derived slope of the phasor:
			_std::swap(u1, u0); U u10 = u1 - u0; u10 *= root_f<-1, 0>(z_(1));

		//	Resets the state to zero if a phase/frequency discontinuity is detected:
		//	u10 *= abs(z10 - z_(1)) < _op::haplo_f(N_zap);
			return u10;
		}
		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_LET method(auto &&u, auto &&v, algebra::d_::circular_q auto &&z_)
		XTAL_0EX -> auto
		{
			using _op = bond::operate<decltype(u)>;
			
			XTAL_USE U = XTAL_ALL_(u);
			XTAL_USE V = XTAL_ALL_(v);
			
			auto [u1, v1] = S_::template cache<U, V>();
			U u0 = XTAL_REF_(u);
			V v0 = XTAL_REF_(v);
			
		//	Divides the difference by the derived slope of the phasor:
			_std::swap(v1, v0); V v10 = v1 - v0; v10 += z_(1);
			_std::swap(u1, u0); U u10 = u1 - u0; u10 *= root_f<-1, 0>(v10);

			return u10;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
