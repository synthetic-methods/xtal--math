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

template <typename ...As> struct   differ;
template <typename ...As> using    differ_t = process::confined_t<differ<As...>>;


////////////////////////////////////////////////////////////////////////////////

template <typename ...As>
struct differ
:	process::link<differ<>, As...>
{
};
template <>
struct differ<>
{
	using alpha_type = typename bond::operating::alpha_type;
	using aphex_type = typename bond::operating::aphex_type;

	using superkind = bond::compose<void
	,	provision::example<>
	>;
	using subkind = bond::compose<void
	,	differ<>
	,	provision::cached<aphex_type[2]>
	>;
	template <class S>
	class subtype : public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:// CONSTRUCT
		using S_::S_;

	};		
	template <class S> requires provision::cached_q<S>
	class subtype<S> : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;
		
	public:// FUNC*

		template <auto ...Is>
		XTAL_DEF_(short)
		XTAL_LET method(auto &&u)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(u)>;
			
			using   U = XTAL_ALL_(u);
			auto [u1] = S_::template cache<U>(); U u0 = XTAL_REF_(u);

			_std::swap(u0, u1);
			return u1 - u0;
		}
		template <auto ...Is>
		XTAL_DEF_(short)
		XTAL_LET method(auto &&u, algebra::phason_q auto &&t_)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(u)>;
			using   U = XTAL_ALL_(u);
			
			auto [u1] = S_::template cache<U>();
			U u0 = XTAL_REF_(u);
			
		//	Divides the difference by the derived slope:
			_std::swap(u1, u0); U u10 = u1 - u0; u10 *= root_f<-1, 0>(t_(1));

		//	Resets the state to zero if a phase/frequency discontinuity is detected:
		//	u10 *= abs(z10 - t_(1)) < _op::haplo_f(N_zap);
			return u10;
		}
		template <auto ...Is>
		XTAL_DEF_(short)
		XTAL_LET method(auto &&u, auto &&v, algebra::phason_q auto &&t_)
		noexcept -> auto
		{
			using _op = bond::operate<decltype(u), decltype(v)>;
			using   U = XTAL_ALL_(u);
			using   V = XTAL_ALL_(v);
			
			auto [u1, v1] = S_::template cache<U, V>();
			U u0 = XTAL_REF_(u);
			V v0 = XTAL_REF_(v);
			
		//	Divides the difference by the derived slope:
			_std::swap(v1, v0); V v10 = v1 - v0; v10 += t_(1);
			_std::swap(u1, u0); U u10 = u1 - u0; u10 *= root_f<-1>(v10);

			return u10;
		}

	public:// BIND*
	//	TODO: Move this to `processor::math`?

		template <class ...Xs>
		struct bind
		{
			using superkind = bond::compose<void
			,	subkind
			,	provision::cached<bond::seek_front_t<Xs...>[2]>
			>;
			template <class R>
			using subtype = bond::compose_s<R, superkind>;

		};
		template <class ...Xs> requires iterated_q<Xs...>
		struct bind<Xs...>
		{
			using superkind = bond::compose<void
			,	subkind
			,	provision::cached<iteratee_t<bond::seek_front_t<Xs...>[2]>>
			>;
			template <class R>
			using subtype = bond::compose_s<R, superkind>;

		};

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
