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
XTAL_DEF_(return,inline,let)
identity_f(auto &&x, auto &&...oo)
noexcept -> decltype(auto)
{
	return identity_t<M_ism, M_car>::method_f(XTAL_REF_(x), XTAL_REF_(oo)...);
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
		XTAL_DEF_(return,inline,set)
		method_f(auto &&x, auto &&...oo)
		noexcept -> auto
		{
			XTAL_IF0
			XTAL_0IF (M_car >= -0) {return XTAL_REF_(x);}
			XTAL_0IF (M_car == -1) {return XTAL_ALL_(x) {one};}
			XTAL_0IF (M_car == -2) {return XTAL_ALL_(x) {zero};}
		}

	public:// CLASS

		template <int M_arg=0>
		struct infix
		{
			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:
				using R_::R_;

			};
			template <class R> requires in_v<M_arg, 1>
			class subtype<R> : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:
				using R_::R_;

				template <auto ...Ns>
				XTAL_DEF_(return,inline,set)
				method_f(auto &&o, auto &&...oo)
				noexcept -> decltype(auto)
				requires      in_v<requires {R ::template method_f<Ns...>  (XTAL_REF_(oo)..., XTAL_REF_(o));}>
				{
					return                    R_::template method_f<Ns...>  (XTAL_REF_(oo)..., XTAL_REF_(o));
				};

				template <auto ...Ns>
				XTAL_DEF_(return,inline,let)
				method(auto &&o, auto &&...oo) const
				noexcept -> decltype(auto)
				requires      un_v<requires {R ::template method_f<Ns...>  (XTAL_REF_(oo)..., XTAL_REF_(o));}>
				and requires (R_ const &s_) {s_ .template method  <Ns...>  (XTAL_REF_(oo)..., XTAL_REF_(o));}
				{
					return                    R_::template method  <Ns...>  (XTAL_REF_(oo)..., XTAL_REF_(o));
				};
				template <auto ...Ns>
				XTAL_DEF_(return,inline,let)
				method(auto &&o, auto &&...oo)
				noexcept -> decltype(auto)
				requires      un_v<requires {R ::template method_f<Ns...>  (XTAL_REF_(oo)..., XTAL_REF_(o));}>
				and requires (R_       &s_) {s_ .template method  <Ns...>  (XTAL_REF_(oo)..., XTAL_REF_(o));}
				{
					return                    R_::template method  <Ns...>  (XTAL_REF_(oo)..., XTAL_REF_(o));
				};

			};
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
