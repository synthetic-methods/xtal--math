#pragma once
#include "./any.hh"

#include "./dot.hh"
#include "./term.hh"
#include "./square.hh"


XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Multiplies the leading argument with the result of the parent-`method` applied to the trailing arguments.
\todo    Implement positional parameter `M_pos`?
*/
template <class ...Ms>	struct  coefficient;
template <class ...Ms>	using   coefficient_t = process::confined_t<coefficient<Ms...>>;

////////////////////////////////////////////////////////////////////////////////

template <occur::any_q M>
struct coefficient<M>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <extent_type M_mask=1>
		struct attach
		{
			using superkind = typename M::template attach<M_mask>;

			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
				using R_ = bond::compose_s<R, superkind>;

			public:
				using R_::R_;

				template <auto ...Ns>
				XTAL_DEF_(inline,let)
				method(auto &&...oo)
				const noexcept -> auto
				{
					return R_::headed()*R_::template method<Ns...>(XTAL_REF_(oo)...);
				}

			};
		};

	};
};
template <>
struct coefficient<>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int M_arg=0>
		struct unfix
		{
			static_assert(M_arg == 0);

			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:
				using R_::R_;

				template <auto ...Ns>
				XTAL_DEF_(return,inline,let)
				method(auto &&o, auto &&...oo)
				const noexcept -> auto
				{
					return R_::template method<Ns...>(XTAL_REF_(oo)...)*XTAL_REF_(o);
				}

			};
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
