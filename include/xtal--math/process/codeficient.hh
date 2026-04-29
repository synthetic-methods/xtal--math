#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Diagonalizes a complex number.
\todo    Implement positional parameter `M_pos`?

Multiplies the imaginary component by the given factor,
diagonalizing the result by `(1 + I)*#&`.
*/
template <class ...Ms>	struct  codeficient;
template <class ...Ms>	using   codeficient_t = process::confined_t<codeficient<Ms...>>;

////////////////////////////////////////////////////////////////////////////////

template <>
struct codeficient<>
{
	template <auto ...Ns>
	XTAL_DEF_(return,inline,let)
	method(auto &&o)
	const noexcept -> auto
	{
		auto const y = _std::imag(o);
		auto const x = _std::real(o);
		return atom::couple_f(x - y, x + y);
	}

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
					auto const z =  R_ ::template method<Ns...>(XTAL_REF_(oo)...);
					auto const y = _std::imag(z)*XTAL_REF_(o);
					auto const x = _std::real(z);
					return atom::couple_f(x - y, x + y);
				}

			};
		};

	};
};
template <occur::any_q M>
struct codeficient<M>
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
					auto const z =  R_ ::template method<Ns...>(XTAL_REF_(oo)...);
					auto const y = _std::imag(z)*R_::headed();
					auto const x = _std::real(z);
					return atom::couple_f(x - y, x + y);
				}

			};
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
