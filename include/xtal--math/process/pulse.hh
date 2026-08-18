#pragma once
#include "./any.hh"

#include "./rate.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Prepends a discrete/non-bandlimited control signal to the `method` arguments.
\tparam  M Specifies the final `stage`.
*/
template <auto  ..._s> XTAL_TYP_(new) pulse;
template <auto  ..._s> XTAL_TYP_(set) pulse_t = confined_t<pulse<_s...>, bond::tab<pulse<>>>;
template <class ..._s> XTAL_TYP_(ask) pulse_q = bond::tab_inner_p<pulse<>, _s...>;


////////////////////////////////////////////////////////////////////////////////

template <occur::stage_q auto M>
struct pulse<M>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;
		XTAL_TYP_(set) U_stage = XTAL_ALL_(M);
		XTAL_VAL_(set) N_stage = M.head();

	public:
		template <int M_arg=0>
		struct transfix
		{
			static_assert(M_arg == 0);// For now...

			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:// CONSTRUCT
				using R_::R_;

			public:// OPERATE

				/*!
				\brief   Scales the input `x` with a gate or `HeavisidePi[# - 1/2] &`.
				*/
				template <auto ...Ns> requires (N_stage != 0)
				XTAL_VAL_(return,inline,let)
				method(auto x, auto &&...oo) const
				noexcept -> decltype(auto)
				{
					auto const     &a = R_::self().template head<U_stage>();
					auto const      i = xtd::unsigned_cast(a.head());
					auto constexpr  I = xtd::unsigned_cast(N_stage);
					x *= static_cast<XTAL_ALL_(x)>(i < I);
					return R_::template method<Ns...>(x, XTAL_REF_(oo)...);
				}
				/*!
				\brief   Scales the input `x` with an impulse or `DiracDelta`,
							with magnitude inversely proportional to the provided frequency.
				*/
				template <auto ...Ns> requires (N_stage == 0)
				XTAL_VAL_(return,inline,let)
				method(auto x, auto &&...oo) const
				noexcept -> decltype(auto)
				{
					auto           &a = const_cast<U_stage &>(R_::template head<U_stage>());
					auto const      n = (signed) !a.head();
				//	x *= rate_f<-1>(u, bond::fit<X>::patio_f(2));
					x *= n;
					a |= n;
					return R_::template method<Ns...>(x, XTAL_REF_(oo)...);
				}

			};
			template <class R> requires incomplete_q<typename R::template head<M>>
			class subtype<R> : public bond::compose_s<R
			,	typename R::stage_type::template attach<>
			,	typename pulse_t<typename R::stage_type(M)>::template transfix<M_arg>
			>
			{
			};
		};

	};
};
template <ordinal_q auto M>
struct pulse<M>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		template <int M_arg=0>
		struct transfix
		{
			template <class R>
			class subtype : public bond::compose_s<R
			,	typename occur::stage_t<>::template attach<>
			,	typename pulse_t<     occur::stage_t<> (M)>::template transfix<M_arg>
			>
			{
			};
			template <class R> requires requires {requires occur::stage_q<typename R::stage_type>;}
			class subtype<R> : public bond::compose_s<R
			,	typename pulse_t<typename R::stage_type(M)>::template transfix<M_arg>
			>
			{
			};
		};

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
