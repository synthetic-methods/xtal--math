#pragma once
#include "./any.hh"

#include "./identity.hh"




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
private:// OPERATE

	XTAL_DEF_(return,inline,set)
	per_f(real_q auto &&o)
	noexcept -> decltype(auto)
	{
		using  U_fit = bond::fit<decltype(o)>;
		return U_fit::template patio_f<-1>(2)/XTAL_REF_(o);
	}
	XTAL_DEF_(return,inline,set)
	per_f(atom::math::phason_q auto &&o)
	noexcept -> decltype(auto)
	{
		return per_f(XTAL_REF_(o) (1));
	}

public:
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;
		XTAL_TYP_(set) U_stage = XTAL_ALL_(M);
		XTAL_DEF_(set) N_stage = M.head();

	public:
		template <int M_arg=0>
		struct infix
		{
			using superkind = identity_t<>::template infix<M_arg>;

			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
				using R_ = bond::compose_s<R, superkind>;

			public:// CONSTRUCT
				using R_::R_;
				using typename R_::data_type;

			public:// OPERATE

				/*!
				\brief   Generates a gate or `HeavisidePi[# - 1/2] &`.
				*/
				template <auto ...Ns> requires un_v<N_stage, 0>
				XTAL_DEF_(return,inline,let)
				method(auto &&o, auto &&...oo)
				noexcept -> decltype(auto)
				{
					using X = XTAL_ALL_(per_f(o));
					auto const    &u =   R_::template head<U_stage>();
					auto const    &v =   u.head();
					auto const     i = _xtd::make_unsigned_f((signed) v);
					auto constexpr I = _xtd::make_unsigned_f((N_stage));
					auto const     x =  static_cast<X>(i < I);
					return two*R_::template method<Ns...>(half*x, XTAL_REF_(o), XTAL_REF_(oo)...);
				}
				/*!
				\brief   Generates an impulse or `DiracDelta`,
							with magnitude inversely proportional to the provided frequency.
				\param   o Frequency scaling.
				*/
				template <auto ...Ns> requires in_v<N_stage, 0>
				XTAL_DEF_(return,inline,let)
				method(auto &&o, auto &&...oo)
				noexcept -> decltype(auto)
				{
					auto  &u =  R_::template head<U_stage>();
					auto  &v =  u.head();
					signed n = !v;
					auto   x =  n*per_f(o);
					u |= n;
					return R_::template method<Ns...>(x, XTAL_REF_(o), XTAL_REF_(oo)...);
				}

			};
			template <class R> requires incomplete_q<typename R::template head<M>>
			class subtype<R> : public bond::compose_s<R
			,	typename R::stage_type::template attach<>
			,	typename pulse_t<typename R::stage_type{M}>::template infix<M_arg>
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
		struct infix
		{
			template <class R>
			class subtype : public bond::compose_s<R
			,	typename occur::stage_t<>::template attach<>
			,	typename pulse_t<     occur::stage_t<> {M}>::template infix<M_arg>
			>
			{
			};
			template <class R> requires requires {requires occur::stage_q<typename R::stage_type>;}
			class subtype<R> : public bond::compose_s<R
			,	typename pulse_t<typename R::stage_type{M}>::template infix<M_arg>
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
