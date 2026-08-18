#pragma once
#include "./any.hh"

#include "./root.hh"
#include "./part.hh"



XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Scales the argument by `M_val`.
*/
template <auto M_val=1>	struct  dilate;
template <auto M_val=1>	using   dilate_t = process::confined_t<dilate<M_val>>;
template <auto M_val=1>
XTAL_VAL_(return,inline,let)
dilate_f(auto &&o)
noexcept -> decltype(auto)
{
	return dilate_t<M_val>{}.method(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <auto M_val>
struct dilate
{
	XTAL_VAL_(return,inline,set)
	after_f(auto &&o)
	noexcept -> decltype(auto)
	{
		using U_fit = bond::fit<decltype(o)>;
		auto constexpr n_val = U_fit::alpha_f(bond::operate_v<M_val>);
		auto constexpr u     =       part_f<unsigned>(n_val);
		auto constexpr v     = (int) part_f<  signed>(n_val);
		return XTAL_REF_(o)*root_f<+v>(u);
	};

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_VAL_(return,inline,let)
		method(auto &&o) const
		noexcept -> auto
		{
			return after_f(XTAL_REF_(o));
		};

		/*!
		\brief   Wraps the super-`function`, unscaling/scaling the domain/codomain by the given `M_val`.
		\note    Must be applied via `{compose,confined}` (etc) rather than `process::{lift,link}`.
		\todo    Incorporate `M_arg`.
		*/
		template <int M_arg=0>
		struct transfix
		{
			static_assert(M_arg == 0);// For now...

			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:
				using R_::R_;

				template <auto ...Ns>
				XTAL_VAL_(return,inline,let)
				method(auto &&...oo) const
				noexcept -> decltype(auto)
				requires
				requires (R_ const &r_)
				{r_.template method  <Ns...>(       XTAL_REF_(oo)...);}
				{return      method_f<Ns...>(*this, XTAL_REF_(oo)...);}

				template <auto ...Ns>
				XTAL_VAL_(return,inline,let)
				method(auto &&...oo)
				noexcept -> decltype(auto)
				requires
				requires (R_       &r_)
				{r_.template method  <Ns...>(       XTAL_REF_(oo)...);}
				{return      method_f<Ns...>(*this, XTAL_REF_(oo)...);}

			private:
				template <auto ...Ns>
				XTAL_VAL_(return,inline,set)
				method_f(auto &&_, auto &&o)
				noexcept -> decltype(auto)
				{
					using U = typename    bond::fit<decltype(o)>::alpha_type;
					U    constexpr n =    bond::operate_v<M_val>;
					U    constexpr u =       part_f<unsigned>(n);
					auto constexpr I = (int) part_f<  signed>(n);
					return qualify_f<R_>(XTAL_REF_(_)).template method<Ns...>(XTAL_REF_(o)*root_f<-I>(u))*root_f<I>(u);
				}

			};
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
