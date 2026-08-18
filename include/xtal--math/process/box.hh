#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Wraps the super-`method`.
*/
template <class ...Ms>	struct  box;
template <class ...Ms>	using   box_t = process::confined_t<box<Ms...>>;


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Wraps the super-`method` return-value with `M_f`.
*/
template <applicable_q M_f>
struct box<M_f>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...Ys>
		struct prepend
		{
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
				method_f(auto &&_, auto &&...oo)
				noexcept -> auto
				{
					auto  &r_ = qualify_f<R_>(XTAL_REF_(_));
					using  Y  = XTAL_ALL_(r_.template method<Ns...>(XTAL_REF_(oo)...));
					return bond::operate<M_f>{}(Y(Ys)..., r_.template method<Ns...>(XTAL_REF_(oo)...));
				}

			};
		};

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
