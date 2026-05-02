#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Divides the sample indexed by `M_arg` by
         the attached `M_att`.
*/
template <class ..._s> XTAL_TYP_(new) per;
template <class ..._s> XTAL_TYP_(set) per_t = confined_t<per<_s...>, bond::tab<per<>>>;
template <class ..._s> XTAL_TYP_(ask) per_q = bond::tab_inner_p<per<>, _s...>;


////////////////////////////////////////////////////////////////////////////////

template <class M_att>
struct per<M_att>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		template <int M_arg=0>
		struct refix
		{
			template <class R>
			class subtype : public bond::compose_s<R>
			{
				using R_ = bond::compose_s<R>;

			public:// CONSTRUCT
				using R_::R_;
				using typename R_::data_type;

			public:// OPERATE

				template <auto ...Ns> requires un_v<M_arg, 0>
				XTAL_DEF_(return,inline,let)
				method(auto ...oo)
				const noexcept -> decltype(auto)
				{
					auto &o = get<M_arg>(_std::tie(oo...)); o *= period();
					return R_::template method<Ns...>(XTAL_MOV_(oo)...);
				}
				template <auto ...Ns> requires in_v<M_arg, 0>
				XTAL_DEF_(return,inline,let)
				method(auto &&o, auto &&...oo)
				const noexcept -> decltype(auto)
				{
					return R_::template method<Ns...>(XTAL_REF_(o)*period(), XTAL_REF_(oo)...);
				}

			private:// OPERATE

				XTAL_DEF_(return,inline,let)
				period()
				const noexcept -> decltype(auto)
				{
					auto const &m = R_::template head<M_att>();
					XTAL_IF0
					XTAL_0IF (occur::resample_q<M_att>) {return     m.period();}
					XTAL_0IF_(else)                     {return one/m.  head();}
				}

			};
			template <class R> requires incomplete_q<typename R::template head_t<M_att>>
			class subtype<R> : public bond::compose_s<R
			,	typename M_att::template attach<>
			,	refix<M_arg>
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
