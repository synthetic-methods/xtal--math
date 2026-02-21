#pragma once
#include "./any.hh"

#include "../../process/dot.hh"
#include "../../process/taylor/logarithm.hh"



XTAL_ENV_(push)
namespace xtal::occur::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class ..._s> struct   dash;
template <class ..._s> using    dash_t = confined_t<dash<_s...>, bond::tag<dash>>;
template <class ..._s> concept  dash_q = bond::tag_in_p<dash, _s...>;


////////////////////////////////////////////////////////////////////////////////

template <class U>
struct dash<U>
{
	using superkind = confer<U>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(any_q<S>);
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;
		
		XTAL_DEF_(set) dot_f = [] XTAL_1FN_(call) (process::math::dot_f);
		XTAL_DEF_(set) exp_f = [] XTAL_1FN_(call) (process::math::taylor::logarithm_t<-1>::template method_f<3>);

	public:
		using S_::S_;
		
		/*!
		\brief  	Attaches `T_`, and appends to the arguments of `method` and `method_f`.
		*/
		template <extent_type N_mask=1>
		struct attend
		{
			using superkind = typename T_::template attach<N_mask>;

			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
				using R_ = bond::compose_s<R, superkind>;
				using Q_ = typename R_::state_type;

			public:// CONSTRUCT
				using R_::R_;

			public:// OPERATE

				template <auto ...Ns>
				XTAL_DEF_(return,inline,let)
				method(auto &&...oo)// NOTE: Non-`const` only, for now...
				noexcept -> decltype(auto)
				{
					auto const &[s_state_] = R_::template memory<Q_>();
					auto const & s_shape_  = R_::template head<T_>().head();
					return R_::template method<Ns...>(XTAL_REF_(oo)..., exp_f(dot_f(s_state_, s_shape_)));
				}

			};
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
