#pragma once
#include "./any.hh"

#include "../../process/dot.hh"
#include "../../process/taylor/logarithm.hh"



XTAL_ENV_(push)
namespace xtal::occur::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Appends the exponentiated-dot-product of the stored value with the internal `state`.
\todo    Consider working with warp and damp in the logarithmic/additive domain.
\todo    Combine with existing argument if possible,
         or include a constant parameter (or the highpass component `1 - s[0] - s[1] - ...`).
*/
template <class ..._s> XTAL_TYP_(new) probe;
template <class ..._s> XTAL_TYP_(set) probe_t = confined_t<probe<_s...>, bond::tag<probe>>;
template <class ..._s> XTAL_TYP_(ask) probe_q = bond::tag_inner_p<probe, _s...>;


////////////////////////////////////////////////////////////////////////////////

template <class U>
struct probe<U>
{
	using superkind = confer<U>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		static_assert(any_q<S>);
		using S_ = bond::compose_s<S, superkind>;
		using T_ = typename S_::self_type;
		using U_ = U;
		
		XTAL_DEF_(set) dot_f = [] XTAL_1FN_(call) (process::math::dot_f);
		XTAL_DEF_(set) exp_f = [] XTAL_1FN_(call) (process::math::taylor::logarithm_t<-1>{}.template method<3>);

	public:
		using S_::S_;
		
		template <extent_type N_mask=1>
		struct affix
		{
			using superkind = typename S_::template attach<N_mask>;

			template <class R>
			class subtype : public bond::compose_s<R, superkind>
			{
				using R_ = bond::compose_s<R, superkind>;

			public:// CONSTRUCT
				using R_::R_;

			public:// OPERATE

				template <auto ...Ns>
				XTAL_DEF_(return,inline,let)
				method(auto &&...oo)
				const noexcept -> decltype(auto)
				{
					auto const &[s_state_] = R_::template memory<typename R_::data_type>();
					auto const & s_shape_  = R_::template head<T_>().template head<U_>();
					return R_::template method<Ns...>(XTAL_REF_(oo)..., exp_f(dot_f(s_state_, s_shape_)));
				}

			};
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
