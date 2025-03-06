#pragma once
#include "./any.hh"

#include "../pade/tangy.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class ...As>	struct  prewarped;
template <class ...As>	using   prewarped_t = confined_t<prewarped<As...>>;
template <class ..._s>	concept prewarped_q = bond::tag_in_p<prewarped, _s...>;


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Prewarps the `gain` parameter of `method`, indexed from zero by the constant `I`.
*/
template <constant_q I, typename ...As>
struct prewarped<I, As...>
{
	using superkind = bond::compose<As..., bond::tag<prewarped>>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:
		using S_::S_;

		template <auto ...Ns>
		XTAL_DEF_(return,inline,let)
		method(auto ...oo)
		noexcept -> decltype(auto)
		{
			if constexpr (0 <= I{}()) {
				using resample_type = occur::resample_t<>;
				auto &o = get<I{}>(_std::tie(oo...));
				o *= S_::template head<resample_type>().period();
				o *= pade::tangy_t<1,-1>::template method_f<6>(o);
				return S_::template method<Ns...>(XTAL_MOV_(oo)...);
			}
		};

	};
};
template <variable_q I, typename ...As>
struct prewarped<I, As...>
:	prewarped<ordinal_constant_t<0>, As...>
{
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
