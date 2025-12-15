#pragma once
#include "./any.hh"

#include "../process/pade/tangy.hh"




XTAL_ENV_(push)
namespace xtal::provision::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_arg=0>	struct  prewarping;
template <class ..._s>	concept prewarping_q = bond::tab_in_p<prewarping<-1>, _s...>;


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Prewarps the `gain` parameter of `method`, indexed from zero by `M_arg`.
*/
template <>
struct prewarping<-1>
:	bond::tab<prewarping<-1>>
{
};
template <int M_arg>
struct prewarping
{
	using superkind = prewarping<-1>;

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
			auto &o = get<M_arg>(_std::tie(oo...));
			o *= S_::template head<occur::resample_t<>>().period();
			o *= process::math::pade::tangy_f<1, -1>(o);
			return S_::template method<Ns...>(XTAL_MOV_(oo)...);
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
