#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::atom::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Extends `couple` with `operator*` redefined by the scalar product. \
Indended to act as a coefficient of a similar type where a scalar result is required. \

///\note\
May get refactored as either rational addition (numerator), \
or complex multiplication (imaginary component). \

///\todo\
Specialize `plus_multiplies` or `fma`? \

template <class ..._s>	struct  fadon;
template <class ..._s>	using   fadon_t = typename fadon<_s...>::type;
template <class ..._s>	concept fadon_q = bond::array_or_any_tags_p<fadon_t, _s...> and fixed_shaped_q<_s...>;

XTAL_DEF_(let) fadon_f = [] XTAL_1FN_(call) (_detail::fake_f<fadon_t>);


////////////////////////////////////////////////////////////////////////////////

template <scalar_q ..._s> requires same_q<_s...>
struct fadon<_s ...>
:	fadon<common_t<_s...>[sizeof...(_s)]>
{
};
template <class ..._s>
struct fadon
{
	template <class T>
	using endotype = typename couple<_s...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<fadon_t>>;

	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;

	public:// CONSTRUCT
		using S_::S_;

	public:// OPERATE

		XTAL_DEF_(return,inline,let)
		operator * (auto const &t) const
		noexcept -> auto
		{
			return S_::product(reinterpret_cast<T const &>(t));
		}

	};
	using type = bond::derive_t<homotype>;

};
template <scalar_q U>
struct fadon<U>
:	fadon<U[2]>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
