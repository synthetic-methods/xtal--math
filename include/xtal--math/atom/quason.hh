#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::atom::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
///\
Extends `multiplicative_group` as a distinct quantity for control-bundling. \

template <class ...Us>	struct  multiplicative_quason;
template <class ...Us>	using   multiplicative_quason_t = typename multiplicative_quason<Us...>::type;
template <class ...Us>	concept multiplicative_quason_q = bond::fixed_tagged_with_p<multiplicative_quason_t, Us...>;

XTAL_DEF_(let) multiplicative_quason_f = [] XTAL_1FN_(call) (_detail::fake_f<multiplicative_quason_t>);


template <scalar_q ...Us> requires common_q<Us...>
struct multiplicative_quason<Us ...>
:	multiplicative_quason<common_t<Us...>[sizeof...(Us)]>
{
};
template <class ...Us>
struct multiplicative_quason
{
	template <class T>
	using endotype = typename multiplicative_group<Us...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<multiplicative_quason_t>>;

	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;

	public:// CONSTRUCT
		using S_::S_;

	};
	using type = bond::derive_t<homotype>;

};


////////////////////////////////////////////////////////////////////////////////
///\
Extends `additive_group` as a distinct quantity for control-bundling. \

template <class ...Us>	struct  additive_quason;
template <class ...Us>	using   additive_quason_t = typename additive_quason<Us...>::type;
template <class ...Us>	concept additive_quason_q = bond::fixed_tagged_with_p<additive_quason_t, Us...>;

XTAL_DEF_(let) additive_quason_f = [] XTAL_1FN_(call) (_detail::fake_f<additive_quason_t>);


template <scalar_q ...Us> requires common_q<Us...>
struct additive_quason<Us ...>
:	additive_quason<common_t<Us...>[sizeof...(Us)]>
{
};
template <class ...Us>
struct additive_quason
{
	template <class T>
	using endotype = typename additive_group<Us...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<additive_quason_t>>;

	template <class T>
	class homotype : public holotype<T>
	{
		using S_ = holotype<T>;

	public:// CONSTRUCT
		using S_::S_;

	};
	using type = bond::derive_t<homotype>;

};


////////////////////////////////////////////////////////////////////////////////

template <class T        > struct quason;

template <class U, auto N> struct quason<_std::plus       <U>   [N]> :       additive_quason<U   [N]> {};
template <class U, auto N> struct quason<_std::plus       <U>(&)[N]> :       additive_quason<U(&)[N]> {};

template <class U, auto N> struct quason<_std::multiplies <U>   [N]> : multiplicative_quason<U   [N]> {};
template <class U, auto N> struct quason<_std::multiplies <U>(&)[N]> : multiplicative_quason<U(&)[N]> {};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
