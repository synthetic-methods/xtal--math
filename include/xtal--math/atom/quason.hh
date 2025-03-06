#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::atom::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Extends `group_multiplication` as a distinct quantity for control-bundling.
*/
template <class ...Us>	struct  quason_multiplication;
template <class ...Us>	using   quason_multiplication_t = typename quason_multiplication<Us...>::type;
template <class ...Us>	concept quason_multiplication_q = bond::tag_infixed_p<quason_multiplication_t, Us...>;

XTAL_DEF_(let) quason_multiplication_f = [] XTAL_1FN_(call) (_detail::fake_f<quason_multiplication_t>);


template <scalar_q ...Us> requires common_q<Us...>
struct quason_multiplication<Us ...>
:	quason_multiplication<common_t<Us...>[sizeof...(Us)]>
{
};
template <class ...Us>
struct quason_multiplication
{
private:
	template <class T>
	using endotype = typename group_multiplication<Us...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<quason_multiplication_t>>;

public:
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
/*!
\brief   Extends `group_addition` as a distinct quantity for control-bundling.
*/
template <class ...Us>	struct  quason_addition;
template <class ...Us>	using   quason_addition_t = typename quason_addition<Us...>::type;
template <class ...Us>	concept quason_addition_q = bond::tag_infixed_p<quason_addition_t, Us...>;

XTAL_DEF_(let) quason_addition_f = [] XTAL_1FN_(call) (_detail::fake_f<quason_addition_t>);


template <scalar_q ...Us> requires common_q<Us...>
struct quason_addition<Us ...>
:	quason_addition<common_t<Us...>[sizeof...(Us)]>
{
};
template <class ...Us>
struct quason_addition
{
private:
	template <class T>
	using endotype = typename group_addition<Us...>::template homotype<T>;

	template <class T>
	using holotype = bond::compose_s<endotype<T>, bond::tag<quason_addition_t>>;

public:
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

template <class U, auto N> struct quason<_std::plus       <U>   [N]> :       quason_addition<U   [N]> {};
template <class U, auto N> struct quason<_std::plus       <U>(&)[N]> :       quason_addition<U(&)[N]> {};

template <class U, auto N> struct quason<_std::multiplies <U>   [N]> : quason_multiplication<U   [N]> {};
template <class U, auto N> struct quason<_std::multiplies <U>(&)[N]> : quason_multiplication<U(&)[N]> {};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
