#pragma once
#include "./any.hh"
#include "../../atom/dot.hh"
#include "../../provision/zavalishin/shaped.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines associated types.
*/
template <class ...As>
struct meta
:	occur::meta<meta<As...>>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////


namespace xtal::occur
{////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <vector_q A, class ..._s>
struct meta<process::math::zavalishin::meta<A, _s...>>
{
private:
	static_assert(incomplete_q<_s...>);
	XTAL_DEF_(set) N_pole =           fixed<A>::extent();
	XTAL_TYP_(set) U_pole = typename  fixed<A>::value_type;
	XTAL_TYP_(set) V_pole = unstruct_t<A>;

public:
	using superkind = meta<>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;

	public:// INTERNAL
		using   data_type      = atom::math::dot_t<U_pole[N_pole]>;
		using codata_type      = atom::math::dot_t<V_pole[N_pole]>;

	public:// DISPATCH
		using  order_attribute = occur::inferred_t<_s..., union ORDER, bond::seek_s<1 + N_pole>>;

	};
};

template <scalar_q A>
struct meta<process::math::zavalishin::meta<A>>
:	meta<process::math::zavalishin::meta<A[2]>>
{
};
template <>
struct meta<process::math::zavalishin::meta< >>
:	meta<process::math::zavalishin::meta<typename bond::fit<>::alpha_type>>
{
};

template <bond::compose_q A, class ..._s>
struct meta<process::math::zavalishin::meta<A, _s...>>
:	bond::compose<A
	,	meta<process::math::zavalishin::meta<_s...>>
	>
{
};
template <incomplete_q A, class ..._s>
struct meta<process::math::zavalishin::meta<A, _s...>>
:	bond::compose<meta<A>
	,	meta<process::math::zavalishin::meta<_s...>>
	>
{
};
//template <template <class ...> class T_, class ..._s>
//struct meta<T_<_s...>> : meta<process::math::zavalishin::meta<_s...>>
//{
//};


///////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
