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
struct base
:	occur::auxiliary<base<As...>>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////


namespace xtal::occur
{////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <vector_q A, class ..._s>
struct auxiliary<process::math::zavalishin::base<A, _s...>>
{
private:
	static_assert(incomplete_q<_s...>);
	XTAL_DEF_(set) N_pole =           fixed<A>::extent();
	XTAL_TYP_(set) U_pole = typename  fixed<A>::value_type;
	XTAL_TYP_(set) V_pole = unstruct_t<A>;

public:
	using superkind = auxiliary<>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;

	public:// INTERNAL
		using   data_type = atom::math::dot_t<U_pole[N_pole]>;
		using codata_type = atom::math::dot_t<V_pole[N_pole]>;

	public:// DISPATCH
		using  order_attribute = occur::inferred_t<_s..., union ORDER, bond::seek_to_t<1 + N_pole>>;

	};
};

template <scalar_q A>
struct auxiliary<process::math::zavalishin::base<A>>
:	auxiliary<process::math::zavalishin::base<A[2]>>
{
};
template <>
struct auxiliary<process::math::zavalishin::base< >>
:	auxiliary<process::math::zavalishin::base<typename bond::fit<>::alpha_type>>
{
};

template <bond::compose_q A, class ..._s>
struct auxiliary<process::math::zavalishin::base<A, _s...>>
:	bond::compose<A
	,	auxiliary<process::math::zavalishin::base<_s...>>
	>
{
};
template <incomplete_q A, class ..._s>
struct auxiliary<process::math::zavalishin::base<A, _s...>>
:	bond::compose<auxiliary<A>
	,	auxiliary<process::math::zavalishin::base<_s...>>
	>
{
};
//template <template <class ...> class T_, class ..._s>
//struct auxiliary<T_<_s...>> : auxiliary<process::math::zavalishin::base<_s...>>
//{
//};


///////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
