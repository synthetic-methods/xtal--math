#pragma once
#include "./any.hh"
#include "../../atom/dot.hh"
#include "../../scheme/zavalishin/shaped.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Defines associated types.
*/
template <class ...As>
struct base
:	process::occurrence<base<As...>>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////


namespace xtal::process
{////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <vector_array_q A, class ..._s>
struct occurrence<process::math::zavalishin::base<A, _s...>>
{
private:
	static_assert(incomplete_q<_s...>);
	XTAL_DEF_(set) N_pole =           fixed<A>::extent();
	XTAL_TYP_(set) U_pole = typename  fixed<A>::value_type;
	XTAL_TYP_(set) V_pole = unstruct_t<A>;

public:
	using superkind = occurrence<>;

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

template <scalar_array_q A>
struct occurrence<process::math::zavalishin::base<A>>
:	occurrence<process::math::zavalishin::base<A[2]>>
{
};
template <>
struct occurrence<process::math::zavalishin::base< >>
:	occurrence<process::math::zavalishin::base<typename bond::fit<>::alpha_type>>
{
};

template <bond::compose_q A, class ..._s>
struct occurrence<process::math::zavalishin::base<A, _s...>>
:	bond::compose<A
	,	occurrence<process::math::zavalishin::base<_s...>>
	>
{
};
template <incomplete_q A, class ..._s>
struct occurrence<process::math::zavalishin::base<A, _s...>>
:	bond::compose<occurrence<A>
	,	occurrence<process::math::zavalishin::base<_s...>>
	>
{
};
//template <template <class ...> class T_, class ..._s>
//struct occurrence<T_<_s...>> : occurrence<process::math::zavalishin::base<_s...>>
//{
//};


///////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
