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
:	occur::codex<meta<As...>>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////


namespace xtal::occur
{////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <vector_q A, class ..._s>
struct codex<process::math::zavalishin::meta<A, _s...>>
{
private:
	static_assert(incomplete_q<_s...>);
	XTAL_DEF_(set) N_pole =           fixed<A>::extent();
	XTAL_TYP_(set) U_pole = typename  fixed<A>::value_type;
	XTAL_TYP_(set) V_pole = unstruct_t<A>;

public:
	using superkind = codex<>;

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
struct codex<process::math::zavalishin::meta<A>>
:	codex<process::math::zavalishin::meta<A[2]>>
{
};
template <>
struct codex<process::math::zavalishin::meta< >>
:	codex<process::math::zavalishin::meta<typename bond::fit<>::alpha_type>>
{
};

template <bond::compose_q A, class ..._s>
struct codex<process::math::zavalishin::meta<A, _s...>>
:	bond::compose<A
	,	codex<process::math::zavalishin::meta<_s...>>
	>
{
};
template <incomplete_q A, class ..._s>
struct codex<process::math::zavalishin::meta<A, _s...>>
:	bond::compose<codex<A>
	,	codex<process::math::zavalishin::meta<_s...>>
	>
{
};
//template <template <class ...> class T_, class ..._s>
//struct codex<T_<_s...>> : codex<process::math::zavalishin::meta<_s...>>
//{
//};


///////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
