#pragma once
#include "./any.hh"
#include "../../provision/saturation.hh"





XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

template <class ...As>
struct scaffold
:	occur::context<scaffold<As...>>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////


namespace xtal::occur
{////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

template <vector_q A, class ..._s>
struct context<process::math::zavalishin::scaffold<A, _s...>>
{
private:
	static_assert(incomplete_q<_s...>);
	XTAL_DEF_(set) N_pole =           fixed<A>::extent();
	XTAL_TYP_(set) U_pole = typename  fixed<A>::value_type;
	XTAL_TYP_(set) V_pole = unstruct_t<A>;

public:
	using superkind = context<>;

	template <class S>
	class subtype : public bond::compose_s<S, superkind>
	{
		using S_ = bond::compose_s<S, superkind>;

	public:// CONSTRUCT
		using S_::S_;

	public:// EXTERNAL
		using    resample_type = occur::resample_t<>;
		using       stage_type = occur::   stage_t<>;

	public:// INTERNAL
		using       state_type = atom::couple_t<U_pole[N_pole]>;
		using       slope_type = atom::couple_t<U_pole[N_pole]>;
		using       shape_type = atom::couple_t<V_pole[N_pole]>;

	public:// ATTEND
		using  gain_parameter  = occur::inferred_t<_s..., union GAIN, V_pole>;
		using  damp_parameter  = occur::inferred_t<_s..., union DAMP, V_pole>;
		using  fade_parameter  = occur::inferred_t<_s..., union FADE, V_pole>;
		using  zoom_parameter  = occur::inferred_t<_s..., union ZOOM, V_pole>;

	public:// ATTACH
		using shape_parameter  = occur::inferred_t<_s..., union SHAPE, shape_type>;

	public:// DISPATCH
		using order_attribute  = occur::inferred_t<_s..., union ORDER, bond::seek_s<1 + N_pole>>;

	};
};

template <scalar_q A>
struct context<process::math::zavalishin::scaffold<A>>
:	context<process::math::zavalishin::scaffold<A[2]>>
{
};
template <>
struct context<process::math::zavalishin::scaffold< >>
:	context<process::math::zavalishin::scaffold<typename bond::fit<>::alpha_type>>
{
};

template <bond::compose_q A, class ..._s>
struct context<process::math::zavalishin::scaffold<A, _s...>>
:	bond::compose<A
	,	context<process::math::zavalishin::scaffold<_s...>>
	>
{
};
template <incomplete_q A, class ..._s>
struct context<process::math::zavalishin::scaffold<A, _s...>>
:	bond::compose<context<A>
	,	context<process::math::zavalishin::scaffold<_s...>>
	>
{
};
//template <template <class ...> class T_, class ..._s>
//struct context<T_<_s...>> : context<process::math::zavalishin::scaffold<_s...>>
//{
//};


///////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
