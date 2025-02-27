#pragma once
#include "../any.hxx"






XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <class ..._s> struct  abstract;
template <class ..._s> struct  any   : process::any<_s...> {};
template <class ..._s> using   any_t = confined_t<any<_s...>>;


////////////////////////////////////////////////////////////////////////////////

template <vector_q A, class ..._s>
struct any<abstract<A, _s...>>
{
	template <class S>
	class subtype : public S
	{
		static_assert(incomplete_q<_s...>);
		XTAL_DEF_(set) N_pole = vector_n<A>;
		using          U_pole = vector_u<A>;

	public:// CONSTRUCT
		using S::S;

	public:// EXTERNAL
		using    resample_type = occur::resample_t<>;
		using       stage_type = occur::   stage_t<>;

	public:// INTERNAL
		using        pole_size = constant_t<N_pole>;
		using        pole_type =            U_pole ;
		using coefficient_type = unstruct_u<U_pole>;

		using       state_type = atom::couple_t<       pole_type[N_pole]>;
		using       slope_type = atom::couple_t<       pole_type[N_pole]>;
		using       shape_type = atom::couple_t<coefficient_type[N_pole]>;

	public:// ATTACH
		using     reshape_type = occur::inferred_t<_s..., union FILTER, union RESHAPE,       shape_type>;
		using      rezoom_type = occur::inferred_t<_s..., union FILTER, union  REZOOM, coefficient_type>;

	public:// ATTEND
		using      regain_type = occur::inferred_t<_s..., union FILTER, union  REGAIN, coefficient_type>;
		using      redamp_type = occur::inferred_t<_s..., union FILTER, union  REDAMP, coefficient_type>;
		using      refade_type = occur::inferred_t<_s..., union FILTER, union  REFADE, coefficient_type>;

	public:// DISPATCH
		using       order_type = occur::inferred_t<_s..., union FILTER, union   ORDER, bond::seek_s<1 + N_pole>>;
		using       patch_type = occur::inferred_t<_s..., union FILTER, union   PATCH, bond::seek_s<2>>;
		using       limit_type = occur::inferred_t<_s..., union FILTER, union   LIMIT, bond::seek_s<2>>;

	};
};


////////////////////////////////////////////////////////////////////////////////

template <bond::compose_q A, class ..._s>
struct any<abstract<A, _s...>>
:	bond::compose<A, any<abstract<_s...>>>
{
};
template <incomplete_q A, class ..._s>
struct any<abstract<A, _s...>>
:	bond::compose<any<A>, any<abstract<_s...>>>
{
};
template <template <class ...> class T_, class ..._s>
struct any<T_<_s...>> : any<abstract<_s...>>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
