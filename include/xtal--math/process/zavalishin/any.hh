#pragma once
#include "../any.hxx"






XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <auto  ..._s> struct  value_template;
template <class ..._s> struct  class_template;
template <class ..._s> struct  any   : process::any<_s...> {};
template <class ..._s> using   any_t = confined_t<any<_s...>>;


////////////////////////////////////////////////////////////////////////////////

template <vector_q A, class ..._s>
struct any<class_template<A, _s...>>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

		static_assert(incomplete_q<_s...>);
		XTAL_DEF_(set) N_pole = vector_n<A>;
		using          U_pole = vector_u<A>;

	public:// CONSTRUCT
		using S_::S_;

	public:// EXTERNAL
		using    resample_type = occur::resample_t<>;
		using       stage_type = occur::   stage_t<>;

	public:// INTERNAL
		using        pole_size = constant_t<N_pole>;
		using        pole_type =            U_pole ;
		using       scale_type = unstruct_u<U_pole>;

		using       state_type = atom::couple_t< pole_type[N_pole]>;
		using       slope_type = atom::couple_t< pole_type[N_pole]>;
		using       shape_type = atom::couple_t<scale_type[N_pole]>;

	public:// ATTEND
		using      regain_type = occur::inferred_t<_s..., union  REGAIN, scale_type>;
		using      redamp_type = occur::inferred_t<_s..., union  REDAMP, scale_type>;
		using      refade_type = occur::inferred_t<_s..., union  REFADE, scale_type>;
		using      rezoom_type = occur::inferred_t<_s..., union  REZOOM, scale_type>;

	public:// ATTACH
		using     reshape_type = occur::inferred_t<_s..., union RESHAPE, shape_type>;

	public:// DISPATCH
		using       order_type = occur::inferred_t<_s..., union   ORDER, bond::seek_s<1 + N_pole>>;
		using       patch_type = occur::inferred_t<_s..., union   PATCH, bond::seek_s<2>>;
		using       limit_type = occur::inferred_t<_s..., union   LIMIT, bond::seek_s<2>>;

	public:// ...

		template <size_type N_mask=1>
		struct   attach {template <class R> using subtype = R;};
		
		template <size_type N_mask>
		requires requires {typename S_::template   attach<N_mask>;}
		struct   attach<N_mask> :   S_::template   attach<N_mask>{};


		template <size_type N_mask=1>
		struct dispatch {template <class R> using subtype = R;};
		
		template <size_type N_mask>
		requires requires {typename S_::template dispatch<N_mask>;}
		struct dispatch<N_mask> :   S_::template dispatch<N_mask>{};

	};
};


////////////////////////////////////////////////////////////////////////////////

template <bond::compose_q A, class ..._s>
struct any<class_template<A, _s...>>
:	bond::compose<A, any<class_template<_s...>>>
{
};
template <incomplete_q A, class ..._s>
struct any<class_template<A, _s...>>
:	bond::compose<any<A>, any<class_template<_s...>>>
{
};
template <template <class ...> class T_, class ..._s>
struct any<T_<_s...>> : any<class_template<_s...>>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
