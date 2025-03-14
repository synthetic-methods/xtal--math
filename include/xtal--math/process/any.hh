#pragma once
#include "../any.hh"

#include "../flow/any.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <auto  ..._s> struct  value_template;
template <class ..._s> struct  class_template;
template <class ..._s> struct  any   : process::any<_s...> {};
template <class ..._s> using   any_t = confined_t<any<_s...>>;


////////////////////////////////////////////////////////////////////////////////

template <class T>
struct any<T>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

	public:// ...

		template <extent_type N_mask=-1>
		struct   attach {template <class R> using subtype = R;};
		
		template <extent_type N_mask>
		requires requires {typename S_::template   attach<N_mask>;}
		struct   attach<N_mask> :   S_::template   attach<N_mask>{};


		template <extent_type N_mask=-1>
		struct dispatch {template <class R> using subtype = R;};
		
		template <extent_type N_mask>
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
template <template <class ...> class T_, class ..._s> requires different_q<T_<_s...>, class_template<_s...>>
struct any<T_<_s...>> : any<class_template<_s...>>
{
};
template <template <auto  ...> class T_, auto  ..._s> requires different_q<T_<_s...>, value_template<_s...>>
struct any<T_<_s...>> : any<value_template<_s...>>
{
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
