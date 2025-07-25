#pragma once
#include "../any.hxx"






XTAL_ENV_(push)
namespace xtal::process::math::pade
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

#include "./_detail.ii"


////////////////////////////////////////////////////////////////////////////////

template <auto  ..._s> struct  value_template;
template <class ..._s> struct  class_template;

template <class ..._s> struct  any   : process::any<_s...> {};
template <class ..._s> using   any_t = confined_t<any<_s...>>;


////////////////////////////////////////////////////////////////////////////////

template <complete_q T>
struct any<T>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;
		using T_ = typename S_::self_type;

	public:// CONSTRUCT
		using S_::S_;

	public:// DISPATCH
		using limit_type = occur::inferred_t<union LIMIT, bond::seek_s<(1<<3)>>;

		template <size_type N_mask=1>
		struct dispatch
		{
			template <class R>
			using subtype = bond::compose_s<R, provision::voiced<void
			,	typename T_::limit_type::template dispatch<N_mask>
			>>;

		};

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
