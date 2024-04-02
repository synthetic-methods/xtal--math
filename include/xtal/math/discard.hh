#pragma once
#include "./any.hh"
#include "./square.hh"





XTAL_ENV_(push)
namespace xtal::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_car>
struct discard;


////////////////////////////////////////////////////////////////////////////////

template <>
struct discard<0>
{
	template <class S>
	using subtype = bond::compose_s<S, bond::tag<process::chain>>;

};
template <>
struct discard<1>
{
	using subkind = bond::tag<process::chain>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_FN2 function(auto &&o)
		XTAL_0EX
		{
			auto const v = S_::template function<Is...>(o); using V = XTAL_TYP_(v);
			if constexpr (complex_field_q<V>) {
				return V {v.real(), v.imag()*XTAL_REF_(o)};
			}
			else {
				return v*XTAL_REF_(o);
			}
		};

	};
};
template <>
struct discard<2>
{
	using subkind = bond::tag<process::chain>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_FN2 function(auto &&o)
		XTAL_0EX
		{
			using U = XTAL_TYP_(o); using re = bond::realize<U>;
			return S_::template function<Is...>(square_f<1>(XTAL_REF_(o)));
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
