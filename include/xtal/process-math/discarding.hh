#pragma once
#include "./any.hh"
#include "./square.hh"





XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_car>
struct discarding;


////////////////////////////////////////////////////////////////////////////////

template <>
struct discarding<0>
{
	template <class S>
	using subtype = bond::compose_s<S, bond::tag<process::link>>;

};
template <>
struct discarding<1>
{
	using subkind = bond::tag<process::link>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_FN1 function(auto &&o)
		XTAL_0EX
		{
			auto z = S_::template function<Is...>(o); using Z = XTAL_TYP_(z);
			XTAL_IF0
			XTAL_0IF (complex_number_q<Z>) {
				//\
				devolved_f(z)[1] *= XTAL_REF_(o);
				get<1>(z) *= XTAL_REF_(o);
				return z;
			}
			XTAL_0IF (complex_field_q<Z>) {
				return complexion_f(z.real(), z.imag()*XTAL_REF_(o));
			}
			XTAL_0IF_(default) {
				return z*XTAL_REF_(o);
			}
		};

	};
};
template <>
struct discarding<2>
{
	using subkind = bond::tag<process::link>;

	template <class S>
	class subtype: public bond::compose_s<S, subkind>
	{
		using S_ = bond::compose_s<S, subkind>;

	public:
		using S_::S_;

		template <auto ...Is>
		XTAL_DEF_(return,inline)
		XTAL_FN1 function(auto &&o)
		XTAL_0EX
		{
			using U = XTAL_TYP_(o); using _op = bond::operate<U>;
			return S_::template function<Is...>(square_f<1>(XTAL_REF_(o)));
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
