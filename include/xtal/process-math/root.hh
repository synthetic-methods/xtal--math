#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_pow=1> struct root;
template <int N_pow=1> using  root_t = process::confined_t<root<N_pow>>;
template <int N_pow=1>
XTAL_FN2 root_f(auto &&o)
XTAL_0EX
{
	return root_t<N_pow>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int N_pow>
struct root//<N_pow>
{
	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <int N_lim=-1>// requires sign_p<N_pow, 0>
		XTAL_DEF_(return,inline)
		XTAL_FN1 function(auto const &o)
		XTAL_0EX
		{
			using Op = bond::operate<decltype(o)>;
			/*/
			return Op::template root_f<N_pow, N_lim>(o);
			/*/
			XTAL_IF0
			XTAL_0IF (N_pow ==  2) {return             sqrt(o);}
			XTAL_0IF (N_pow ==  1) {return                 (o);}
			XTAL_0IF (N_pow == -1) {return Op::alpha_1/    (o);}
			XTAL_0IF (N_pow == -2) {return Op::alpha_1/sqrt(o);}
			/***/
		}

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
