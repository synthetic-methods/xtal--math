#pragma once
#include "./any.cc"





#include "./cosy.hh"
XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("cosy")
{
	using U_fit   = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;

	TRY_("task")
	{
		double constexpr N_1pi = one*_std::numbers::pi_v<U_alpha>;
		double constexpr N_2pi = two*_std::numbers::pi_v<U_alpha>;
		
		double const t9 = 0.99999;
		double const t8 = 0.88888;
		double const t7 = 0.77777;
		double const t6 = 0.66666;
		double const t5 = 0.55555;
		double const t4 = 0.44444;
		double const t3 = 0.33333;
		double const t2 = 0.22222;
		double const t1 = 0.11111;

		TRUE_(check_f<- 1>(cosh(N_1pi*t9), cosy_t<2>{}.template method<~0>(t9)));
		TRUE_(check_f<-27>(cosh(N_1pi*t9), cosy_t<2>{}.template method< 3>(t9)));
		TRUE_(check_f<-36>(cosh(N_1pi*t9), cosy_t<2>{}.template method< 2>(t9)));
		TRUE_(check_f<-45>(cosh(N_1pi*t9), cosy_t<2>{}.template method< 1>(t9)));
	//	TRUE_(check_f<-00>(cosh(N_1pi*t9), cosy_t<2>{}.template method< 0>(t9)));

		TRUE_(check_f<- 1>(cos(N_1pi*t1), cosy_t<1>{}.template method<~0>(t1)));
		TRUE_(check_f<- 1>(cos(N_1pi*t1), cosy_t<1>{}.template method< 3>(t1)));
		TRUE_(check_f<-13>(cos(N_1pi*t1), cosy_t<1>{}.template method< 2>(t1)));
		TRUE_(check_f<-27>(cos(N_1pi*t1), cosy_t<1>{}.template method< 1>(t1)));
		TRUE_(check_f<-38>(cos(N_1pi*t1), cosy_t<1>{}.template method< 0>(t1)));

		TRUE_(check_f<- 1>(square_f(cos(N_1pi*t1)), cosy_t<1, 2>{}.template method<~0>(t1)));
		TRUE_(check_f<- 1>(square_f(cos(N_1pi*t1)), cosy_t<1, 2>{}.template method< 3>(t1)));
		TRUE_(check_f<-15>(square_f(cos(N_1pi*t1)), cosy_t<1, 2>{}.template method< 2>(t1)));
		TRUE_(check_f<-28>(square_f(cos(N_1pi*t1)), cosy_t<1, 2>{}.template method< 1>(t1)));
		TRUE_(check_f<-38>(square_f(cos(N_1pi*t1)), cosy_t<1, 2>{}.template method< 0>(t1)));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
