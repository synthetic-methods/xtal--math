#pragma once
#include "./any.cc"
#include "./sine.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::taylor::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("sine")
{
	using _op = bond::operate<size_type>;
	using T_sigma = typename _op::sigma_type;
	using T_delta = typename _op::delta_type;
	using T_alpha = typename _op::alpha_type;
	static constexpr T_alpha one =  1;
	static constexpr T_alpha two =  2;
	static constexpr T_alpha ten = 10;

	TRY_("task")
	{
		double constexpr N_2pi = 6.283185307179586476925286766559005768L;
		int constexpr N_lim = 0b111;
		int constexpr N_bas = 0;// e
		
		double const t = 0.11111, t_2pi = N_2pi*t;

		TRUE_(check_f<15>(std:: sin (t)/(1), sine_t<+1, -0>::template function<N_lim>(t)));
		TRUE_(check_f<15>(std:: sin (t)/(t), sine_t<+1, -1>::template function<N_lim>(t)));
		TRUE_(check_f<15>(std:: sinh(t)/(1), sine_t<+2, -0>::template function<N_lim>(t)));
		TRUE_(check_f<15>(std:: sinh(t)/(t), sine_t<+2, -1>::template function<N_lim>(t)));
		TRUE_(check_f<15>(std::asin (t)/(1), sine_t<-1, -0>::template function<N_lim>(t)));
		TRUE_(check_f<15>(std::asin (t)/(t), sine_t<-1, -1>::template function<N_lim>(t)));
		TRUE_(check_f<15>(std::asinh(t)/(1), sine_t<-2, -0>::template function<N_lim>(t)));
		TRUE_(check_f<15>(std::asinh(t)/(t), sine_t<-2, -1>::template function<N_lim>(t)));

		TRUE_(check_f<15>(std:: sin (t_2pi)/(N_2pi), sine_t<+1, +1>::template function<N_lim>(t)));
		TRUE_(check_f<15>(std:: sinh(t_2pi)/(N_2pi), sine_t<+2, +1>::template function<N_lim>(t)));
		TRUE_(check_f<15>(std::asin (t_2pi)/(N_2pi), sine_t<-1, +1>::template function<N_lim>(t)));
		TRUE_(check_f<15>(std::asinh(t_2pi)/(N_2pi), sine_t<-2, +1>::template function<N_lim>(t)));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
