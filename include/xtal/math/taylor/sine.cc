#pragma once
#include "./any.cc"
#include "./sine.ii"// testing...






XTAL_ENV_(push)
namespace xtal::math::taylor::__test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("sine")
{
	TRY_("task")
	{
		double constexpr N_2pi = 6.283185307179586476925286766559005768L;
		int constexpr N_lim = 0b1111;
		int constexpr N_bas = 0;// e
		
		double const t = 0.11111, t_2pi = N_2pi*t;

		TRUE_(check_f<22>(std:: sin (t_2pi)/(N_2pi), sine_t<+1, 0, 0, (-0)>::template function<N_lim>(t)));
		TRUE_(check_f<22>(std:: sin (t_2pi)/(t_2pi), sine_t<+1, 0, 0, (-1)>::template function<N_lim>(t)));
		TRUE_(check_f<22>(std:: sinh(t_2pi)/(N_2pi), sine_t<+1, 1, 0, (-0)>::template function<N_lim>(t)));
		TRUE_(check_f<22>(std:: sinh(t_2pi)/(t_2pi), sine_t<+1, 1, 0, (-1)>::template function<N_lim>(t)));
		TRUE_(check_f<22>(std::asin (t_2pi)/(N_2pi), sine_t<-1, 0, 0, (-0)>::template function<N_lim>(t)));
		TRUE_(check_f<22>(std::asin (t_2pi)/(t_2pi), sine_t<-1, 0, 0, (-1)>::template function<N_lim>(t)));
		TRUE_(check_f<22>(std::asinh(t_2pi)/(N_2pi), sine_t<-1, 1, 0, (-0)>::template function<N_lim>(t)));
		TRUE_(check_f<22>(std::asinh(t_2pi)/(t_2pi), sine_t<-1, 1, 0, (-1)>::template function<N_lim>(t)));

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
