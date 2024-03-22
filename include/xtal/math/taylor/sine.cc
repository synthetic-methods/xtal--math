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
		int constexpr N_lim = 0b11;
		int constexpr N_bas = 0;// e
		double w = 0.123;

		xtal::echo(std::sin(w));
		xtal::echo(sine_t< 1, 0>::template function<N_lim>(w));
		xtal::echo();

		xtal::echo(std::sin(w)/w);
		xtal::echo(sine_t< 1, 1>::template function<N_lim>(w));
		xtal::echo();

		xtal::echo(std::sinh(w));
		xtal::echo(sine_t< 2, 0>::template function<N_lim>(w));
		xtal::echo();

		xtal::echo(std::sinh(w)/w);
		xtal::echo(sine_t< 2, 1>::template function<N_lim>(w));
		xtal::echo();
		
		xtal::echo();

		xtal::echo(std::asin(w));
		xtal::echo(sine_t<-1, 0>::template function<N_lim>(w));
		xtal::echo();

		xtal::echo(std::asin(w)/w);
		xtal::echo(sine_t<-1, 1>::template function<N_lim>(w));
		xtal::echo();

		xtal::echo(std::asinh(w));
		xtal::echo(sine_t<-2, 0>::template function<N_lim>(w));
		xtal::echo();

		xtal::echo(std::asinh(w)/w);
		xtal::echo(sine_t<-2, 1>::template function<N_lim>(w));
		xtal::echo();


		TRUE_(true);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
