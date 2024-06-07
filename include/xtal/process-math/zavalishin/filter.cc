#pragma once
#include "./any.cc"
#include "./filter.hh"// testing...

#include "./prewarping.hh"



XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("filter")
{
	TRY_("instantiation")
	{
		using U = double;
		//\
		using Y = filter_t<U[2]>;
		using Y = filter_t<prewarping<>, U[2]>;
		//\
		using Z = processor::monomer_t<prewarping<>, Y>;
		using Z = processor::monomer_t<Y>;

		Y y{}; y <<= occur::sample_t<>{44100};

		echo(y(0.5, 1.0));

		_std::array<U, 4> u_{1, 2, 0, 1};
		_std::array<U, 4> f_{3, 3, 3, 3};

		//\
		auto z = Z::bind_f(let_f(u_), let_f(f_));
		auto z = Z::bind_f(u_, f_);
//		TRUE_(true);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
