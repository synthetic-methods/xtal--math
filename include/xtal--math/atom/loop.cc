#pragma once
#include "./any.cc"

#include "../process/all.hh"



#include "./loop.hh"
XTAL_ENV_(push)
namespace xtal::atom::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("loop")
{
	using U_fit = bond::fit<>;
	using U_delta = typename U_fit::delta_type;
	using U_sigma = typename U_fit::sigma_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	TRY_("creation")
	{
		using W = loop_t<U_alpha[4]>;
		W w{};

		w.push(0.); TRUE_(0 == w[0] and 0 == w[1] and 0 == w[2] and 0 == w[3]);
		w.push(1.); TRUE_(1 == w[0] and 0 == w[1] and 0 == w[2] and 0 == w[3]);
		w.push(2.); TRUE_(2 == w[0] and 1 == w[1] and 0 == w[2] and 0 == w[3]);
		w.push(3.); TRUE_(3 == w[0] and 2 == w[1] and 1 == w[2] and 0 == w[3]);
		w.push(4.); TRUE_(4 == w[0] and 3 == w[1] and 2 == w[2] and 1 == w[3]);

		TRUE_(check_f<-1>(w[0], w.template element<0>()));
		TRUE_(check_f<-1>(w[1], w.template element<1>()));
		TRUE_(check_f<-1>(w[2], w.template element<2>()));
		TRUE_(check_f<-1>(w[3], w.template element<3>()));
		TRUE_(check_f<-1>(w[4], w.template element<4>()));

		TRUE_(check_f<-1>(w[0], get<0>(w)));
		TRUE_(check_f<-1>(w[1], get<1>(w)));
		TRUE_(check_f<-1>(w[2], get<2>(w)));
		TRUE_(check_f<-1>(w[3], get<3>(w)));
		TRUE_(check_f<-1>(w[4], get<4>(w)));

		TRUE_(check_f<-1>(w(0.0), 4.00));
		TRUE_(check_f<-1>(w(0.5), 3.75));
		TRUE_(check_f<-1>(w(1.0), 3.00));
		TRUE_(check_f<-1>(w(1.5), 2.50));
		TRUE_(check_f<-1>(w(2.0), 2.00));
		TRUE_(check_f<-1>(w(2.5), 1.25));
		TRUE_(check_f<-1>(w(3.0), 1.00));
		TRUE_(check_f<-1>(w(3.5), 2.50));

	}
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
