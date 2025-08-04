#pragma once
#include "./any.cc"
#include "./wniplex.hh"// testing...

#include "../../process/all.hh"



XTAL_ENV_(push)
namespace xtal::atom::math::pade::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("wniplex")
{
	using _fit = bond::fit<>;
	using T_delta = typename _fit::delta_type;
	using T_sigma = typename _fit::sigma_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	auto constexpr pi = _fit::patio_1;
	/*/
	TRY_("operation")
	{
		T_aphex const x{-0.125, -0.125};
		//\
		wniplex_t<T_alpha> u(x);
		auto u = wniplex_f(x);
		auto const [v_re, v_im] = u.reflexed();
		auto const [w_re, w_im] = u.resolved();
		TRUE_(check_f<-30>(v_re, process::math::imagine_f<0>(two*cos(two*pi*x))));
		TRUE_(check_f<-30>(v_im, process::math::imagine_f<1>(two*sin(two*pi*x))));

		u *= u;

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
