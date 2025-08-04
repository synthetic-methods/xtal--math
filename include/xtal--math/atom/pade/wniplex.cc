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
	using U_fit = bond::fit<>;
	using U_delta = typename U_fit::delta_type;
	using U_sigma = typename U_fit::sigma_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	using W_alpha = atom::couple_t<U_alpha[2]>;
	using W_aphex = atom::couple_t<U_aphex[2]>;

	auto constexpr pi = U_fit::patio_1;
	/**/
	TRY_("operation")
	{
		U_aphex const x{-0.125, -0.125};
		//\
		wniplex_t<U_alpha> u(x);
		auto u = wniplex_f(x);
		auto const [v_re, v_im] = u.reflection();
		auto const [w_re, w_im] = u.resolution();
		TRUE_(check_f<-30>(v_re, process::math::imagine_f<0>(two*cos(two*pi*x))));
		TRUE_(check_f<-30>(v_im, process::math::imagine_f<1>(two*sin(two*pi*x))));

		u *= u;
		TRUE_(check_f<-2>(u, u.flipped(U_alpha{1})));

		u *= u.flipped(U_alpha{-1});

		auto const [u_re, u_im] = destruct_f(u.signum());
		auto const [u_dn, u_up] =           (u.magnum());

		TRUE_(check_f<-2>(u_im, U_alpha{0}));
		TRUE_(check_f<-2>(u_re, U_alpha{1}));
		TRUE_(check_f<-2>(u_up, U_alpha{1}));
		TRUE_(check_f<-2>(u_dn, U_alpha{1}));

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
