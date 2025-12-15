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

	//\
	using W_alpha = atom::couple_t<U_alpha[2]>;
	using W_alpha = field_t<_xtd::plus_multiplies<U_alpha>[2]>;
	using W_aphex = _std::complex<W_alpha>;

	auto constexpr pi = U_fit::patio_1;
	auto constexpr qi = U_fit::patio_2;
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
	/**/
	TRY_("`wniplex_t<couple_t<...>>` operation")
	{
		U_aphex const x0{ 0.11, 0.33};
		U_aphex const x1{ 0.22, 0.44};
		W_aphex const xA{{ 0.11, 0.22}, { 0.33, 0.44}};
		//\
		wniplex_t<U_alpha> u(xA);
		auto u0 = wniplex_f(x0);
		auto u1 = wniplex_f(x1);
		auto uA = wniplex_f(xA);
		auto const [v0_up, v0_dn] = u0.reflection();
		auto const [v1_up, v1_dn] = u1.reflection();
		auto const [vA_up, vA_dn] = uA.reflection();
		auto const [w0_up, w0_dn] = u0.resolution();
		auto const [w1_up, w1_dn] = u1.resolution();
		auto const [wA_up, wA_dn] = uA.resolution();

		TRUE_(check_f<-1>(v0_up, complexion_f(get<0>(vA_up.real()), get<0>(vA_up.imag()))));
		TRUE_(check_f<-1>(v0_dn, complexion_f(get<0>(vA_dn.real()), get<0>(vA_dn.imag()))));

		TRUE_(check_f<-1>(v1_up, complexion_f(get<1>(vA_up.real()), get<1>(vA_up.imag()))));
		TRUE_(check_f<-1>(v1_dn, complexion_f(get<1>(vA_dn.real()), get<1>(vA_dn.imag()))));

	}
	/***/
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
