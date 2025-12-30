#pragma once
#include "./any.cc"





#include "./dot.hh"// testing...
XTAL_ENV_(push)
namespace xtal::atom::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("dot")
{
	using U_fit = bond::fit<>;
	using U_delta = typename U_fit::delta_type;
	using U_sigma = typename U_fit::sigma_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;
	using W_alpha = atom::field_t<_xtd::plus_multiplies<U_alpha>>;

	TRY_("dot (different scalar type)")
	{
		couple_t<U_alpha[2]> u0{2, 3};
		dot_t   <U_alpha[2]> u1{4, 9};
		TRUE_(u0*u1 == 35.);
		TRUE_(u1*u0 == 35.);

	}
	TRY_("dot (same scalar type)")
	{
		dot_t<U_alpha[2]> u0{2, 3};
		dot_t<U_alpha[2]> u1{4, 9};
		TRUE_(u0*u1 == 35.);
		TRUE_(u1*u0 == 35.);

	}
	TRY_("dot (same vector type)")
	{
		dot_t<W_alpha[2]> u0{{1, 2}, {2, 3}};
		dot_t<W_alpha[2]> u1{{1, 4}, {4, 9}};
		TRUE_(u0*u1 == W_alpha{9., 35.});
		TRUE_(u1*u0 == W_alpha{9., 35.});

	}
	TRY_("dot (same vector type via matrix)")
	{
		//\
		couple_t<dot_t<W_alpha[2]>[2]> u_{{{1, 2}, {2, 3}}, {{1, 4}, {4, 9}}};
		dot_t<W_alpha[2][2]> u_{{{1, 2}, {2, 3}}, {{1, 4}, {4, 9}}};
		auto const &u0 = u_[0];
		auto const &u1 = u_[1];
		TRUE_(u0*u1 == W_alpha{9., 35.});
		TRUE_(u1*u0 == W_alpha{9., 35.});
		TRUE_(u_.product() == W_alpha{9., 35.});

	}
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
