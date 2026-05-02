#pragma once
#include "./any.cc"

#include "./phason.hh"



#include "./quason.hh"
XTAL_ENV_(push)
namespace xtal::atom::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/*/
TAG_("atom", "quason")
{
	using U_fit = bond::fit<>;
	using U_delta = typename U_fit::delta_type;
	using U_sigma = typename U_fit::sigma_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	using U_phi = phason_t<U_alpha[2]>;

	using B1  = block_t<U_alpha[1]>;
	using B2  = block_t<U_alpha[2]>;
	using B3  = block_t<U_alpha[3]>;
	using B4  = block_t<U_alpha[4]>;

	using Q1  = quason_t<U_alpha[1]>;
	using Q2  = quason_t<U_alpha[2]>;
	using Q3  = quason_t<U_alpha[3]>;
	using Q4  = quason_t<U_alpha[4]>;
	using Q4_ = quason_t<U_phi, U_alpha, U_alpha, U_alpha>;
	
	TRY_("partial construction")
	{
		Q4  q4 {1000};
		TRUE_(q4 == Q4{1000, 0, 0, 0});

		Q4_ q4_{{0, U_fit::dnsilon_f(45)}, 1, 2, 3};
		//\
		B4  b4 = q4_;
		B4  b4(q4_);
	//	echo_(b4[0], b4[1], b4[2], b4[3]);//NOTE: Unordered!


	}
	TRY_("quason addition")
	{
		Q2 d2_0{2, 2};
		Q2 d2_1{5, 7};

		TRUE_(d2_0+d2_1 == Q2{ 7,  9});
		d2_0 += d2_1;
		TRUE_(d2_0      == Q2{ 7,  9});

	}
	TRY_("quason multiplication")
	{
		Q2 d2_0{2, 2};
		Q2 d2_1{5, 7};

		TRUE_(d2_0*d2_1 == Q2{10, 14});
		d2_0 *= d2_1;
		TRUE_(d2_0      == Q2{10, 14});

		using W =  quason_t<U_aphex, U_alpha>;
		auto  x =  W{2, 3};
		auto  y =  W{4, 9};
		auto  z =  y*U_alpha{3};
		TRUE_(z == W{12, 27});

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
