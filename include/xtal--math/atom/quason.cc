#pragma once
#include "./any.cc"

#include "./phason.hh"



#include "./quason.hh"// testing...
XTAL_ENV_(push)
namespace xtal::atom::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/*/
TAG_("atom", "quason")
{
	using _fit = bond::fit<>;
	using T_delta = typename _fit::delta_type;
	using T_sigma = typename _fit::sigma_type;
	using T_alpha = typename _fit::alpha_type;
	using T_aphex = typename _fit::aphex_type;

	using T_phi = phason_t<T_alpha[2]>;

	using B1  = block_t<T_alpha[1]>;
	using B2  = block_t<T_alpha[2]>;
	using B3  = block_t<T_alpha[3]>;
	using B4  = block_t<T_alpha[4]>;

	using Q1  = quason_t<T_alpha[1]>;
	using Q2  = quason_t<T_alpha[2]>;
	using Q3  = quason_t<T_alpha[3]>;
	using Q4  = quason_t<T_alpha[4]>;
	using Q4_ = quason_t<T_phi, T_alpha, T_alpha, T_alpha>;
	
	TRY_("partial construction")
	{
		Q4  q4 {1000};
		TRUE_(q4 == Q4{1000, 0, 0, 0});

		Q4_ q4_{{0, _fit::dnsilon_f(45)}, 1, 2, 3};
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

		using W =  quason_t<T_aphex, T_alpha>;
		auto  x =  W{2, 3};
		auto  y =  W{4, 9};
		auto  z =  y*T_alpha{3};
		TRUE_(z == W{12, 27});

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
