#pragma once
#include "./any.cc"

#include "../atom/dot.hh"



#include "./dent.hh"
XTAL_ENV_(push)
namespace xtal::flow::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("dent")
{
	using U_alpha = typename bond::fit<>::alpha_type;
	//\
	using W_alpha = atom::couple_t<U_alpha[2]>;
	using W_alpha = atom::math::dot_t<U_alpha[2]>;

	TRY_("task")
	{
		using X_vector = atom::bundle_t<U_alpha, W_alpha>;
		using X_matrix = atom::bundle_t<X_vector[2]>;

		using D22   = X_matrix;
		using D22_0 = dent_s<X_matrix, 0>;
		using D22_1 = dent_s<X_matrix, 1>;

	//	TRUE_(same_q<ordinal_constant_t<0>, typename D22_0           ::head_type>);
	//	TRUE_(same_q<ordinal_constant_t<1>, typename D22_1           ::head_type>);
	//	TRUE_(same_q<X_vector             , typename D22_0::tail_type::head_type>);
	//	TRUE_(same_q<X_vector             , typename D22_1::tail_type::head_type>);

	//	using E22_1 = typename D22_1::data_type;

	//	static_assert( xtd::fungible_with<typename D22_1::data_type, X_matrix>);
	//	static_assert(flow::math::dent_q<         D22_1           , X_matrix>);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
