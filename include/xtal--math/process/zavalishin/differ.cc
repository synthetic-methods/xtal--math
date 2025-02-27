#pragma once
#include "./any.cc"





#include "./differ.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::zavalishin::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

namespace _differ
{///////////////////////////////////////////////////////////////////////////////

struct reparameterized
{
	template <class S>
	class subtype : public S
	{
	public:
		using S::S;

		using typename S::pole_size;
		using order_type = occur::inferred_t<union ORDER, bond::seek_s<1 + pole_size{}()>>;

	};
};


}///////////////////////////////////////////////////////////////////////////////

/**/
TAG_("differ")
{
	using A_pole = typename bond::fit<>::alpha_type[1];
	auto constexpr N_pole = fixed_n<A_pole>;

	using R_order = occur::inferred_t<union ORDER, bond::seek_s<1 + N_pole>>;

	//\
	using A_differ = differ<A_pole>;
	using A_differ = differ<_differ::reparameterized, A_pole>;
	using W_differ = any_t<A_differ>;

	using Y_differ = process::confined_t<void
	,	typename W_differ::template dispatch<>
	,	A_differ
	>;

	using U_order = typename W_differ::order_type;
	static_assert(same_q<U_order, typename Y_differ::order_type>);
	static_assert(same_q<R_order, typename Y_differ::order_type>);

	TRY_("instantiation")
	{
		Y_differ diff{};
		//\
		diff <<= typename Y_differ::order_type{1};
		diff <<= typename W_differ::order_type{1};

		TRUE_(check_f<19>(diff(1), 1));
		TRUE_(check_f<19>(diff(2), 1));
		TRUE_(check_f<19>(diff(4), 2));
		TRUE_(check_f<19>(diff(8), 4));

	}
//	TRY_("scaling")
//	{
//		confined_t<dilate<[] XTAL_1FN_(to) (bond::fit<>::patio_2)>, Y_differ> diff{};
//		diff <<= typename Y_differ::order_type{1};
//
//		TRUE_(check_f<19>(0.15915493667125702, diff(1)));
//		TRUE_(check_f<19>(0.15915493667125702, diff(2)));
//
//	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
