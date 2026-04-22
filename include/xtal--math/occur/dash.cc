#pragma once
#include "./any.cc"





#include "./dash.hh"// testing...
XTAL_ENV_(push)
namespace xtal::occur::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

TAG_("dash")
{
	TRY_("construct")
	{
		using U_source = flow::conferred_t<counted_t<>>;
		using U_target = dash_s<U_source>;
		using V_target = dash_s<>;

		V_target v_target(99);
		U_target u_target(99, U_source(counted_t<>(11, 22)));

		TRUE_(99 == v_target.template head<constant_t<0>>());
		TRUE_(99 == u_target.template head<constant_t<0>>());
		TRUE_(equal_f(counted_t<>(11, 22), u_target.tail()));

		TRUE_(same_q<decltype(XTAL_ANY_(U_target).head()), typename V_target::head_type>);
		TRUE_(same_q<decltype(XTAL_ANY_(U_target).tail()), U_source>);

	}
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
