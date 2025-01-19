#pragma once
#include "./any.cc"
#include "./imagine.hh"// testing...





XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("imagine")
{
	using _fix = bond::fixture<>;
	using T_sigma = typename _fix::sigma_type;
	using T_delta = typename _fix::delta_type;
	using T_alpha = typename _fix::alpha_type;
	using T_aphex = typename _fix::aphex_type;

	TRY_("imagine `complex_field_q auto const &`")
	{
		T_alpha x{1};
		T_alpha y{2};
		T_aphex o{x, y};

		TRUE_(imagine_f<0, 0>(o) == T_aphex{ x,  y});
		TRUE_(imagine_f<0, 1>(o) == T_aphex{ x, -y});

		TRUE_(imagine_f<1, 0>(o) == T_aphex{-y,  x});
		TRUE_(imagine_f<1, 1>(o) == T_aphex{-y, -x});

		TRUE_(imagine_f<2, 0>(o) == T_aphex{-x, -y});
		TRUE_(imagine_f<2, 1>(o) == T_aphex{-x,  y});

		TRUE_(imagine_f<3, 0>(o) == T_aphex{ y, -x});
		TRUE_(imagine_f<3, 1>(o) == T_aphex{ y,  x});

	}
	TRY_("imagine `complex_variable_q auto &&`")
	{
		T_alpha x{1};
		T_alpha y{2};

		TRUE_(imagine_f<0, 0>(T_aphex{x, y}) == T_aphex{ x,  y});
		TRUE_(imagine_f<0, 1>(T_aphex{x, y}) == T_aphex{ x, -y});

		TRUE_(imagine_f<1, 0>(T_aphex{x, y}) == T_aphex{-y,  x});
		TRUE_(imagine_f<1, 1>(T_aphex{x, y}) == T_aphex{-y, -x});

		TRUE_(imagine_f<2, 0>(T_aphex{x, y}) == T_aphex{-x, -y});
		TRUE_(imagine_f<2, 1>(T_aphex{x, y}) == T_aphex{-x,  y});

		TRUE_(imagine_f<3, 0>(T_aphex{x, y}) == T_aphex{ y, -x});
		TRUE_(imagine_f<3, 1>(T_aphex{x, y}) == T_aphex{ y,  x});

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
