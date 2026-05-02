#pragma once
#include "./any.cc"





#include "./imagine.hh"
XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/**/
TAG_("imagine")
{
	using U_fit   = bond::fit<>;
	using U_sigma = typename U_fit::sigma_type;
	using U_delta = typename U_fit::delta_type;
	using U_alpha = typename U_fit::alpha_type;
	using U_aphex = typename U_fit::aphex_type;

	TRY_("imagine `complex_field_q auto const &`")
	{
		U_alpha x{1};
		U_alpha y{2};
		U_aphex o{x, y};

		TRUE_(imagine_f<0, 0>(o) == U_aphex{ x,  y});
		TRUE_(imagine_f<0, 1>(o) == U_aphex{ x, -y});

		TRUE_(imagine_f<1, 0>(o) == U_aphex{-y,  x});
		TRUE_(imagine_f<1, 1>(o) == U_aphex{-y, -x});

		TRUE_(imagine_f<2, 0>(o) == U_aphex{-x, -y});
		TRUE_(imagine_f<2, 1>(o) == U_aphex{-x,  y});

		TRUE_(imagine_f<3, 0>(o) == U_aphex{ y, -x});
		TRUE_(imagine_f<3, 1>(o) == U_aphex{ y,  x});

	}
	TRY_("imagine `complex_variable_q auto &&`")
	{
		U_alpha x{1};
		U_alpha y{2};

		TRUE_(imagine_f<0, 0>(U_aphex{x, y}) == U_aphex{ x,  y});
		TRUE_(imagine_f<0, 1>(U_aphex{x, y}) == U_aphex{ x, -y});

		TRUE_(imagine_f<1, 0>(U_aphex{x, y}) == U_aphex{-y,  x});
		TRUE_(imagine_f<1, 1>(U_aphex{x, y}) == U_aphex{-y, -x});

		TRUE_(imagine_f<2, 0>(U_aphex{x, y}) == U_aphex{-x, -y});
		TRUE_(imagine_f<2, 1>(U_aphex{x, y}) == U_aphex{-x,  y});

		TRUE_(imagine_f<3, 0>(U_aphex{x, y}) == U_aphex{ y, -x});
		TRUE_(imagine_f<3, 1>(U_aphex{x, y}) == U_aphex{ y,  x});

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
