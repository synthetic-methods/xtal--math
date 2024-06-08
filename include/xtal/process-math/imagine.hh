#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int M_rot=0, int M_con=0> XTAL_TYP imagine;
template <int M_rot=0, int M_con=0> XTAL_USE imagine_t = process::confined_t<imagine<M_rot, M_con>>;
template <int M_rot=0, int M_con=0>
XTAL_FN2 imagine_f(complex_field_q auto &&o)
XTAL_0EX
{
	return imagine_t<M_rot, M_con>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int M_rot, int M_con>
struct imagine
{
	XTAL_LET N_rot = M_rot&0b11;
	XTAL_LET N_con = M_con&1;

	template <class S>
	class subtype: public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_FN2 function(complex_number_q auto &&o)
		XTAL_0EX
		{
			auto &xy = devolved_f(XTAL_REF_(o));
			auto const &x = xy[0];
			auto const &y = xy[1];
			XTAL_IF0

			XTAL_0IF (N_rot == 0b00 and N_con == 0) {return complexion_f( x,  y);}
			XTAL_0IF (N_rot == 0b00 and N_con == 1) {return complexion_f( x, -y);}

			XTAL_0IF (N_rot == 0b01 and N_con == 0) {return complexion_f(-y,  x);}
			XTAL_0IF (N_rot == 0b01 and N_con == 1) {return complexion_f(-y, -x);}

			XTAL_0IF (N_rot == 0b10 and N_con == 0) {return complexion_f(-x, -y);}
			XTAL_0IF (N_rot == 0b10 and N_con == 1) {return complexion_f(-x,  y);}

			XTAL_0IF (N_rot == 0b11 and N_con == 0) {return complexion_f( y, -x);}
			XTAL_0IF (N_rot == 0b11 and N_con == 1) {return complexion_f( y,  x);}
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
