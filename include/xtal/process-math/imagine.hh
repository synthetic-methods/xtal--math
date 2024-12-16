#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///\
Applies complex quarter-rotation and conjugation. \

template <int M_rot=0, int M_con=0> struct   imagine;
template <int M_rot=0, int M_con=0> using    imagine_t = process::confined_t<imagine<M_rot, M_con>>;
template <int M_rot=0, int M_con=0>
XTAL_DEF_(short)
XTAL_LET imagine_f(auto &&o)
noexcept -> decltype(auto)
{
	return imagine_t<M_rot, M_con>::function(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int M_rot, int M_con>
struct imagine
{
	static constexpr size_type N_rot = M_rot&0b11U;
	static constexpr size_type N_con = M_con&0b01U;
	static constexpr size_type X_con = N_con << 1U;
	static constexpr size_type X_rot = N_rot^X_con;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(short,static)
		XTAL_LET function(auto &&o)
		noexcept -> decltype(auto)
			requires un_q<complex_field_q<decltype(o)>>
		{
			if constexpr (N_rot == 0 and N_con == 0) {
				return XTAL_REF_(o);
			}
			else {
				return function(complexion_f(o));
			}
		};
		template <auto ...>
		XTAL_DEF_(short,static)
		XTAL_LET function(complex_field_q auto const &o)
		noexcept -> decltype(auto)
			requires un_q<complex_number_q<decltype(o)>>
		{
			auto const x = o.real();
			auto const y = o.imag();
			XTAL_IF0

			XTAL_0IF (N_rot == 0b00 and N_con == 0) {return                    o;}
			XTAL_0IF (N_rot == 0b00 and N_con == 1) {return complexion_f( x, -y);}

			XTAL_0IF (N_rot == 0b01 and N_con == 0) {return complexion_f(-y,  x);}
			XTAL_0IF (N_rot == 0b01 and N_con == 1) {return complexion_f(-y, -x);}

			XTAL_0IF (N_rot == 0b10 and N_con == 0) {return complexion_f(-x, -y);}
			XTAL_0IF (N_rot == 0b10 and N_con == 1) {return complexion_f(-x,  y);}

			XTAL_0IF (N_rot == 0b11 and N_con == 0) {return complexion_f( y, -x);}
			XTAL_0IF (N_rot == 0b11 and N_con == 1) {return complexion_f( y,  x);}
		};
		template <auto ...>
		XTAL_DEF_(short,static)
		XTAL_LET function(complex_number_q auto o)
		noexcept -> decltype(auto)
		{
			auto &[x, y] = destruct_f(o);
			if constexpr (N_rot == 0b01 or N_rot == 0b11) {_std::swap(x, y);}
			if constexpr (N_rot == 0b01 or N_rot == 0b10) {x = -x;}
			if constexpr (X_rot == 0b11 or X_rot == 0b10) {y = -y;}
			return o;
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
