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
XTAL_DEF_(return,inline,let)
imagine_f(auto &&o)
noexcept -> decltype(auto)
{
	return imagine_t<M_rot, M_con>::method_f(XTAL_REF_(o));
}


////////////////////////////////////////////////////////////////////////////////

template <int M_rot, int M_con>
struct imagine
{
	static constexpr unsigned N_rot = M_rot&0b11U;
	static constexpr unsigned N_con = M_con&0b01U;
	static constexpr unsigned X_con = N_con << 1U;
	static constexpr unsigned X_rot = N_rot^X_con;

	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:
		using S_::S_;

		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(auto &&o)
		noexcept -> decltype(auto)
			requires un_n<complex_field_q<decltype(o)>>
		{
			if constexpr (N_rot == 0 and N_con == 0) {
				return XTAL_REF_(o);
			}
			else {
				return method_f(complexion_f(o));
			}
		};
		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(complex_field_q auto const &o)
		noexcept -> decltype(auto)
			requires un_n<complex_variable_q<decltype(o)>>
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
		XTAL_DEF_(return,inline,set)
		method_f(complex_variable_q auto o)
		noexcept -> decltype(auto)
		{
			auto &[x, y] = destruct_f(o);
			if constexpr (N_rot == 0b01 or N_rot == 0b11) {_std::swap(x, y);}
			if constexpr (N_rot == 0b01 or N_rot == 0b10) {x = -XTAL_MOV_(x);}
			if constexpr (X_rot == 0b11 or X_rot == 0b10) {y = -XTAL_MOV_(y);}
			return o;
		};

	};
};


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
