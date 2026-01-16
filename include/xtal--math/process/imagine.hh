#pragma once
#include "./any.hh"

//#import "../atom/pade/uniplex.hh"




XTAL_ENV_(push)
namespace xtal::process::math
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Applies complex quarter-rotation and conjugation.
*/
template <int M_rot=0, int M_con=0>
struct  imagine;


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
		//\
		method_f(atom::math::pade::uniplex_q auto &&z)
		method_f(atom::couple_q auto &&z)
		noexcept -> decltype(auto)
		requires complex_variable_q<decltype(get<0>(z))> and atom::couple_q<decltype(get<1>(z))>
		{
			return z.form(method_f(z.signum()), z.magnum());
		};
		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(complex_field_q auto &&z)
		noexcept -> auto
		{
			if constexpr (complex_variable_q<decltype(z)>) {
				if constexpr (N_rot&1) {
					auto o = XTAL_REF_(z); auto &[x, y] = destruct_f(o);
					if constexpr (N_rot == 0b01 or N_rot == 0b11) {_std::swap(x, y);}
					if constexpr (N_rot == 0b01 or N_rot == 0b10) {x = -XTAL_MOV_(x);}
					if constexpr (X_rot == 0b11 or X_rot == 0b10) {y = -XTAL_MOV_(y);}
					return o;
				}
				else {
					int constexpr N_ron = N_rot|N_con;
					if constexpr (N_ron == 0b0'0) {return      (XTAL_REF_(z));}
					if constexpr (N_ron == 0b1'0) {return -    (XTAL_REF_(z));}
					if constexpr (N_ron == 0b0'1) {return  conj(XTAL_REF_(z));}
					if constexpr (N_ron == 0b1'1) {return -conj(XTAL_REF_(z));}
				}
			}
			else {
				return method_f(z.real(), z.imag());
			}
		};
		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(simplex_field_q auto &&x)
		noexcept -> decltype(auto)
		{
			return method_f(XTAL_REF_(x), XTAL_ALL_(x){});
		};
		template <auto ...>
		XTAL_DEF_(return,inline,set)
		method_f(simplex_field_q auto &&x, simplex_field_q auto &&y)
		noexcept -> decltype(auto)
		{
			XTAL_IF0

			XTAL_0IF (N_rot == 0b00 and N_con == 0) {return complexion_f( XTAL_REF_(x),  XTAL_REF_(y));}
			XTAL_0IF (N_rot == 0b00 and N_con == 1) {return complexion_f( XTAL_REF_(x), -XTAL_REF_(y));}

			XTAL_0IF (N_rot == 0b01 and N_con == 0) {return complexion_f(-XTAL_REF_(y),  XTAL_REF_(x));}
			XTAL_0IF (N_rot == 0b01 and N_con == 1) {return complexion_f(-XTAL_REF_(y), -XTAL_REF_(x));}

			XTAL_0IF (N_rot == 0b10 and N_con == 0) {return complexion_f(-XTAL_REF_(x), -XTAL_REF_(y));}
			XTAL_0IF (N_rot == 0b10 and N_con == 1) {return complexion_f(-XTAL_REF_(x),  XTAL_REF_(y));}

			XTAL_0IF (N_rot == 0b11 and N_con == 0) {return complexion_f( XTAL_REF_(y), -XTAL_REF_(x));}
			XTAL_0IF (N_rot == 0b11 and N_con == 1) {return complexion_f( XTAL_REF_(y),  XTAL_REF_(x));}
		};

	};
};


////////////////////////////////////////////////////////////////////////////////

template <int M_rot=0, int M_con=0>
XTAL_TYP_(let) imagine_t = process::confined_t<imagine<M_rot, M_con>>;

template <int M_rot=0, int M_con=0, auto ...Ns>
XTAL_DEF_(let) imagine_f = [] XTAL_1FN_(call) (imagine_t<M_rot, M_con>::template method_f<Ns...>);


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
