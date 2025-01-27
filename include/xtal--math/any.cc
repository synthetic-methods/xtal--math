#pragma once
#include "./any.hh"// testing...

#include "./etc.cc"


namespace xtal
{
template <class U, class V=U> struct   complexion;
template <class        ...Ts> concept  complexion_q = bond::tag_p<complexion, Ts...>;

XTAL_DEF_(return,inline,let)
objective_f(complexion_q auto &&o)
noexcept -> decltype(auto)
{
	using X = XTAL_ALL_(o);
	using Y = typename X::target_type;
	return Y{objective_f(o.real()), objective_f(o.imag())};
}
XTAL_DEF_(return,inline,let)
injective_f(auto &&o)
noexcept -> decltype(auto)
requires requires {o.real(); o.imag();}
{
	return _std::make_pair(objective_f(o.real()), objective_f(o.imag()));
}

template <class U, class V>
struct complexion
{
	using        type = complexion;
	using  value_type = U;
	using  valve_type = V;
	using source_type = _std::complex<U>;
	using target_type = _std::complex<V>;

	using _fix = xtal::bond::fixture<value_type>;

	value_type re;
	value_type im;

	XTAL_DEF_(return,inline) auto real() const noexcept -> value_type const & {return re;}
	XTAL_DEF_(return,inline) auto imag() const noexcept -> value_type const & {return im;}

	XTAL_DEF_(return,inline) auto real(auto &&...oo) const noexcept -> value_type const & {return re = value_type{XTAL_REF_(oo)...};}
	XTAL_DEF_(return,inline) auto imag(auto &&...oo) const noexcept -> value_type const & {return im = value_type{XTAL_REF_(oo)...};}

	
	template <class T> requires un_n<isotropic_q<T, source_type>> and un_n<isotropic_q<T, target_type>>
	XTAL_DEF_(return,inline,friend,let)
	operator/ (T const &t, source_type const &s)
	noexcept -> target_type
	requires un_n<requires {t.real(); t.imag();}>
	{
		auto const [s_re, s_im] = injective_f(s);
		auto const s_abs = t/(s_re*s_re + s_im*s_im);
		return {s_abs*s_re, -s_abs*s_im};
	}
	XTAL_DEF_(return,inline,friend,let)
	operator/ (source_type const &s, source_type const &t)
	noexcept -> target_type
	{
		return s*(one/t);
	}
	XTAL_DEF_(return,inline,friend,let)
	operator/ (source_type const &s, target_type const &t)
	noexcept -> target_type
	requires un_n<isotropic_q<source_type, target_type>>
	{
		return s*(one/t);
	}
	XTAL_DEF_(return,inline,friend,let)
	operator/ (source_type const &s, auto const &t)
	noexcept -> target_type
	{
		return s*(one/t);
	}
	
//	Complex multiplication:
	XTAL_DEF_(return,inline,friend,let)
	operator* (source_type const &s, source_type const &t)
	noexcept -> target_type
	{
		auto const [s_re, s_im] = injective_f(s);
		auto const [t_re, t_im] = injective_f(t);
		return {s_re*t_re - s_im*t_im, s_im*t_re + s_re*t_im};
	}
	XTAL_DEF_(return,inline,friend,let)
	operator* (source_type const &s, target_type const &t)
	noexcept -> target_type
	requires un_n<isotropic_q<source_type, target_type>>
	{
		auto const [s_re, s_im] = injective_f(s);
		auto const [t_re, t_im] = injective_f(t);
		return {s_re*t_re - s_im*t_im, s_im*t_re + s_re*t_im};
	}
	template <class T> requires un_n<isotropic_q<T, source_type>> and un_n<isotropic_q<T, target_type>>
	XTAL_DEF_(return,inline,friend,let)
	operator* (source_type const &s, T const &t)
	noexcept -> target_type
	requires in_n<requires {t.real()*t.imag();}>
	{
		auto const [s_re, s_im] = injective_f(s);
		auto const [t_re, t_im] = injective_f(t);
		return {s_re*t_re - s_im*t_im, s_im*t_re + s_re*t_im};
	}
//	Scalar multiplication:
	template <class T> requires un_n<isotropic_q<T, source_type>> and un_n<isotropic_q<T, target_type>>
	XTAL_DEF_(return,inline,friend,let)
	operator* (source_type const &s, T const &t)
	noexcept -> target_type
	requires un_n<requires {t.real()*t.imag();}>
	{
		return {s.re*t, s.im*t};
	}
//	Associative multiplication:
	template <class T> requires un_n<isotropic_q<T, source_type>>
	XTAL_DEF_(return,inline,friend,let)
	operator* (T const &t, source_type const &s)
	noexcept -> auto
	{
		return s*t;
	}

//	Complex addition:
	XTAL_DEF_(return,inline,friend,let)
	operator+ (source_type const &s, source_type const &t)
	noexcept -> target_type
	{
		auto const [s_re, s_im] = reinterpret_cast<value_type const(&)[2]>(s);
		auto const [t_re, t_im] = reinterpret_cast<value_type const(&)[2]>(t);
		return {s_re + t_re, s_re + t_im};
	}
	XTAL_DEF_(return,inline,friend,let)
	operator+ (source_type const &s, target_type const &t)
	noexcept -> target_type
	requires un_n<isotropic_q<source_type, target_type>>
	{
		auto const [s_re, s_im] = reinterpret_cast<value_type const(&)[2]>(s);
		auto const [t_re, t_im] = reinterpret_cast<valve_type const(&)[2]>(t);
		return {s_re + t_re, s_re + t_im};
	}
	template <class T> requires un_n<isotropic_q<T, source_type>> and un_n<isotropic_q<T, target_type>>
	XTAL_DEF_(return,inline,friend,let)
	operator+ (source_type const &s, T const &t)
	noexcept -> target_type
	requires in_n<requires {t.real() + t.imag();}>
	{
		auto const [s_re, s_im] = reinterpret_cast<value_type const(&)[2]>(s);
		auto const [t_re, t_im] = injective_f(t);
		return {s_re + t_re, s_re + t_im};
	}
//	Scalar addition:
	template <class T> requires un_n<isotropic_q<T, source_type>> and un_n<isotropic_q<T, target_type>>
	XTAL_DEF_(return,inline,friend,let)
	operator+ (source_type const &s, T const &t)
	noexcept -> target_type
	requires un_n<requires {t.real() + t.imag();}>
	{
		return {s.re+t, s.im};
	}
//	Associative addition:
	template <class T> requires un_n<isotropic_q<T, source_type>>
	XTAL_DEF_(return,inline,friend,let)
	operator+ (T const &t, source_type const &s)
	noexcept -> auto
	{
		return s + t;
	}

	template <class T>
	XTAL_DEF_(return,inline,friend,let)
	operator- (source_type const &s)
	noexcept -> target_type
	{
		return {-s.re, -s.im};
	}

};

}
namespace std
{
XTAL_DEF_(return,inline)
auto exp(xtal::complexion_q auto &&x)
{
	using X = XTAL_ALL_(x);
	using Y = typename  X::target_type;
	//\
	return Y{cos(x.imag()), sin(x.imag())}*exp(x.real());
	return exp(x.real())*Y{cos(x.imag()), sin(x.imag())};
}
XTAL_DEF_(return,inline)
auto conj(xtal::complexion_q auto &&x)
{
	using X = XTAL_ALL_(x);
	using Y = typename  X::target_type;
	return Y{x.re, -x.im};
}


template <> struct complex<           Eigen::ArrayXd > : xtal::complexion<           Eigen::ArrayXd , Eigen::ArrayXd> {};
template <> struct complex<Eigen::Map<Eigen::ArrayXd>> : xtal::complexion<Eigen::Map<Eigen::ArrayXd>, Eigen::ArrayXd> {};

}


XTAL_ENV_(push)
namespace xtal::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <int N_index>
XTAL_DEF_(return,let)
check_f(auto const &u, auto const &v)
noexcept -> bool
{
	using _fix = bond::fixture<decltype(u), decltype(v)>;
	return _fix::template trim_f<N_index>(u) == _fix::template trim_f<N_index>(v);
}
template <int N_index, int N_limit>
XTAL_DEF_(return,let)
check_f(auto const &u, auto const &v)
noexcept -> int
{
	auto constexpr Z_index = sign_n<N_index>;
	auto constexpr Z_limit = sign_n<N_limit>;
	static_assert(Z_index == Z_limit);

	XTAL_IF0
	XTAL_0IF (Z_limit*N_limit <  Z_index*N_index) {
		return check_f<N_limit, N_index>(u, v);
	}
	XTAL_0IF (Z_limit*N_limit == Z_index*N_index) {
		return 0;
	}
	XTAL_0IF (Z_index == -1) {
		return check_f<N_index>(u, v)? N_index: check_f<Z_index + N_index, N_limit>(u, v);
	}
	XTAL_0IF (Z_index ==  1) {
		return check_f<N_limit>(u, v)? N_limit: check_f<N_index, N_limit - Z_limit>(u, v);
	}
}
XTAL_DEF_(return,let)
check_f(auto const &u, auto const &v)
noexcept -> int
{
	return check_f<-1, 1 - (int) bond::fixture<>::fraction.depth>(u, v);
}


////////////////////////////////////////////////////////////////////////////////
/**/
template <int N, int I>
auto plot_f(int const i)
{
	int constexpr Z = sizeof(int) << 2;

	int constexpr  I_pos = 0 < I, I_neg = I < 0, I_abs = (I^-I_neg) + I_neg;
	int const      i_pos = 0 < i, i_neg = i < 0, i_abs = (i^-i_neg) + i_neg;
	int const      q_pos = i_pos == I_pos;
	int const      q_neg = i_neg == I_neg;

	int k = i_pos - i_neg;

	XTAL_IF0
	XTAL_0IF (I == 0) {            k += 2;}
	XTAL_0IF (I != 0) {if (I != i) k  = 2;}

	k &= -(I_abs <= i_abs);
	k &= -(q_pos or q_neg);
	k &= 0b11;

	_std::array<char, 4> constexpr table{' ', ',', '|', '\''};
	return table[k];
}
template <int N, int I>
auto plot_f(double i)
{
	return plot_f<N, I>(static_cast<int>(round(i)));
}
template <int N, int M=N - 1>
void plot(iterated_q auto const o)
{
	int constexpr L = N*2 + 1;
	bond::seek_backward_f<L>([&] (auto i) XTAL_0FN {
		for (auto &&e:o) {
			_std::cout << plot_f<N, i - N>(M*e);
		};
		_std::cout << _std::endl;
	});
}
/***/

////////////////////////////////////////////////////////////////////////////////
/*/
TAG_("any")
{
	using namespace Eigen;

	using Y_alpha =     ArrayXd ;
	using X_alpha = Map<ArrayXd>;
	
	using Y_aphex = _std::complex<Y_alpha>;
	using X_aphex = _std::complex<X_alpha>;

	static_assert(same_q<X_aphex, typename X_aphex::source_type>);
	static_assert(same_q<Y_aphex, typename X_aphex::target_type>);

	TRY_("task")
	{
		double foo[2][2] {{1, 11}, {2, 22}};
		double goo[2][2] {{3, 33}, {4, 44}};

		X_aphex bar{X_alpha(foo[0], 2), X_alpha(foo[1], 2)};
		X_aphex car{X_alpha(goo[0], 2), X_alpha(goo[1], 2)};
		//\
		Y_aphex baz = car*car + bar;
		Y_aphex baz = _std::exp(bar) * car;

	//	echo(baz.real(0), baz.real(1));
		TRUE_(true);

	}
}
/***/

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
