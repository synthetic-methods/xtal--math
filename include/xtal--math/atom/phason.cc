#pragma once
#include "./any.cc"
#include "./phason.hh"// testing...





XTAL_ENV_(push)
namespace xtal::atom::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

//atic_assert(_xtd::trivially_initializable<phason_t<float[2]>>);
static_assert(_xtd::trivially_destructible <phason_t<float[2]>>);
static_assert(_xtd::trivially_copyable     <phason_t<float[2]>>);
static_assert(_xtd::trivially_movable      <phason_t<float[2]>>);
//atic_assert(                     atomic_q<phason_t<float[2]>>);


////////////////////////////////////////////////////////////////////////////////

TAG_("phason")
{
	using V_fit = bond::fit<int>;
	using W_fit = bond::fit<double>;
	
	using V_delta = typename V_fit::delta_type;
	using W_delta = typename W_fit::delta_type;
	using V_sigma = typename V_fit::sigma_type;
	using W_sigma = typename W_fit::sigma_type;
	using V_alpha = typename V_fit::alpha_type;
	using W_alpha = typename W_fit::alpha_type;
	using V_aphex = typename V_fit::aphex_type;
	using W_aphex = typename W_fit::aphex_type;

	static constexpr W_alpha two =  2;
	static constexpr W_alpha six =  6;
	static constexpr W_alpha ten = 10;

	auto mt19937_f = W_fit::mt19937_t(Catch::rngSeed());

	using V_phi = W_alpha;
	using U_phi = phason_t<V_phi[2]>;
	using W_phi = _std::complex<phason_t<W_alpha[2]>>;
	using A_phi = _std::array<V_phi, 2>;

	using _qp = bond::template fit<typename U_phi::value_type>;

	using D1 = phason_t<W_alpha[1]>;
	using D2 = phason_t<W_alpha[2]>;
	using D3 = phason_t<W_alpha[3]>;
	using D4 = phason_t<W_alpha[4]>;
	
	TRY_("phason scaling")
	{
		U_phi u_phi{0.2, 0.1}, v_phi = u_phi.scaled(0.3);
		u_phi.scale(0.3);

		TRUE_(u_phi == v_phi);
		TRUE_(check_f<-20>(u_phi(1), 0.03));
		TRUE_(check_f<-20>(u_phi(0), 0.20));
		//\
		u_phi +=      {0.1};
		u_phi += U_phi{0.1};
		TRUE_(check_f<-20>(u_phi(0), 0.30));

	//	u_phi = U_phi{ 0.5, 0.5*W_fit::dnsilon_f(44)}; TRUE_(u_phi == u_phi);
	//	u_phi = U_phi{-0.5, 0.5*W_fit::dnsilon_f(44)}; TRUE_(u_phi == u_phi);

	}
	/**/
	TRY_("phason of complex")
	{
		using M_phi = phason_t<_std::complex<W_alpha>[2]>;
		phason_t<_std::complex<W_alpha>[2]> a{{0.1, 0.2}, {0.3, 0.3}};
		phason_t<_std::complex<W_alpha>[2]> b{{0.1, 0.2}, {0.3, 0.3}};

		++a;
		++a;
		++a;
		TRUE_(check_f<12>(a(0), W_aphex{0.0, 0.1}));

	//	W_phi z = a;
	//	TRUE_(check_f<19>(z.real() (0), 0.0));
	//	TRUE_(check_f<19>(z.imag() (0), 0.1));
	//	TRUE_(check_f<19>(z.real() (1), 0.3));
	//	TRUE_(check_f<19>(z.imag() (1), 0.3));

	}
	/***/
	/**/
	TRY_("real of phason")
	{
		U_phi a{0.1, 0.2};
		U_phi b{0.1, 0.2};
		//\
		U_phi c = 2.0 * a;
		U_phi c = b + a;

		TRUE_(check_f(c(0), 0.2));

	}
	/***/
	/**/
	TRY_("complex of phason")
	{
		W_phi a{{0.1, 0.0}, {0.2, 0.0}};
		W_phi b{{0.1, 0.0}, {0.2, 0.0}};
		//\
		W_phi c = 2.0 * a;
		W_phi c = b + a;

		TRUE_(check_f(a.real() (0), 0.1));
		TRUE_(check_f(b.real() (0), 0.1));
		TRUE_(check_f(c.real() (0), 0.2));

		TRUE_(check_f(a.imag() (0), 0.2));
		TRUE_(check_f(b.imag() (0), 0.2));
		TRUE_(check_f(c.imag() (0), 0.4));

	}
	TRY_("tuple with phason")
	{
		using U_psi = group_addition_t<W_alpha[2]>;
		using W     = group_addition_t<U_phi, U_psi>;

		TRUE_(sizeof(W) == sizeof(U_phi) + sizeof(U_psi));

		W w0{{0.1, 0.2}, {0.3, 0.4}};
		W w1{{0.1, 0.2}, {0.3, 0.4}};
		w0 *= 2.0;
		TRUE_(check_f<-22>( 0.2, get<0>(w0)(0)));
		TRUE_(check_f<-22>( 0.4, get<0>(w0)(1)));
		TRUE_(check_f<- 1>( 0.6, get<1>(w0)(0)));
		TRUE_(check_f<- 1>( 0.8, get<1>(w0)(1)));
		w0 += w1;
		TRUE_(check_f<-22>( 0.3, get<0>(w0)(0)));
		TRUE_(check_f<-22>(-0.4, get<0>(w0)(1)));
		TRUE_(check_f<- 1>( 0.9, get<1>(w0)(0)));
		TRUE_(check_f<- 1>( 1.2, get<1>(w0)(1)));

	}
	/***/
	/**/
	TRY_("construction")
	{
		D2 a_d2{0.250, 0.250};
		D2 b_d2{0.125, 0.125};
		D2 y_d2{0.000, 0.125};
		//\
		D2 z_d2(0.125);
		D2 z_d2{0.125};

		TRUE_(z_d2 == D2{0.125, 0.000});

		TRUE_(y_d2 == D2{0.000, 0.125}); y_d2 <<= {0.250};
		TRUE_(y_d2 == D2{0.000, 0.250}); y_d2 <<= _std::array<W_alpha, 1>{0.500};
	//	TRUE_(y_d2 == D2{0.000, 0.500}); y_d2 >>= _std::array<W_alpha, 1>{0.333};
	//	TRUE_(y_d2 == D2{0.333, 0.500});

		y_d2 = D2{0.333, 0.500};
		TRUE_(check_f<8>((2.0*y_d2)(0), -0.334));
		TRUE_(check_f<8>((2.0*y_d2)(1),  0.000));
		
		TRUE_(a_d2 == D2{0.250, 0.250});

		a_d2 = b_d2;
		TRUE_(a_d2 == b_d2); b_d2 = a_d2; ++a_d2[0];
		TRUE_(a_d2 != b_d2); b_d2 = a_d2; ++a_d2[0];
		TRUE_(a_d2 != b_d2); b_d2 = a_d2; ++a_d2[0];
		TRUE_(a_d2 != b_d2); b_d2 = a_d2; ++a_d2[0];
	//	TRUE_(a_d2 == D2{0x40000000, 0x40000000});

		TRUE_(check_f<4>(0.1, phason_t<W_alpha[2]>{0.1, 0.1}(0)));

	}
	/***/
	/**/
	TRY_("iteration")
	{
		U_phi   phi{0, W_fit::haplo_f(4)};
		auto constexpr N_tip = _qp::sign.mask;
		auto           n_tip = N_tip;

		n_tip ^= N_tip;
		for (int i = 0; i < 0x8; ++i) {
			TRUE_((N_tip&phi[0]) == n_tip);
			++phi;
		}
		n_tip ^= N_tip;
		for (int i = 0; i < 0x8; ++i) {
			TRUE_((N_tip&phi[0]) == n_tip);
			++phi;
		}
		n_tip ^= N_tip;
		for (int i = 0; i < 0x8; ++i) {
			TRUE_((N_tip&phi[0]) == n_tip);
			++phi;
		}
		n_tip ^= N_tip;
		for (int i = 0; i < 0x8; ++i) {
			TRUE_((N_tip&phi[0]) == n_tip);
			++phi;
		}

	}
	/***/
	/**/
	TRY_("addition")
	{
		W_alpha x =  0.33, x_dt = W_fit::haplo_f(4);
		W_alpha y =  5.55;
		W_alpha z =  x+y; z -= _std::round(z);

		U_phi phi{x, x_dt};

		phi += y;
	//	TRUE_(check_f<8>(phi[0], z));
		TRUE_(check_f<8>(phi(0), z));

	}
	/***/
	/**/
	TRY_("multiplication")
	{
		W_alpha x = 1/3.L, x_dt = W_fit::haplo_f(4);
		U_phi phi{x, x_dt};
		V_phi foo{x};

		for (W_sigma i = 0; i < 0x4; ++i) {
			auto const v = static_cast<W_alpha>(i + 1);
			auto const u{exp(asinh(v*W_fit::haplo_1))};
			phi *= u;
			foo *= u;
			foo -= _std::round(foo);
			TRUE_(check_f<8>(phi(0), foo));
		}

	}
	/***/
}
TAG_("phason trials")
{
	using V_fit = bond::fit<int>;
	using W_fit = bond::fit<double>;
	using W_delta = typename W_fit::delta_type;
	using W_sigma = typename W_fit::sigma_type;
	using W_alpha = typename W_fit::alpha_type;
	using W_aphex = typename W_fit::aphex_type;
	static constexpr W_alpha two =  2;
	static constexpr W_alpha six =  6;
	static constexpr W_alpha ten = 10;

	auto mt19937_o = typename W_fit::mt19937_t{}; mt19937_o.seed(Catch::rngSeed());
	auto mt19937_f = [&] XTAL_1FN_(to) (W_fit::mantissa_f(mt19937_o));

	using V_phi = W_alpha;
	using U_phi = phason_t<V_phi[2]>;
	using W_phi = _std::complex<phason_t<W_alpha[2]>>;
	using A_phi = _std::array<V_phi, 2>;

	using _qp = bond::template fit<typename U_phi::value_type>;

	using D1 = phason_t<W_alpha[1]>;
	using D2 = phason_t<W_alpha[2]>;
	using D3 = phason_t<W_alpha[3]>;
	using D4 = phason_t<W_alpha[4]>;
	
	EST_("atom::phason_t\n   #1*#2&\n   (*fixed-point*)")
	{
		W_alpha x = 0.33, x_dt = W_fit::haplo_f(4);

		U_phi   phi{x, x_dt};
	//	W_alpha foo{x};

		for (W_sigma i = 0x100; ~--i;) {
			W_alpha const u = _std::pow(two, 1 + mt19937_f());
			phi *= u;
		}
		return phi;
	};
	EST_("atom::phason_t\n   #1*#2&\n   (*floating-point*)")
	{
		W_alpha x = 0.33, x_dt = W_fit::haplo_f(4);

	//	U_phi   phi{x, x_dt};
		W_alpha foo{x};

		for (W_sigma i = 0x100; ~--i;) {
			W_alpha const u = _std::pow(two, 1 + mt19937_f());
			foo *= u;
			foo -= _std::round(foo);
		}
		return foo;
	};
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
