#pragma once
#include "./any.cc"
#include "./phasor.hh"// testing...

#include "../atom/phason.hh"



XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <size_type N_columns>
XTAL_DEF_(return,let)
ConvertToEigenMatrix(double **data, size_type n_rows)
noexcept
{
	using namespace Eigen;
	Array<double, Dynamic, N_columns, RowMajor> m(n_rows, N_columns);
	#pragma unroll
	for (int i = 0; i < N_columns; ++i) {
		//\
		m.col(i) = VectorXd::Map(data[i], n_rows);
		m.col(i) = VectorXd::Map(&data[i][0], n_rows);
	}
	return m;
}

////////////////////////////////////////////////////////////////////////////////

TAG_("phasor")
{
	using namespace Eigen;
	
	using U_stored  = provision::stored<unit_type[0x1000]>;
	using U_sampled = occur::resample_t<>::template attach<>;

	using _Op = bond::fit<>;
	using T_sigma = typename _Op::sigma_type;
	using T_delta = typename _Op::delta_type;
	using T_alpha = typename _Op::alpha_type;
	
	using T_cell = T_alpha;
	using T_eigencolumns = Array<T_cell, Dynamic, 2, ColMajor>;
	using T_eigenrows    = Array<T_cell, Dynamic, 2, RowMajor>;
	using T_eigenrow     = Array<T_cell,       1, 2, RowMajor>;


	using  _phi = T_cell[2];
	using W_phi = bond::repack_t<_phi>;
	using X_phi = atom::math::phason_t<_phi>;
	
	using Y_chi = process::conveyor_t<phasor<_phi, U_sampled>>;
//	using Y_chi = process::conveyor_t<phasor<_phi>>;
	using Y_phi = phasor_t<_phi>;

	//\
	using Y_psi = phasor_t<_phi, U_sampled>;
	using Y_psi = process::lift_t<bond::repack_t<_phi>, phasor<_phi>>;
	//\
	using Y_eig = process::link_t<T_eigenrow, phasor<_phi>>;
	using Y_eig = process::lift_t<T_eigenrow, phasor<_phi>>;

	using Z_chi = processor::monomer_t<Y_chi, U_stored>;
	using Z_phi = processor::monomer_t<Y_phi, U_stored>;
	using Z_psi = processor::monomer_t<Y_psi, U_stored>;
	//\
	using Z_eig = processor::monomer_t<process::lift<bond::operate<T_eigenrow>>, Y_chi>;
//	using Z_eig = processor::monomer_t<confined_t<lift<bond::operate<_std::array<T_cell, 2>>>, phasor<_phi, U_sampled>>>;

	using _fit = bond::template fit<typename X_phi::value_type>;

	/**/
	TRY_("trial")
	{
		TRUE_(sizeof(X_phi) == sizeof(Y_phi));

		static constexpr T_alpha x_delta  = _fit::ratio_f(7);
		
		T_sigma constexpr N_data = 0x1000;
		T_alpha   z_data[2][N_data]{};
		T_alpha  *y_data   [N_data]{z_data[0], z_data[1]};
		T_alpha **x_data = y_data;
		for (int i = 0; i < N_data; ++i) {
			z_data[0][i] =  i;
			z_data[1][i] = -1;
		}


		auto w_data  = bond::transpack_f<void_type[2]>(N_data, z_data);
		//\
		auto e_data  = Map<T_eigencolumns>(*z_data, N_data, 2).rowwise();
		auto e_data  = ConvertToEigenMatrix<2>(y_data, N_data).rowwise();
		
		auto x_phi = X_phi{}; x_phi <<=                                {_fit::ratio_f(7)};
		auto y_phi = Y_phi{}; y_phi <<= occur::math::indent_s<X_phi, 1>{_fit::ratio_f(7)}; y_phi <<= occur::resize_t<>(N_data);
		
		auto z_chi = Z_chi::bind_f(); z_chi <<= occur::math::indent_s<X_phi, 1>{_fit::ratio_f(7)}; z_chi <<= occur::resize_t<>(N_data);
		auto z_phi = Z_phi::bind_f(); z_phi <<= occur::math::indent_s<X_phi, 1>{_fit::ratio_f(7)}; z_phi <<= occur::resize_t<>(N_data);
		auto z_psi = Z_psi::bind_f(); z_psi <<= occur::math::indent_s<X_phi, 1>{_fit::ratio_f(7)}; z_psi <<= occur::resize_t<>(N_data);
	//	auto z_eig = Z_eig::bind_f(); z_eig <<= occur::math::indent_s<X_phi, 1>{_fit::ratio_f(7)}; z_eig <<= occur::resize_t<>(N_data);

		occur::cursor_t<>               z_cursor(N_data);
		occur::math::indent_s<X_phi, 1> z_indent{x_delta};
		
		z_phi <<= z_indent;
		z_psi <<= z_indent;
	//	z_eig <<= z_indent;

		EST_("procession (process in-place)")
		{
			auto &z_d0 = z_data[0];
			auto &z_d1 = z_data[1];
			for (int i = 0; i < N_data; ++i) {
				auto const &y = y_phi();
				z_d0[i] = y(0);
				z_d1[i] = y(1);
			}

		};
		EST_("procession (processor in-place: `ranges::copy...`)")
		{
			//\
			z_psi >>= z_cursor++ >> occur::review_f(w_data);
			z_chi >>= z_cursor++ >> occur::review_f(w_data);

		};
	}
	/***/
	/**/
	TRY_("reprogression")
	{
		using Y_source = phasor_t<_phi>;
		using Y_target = confined_t<phasor<_phi>, phasor<_phi>>;

		Y_source y_source{0, _fit::haplo_f(5)}; X_phi x_source;
		Y_target y_target{0, _fit::haplo_f(4)}; X_phi x_target;
		
		x_source = y_source(); x_target = y_target(x_source, 2.0); TRUE_(check_f<-1>(2.0*x_source(0), x_target(0)));
		x_source = y_source(); x_target = y_target(x_source, 2.0); TRUE_(check_f<-1>(2.0*x_source(0), x_target(0)));
		x_source = y_source(); x_target = y_target(x_source, 2.0); TRUE_(check_f<-1>(2.0*x_source(0), x_target(0)));
		x_source = y_source(); x_target = y_target(x_source, 2.0); TRUE_(check_f<-1>(2.0*x_source(0), x_target(0)));
	//	y_source <<= occur::math::indent_s<X_phi, 0>{0.1};
	//	x_source = y_source(); x_target = y_target(x_source, 2.0); TRUE_(check_f<-18>(2.0*x_source(0), x_target(0)));
	//	x_source = y_source(); x_target = y_target(x_source, 2.0); TRUE_(check_f<-18>(2.0*x_source(0), x_target(0)));
	//	x_source = y_source(); x_target = y_target(x_source, 2.0); TRUE_(check_f<-18>(2.0*x_source(0), x_target(0)));
	//	x_source = y_source(); x_target = y_target(x_source, 2.0); TRUE_(check_f<-18>(2.0*x_source(0), x_target(0)));
	}
	/***/
	/**/
	TRY_("progression")
	{
		T_alpha x_d4 = _fit::haplo_f(4);
		T_alpha x_d3 = _fit::haplo_f(3);
		Y_phi y_phi{0, x_d4};
		X_phi x_phi;

		TRUE_(y_phi() == X_phi{ 1*x_d4, x_d4});
		TRUE_(y_phi() == X_phi{ 2*x_d4, x_d4});
		TRUE_(y_phi() == X_phi{ 3*x_d4, x_d4});
		TRUE_(y_phi() == X_phi{ 4*x_d4, x_d4});
		TRUE_(y_phi() == X_phi{ 5*x_d4, x_d4});
		TRUE_(y_phi() == X_phi{ 6*x_d4, x_d4});
		TRUE_(y_phi() == X_phi{ 7*x_d4, x_d4});
		y_phi <<= occur::math::indent_s<X_phi, 1>{x_d3};
	//	TRUE_(y_phi() == X_phi{-8*x_d4, x_d3});
		TRUE_(y_phi() == X_phi{-7*x_d4, x_d3});
	//	TRUE_(y_phi() == X_phi{-6*x_d4, x_d3});
		TRUE_(y_phi() == X_phi{-5*x_d4, x_d3});
	//	TRUE_(y_phi() == X_phi{-4*x_d4, x_d3});
		TRUE_(y_phi() == X_phi{-3*x_d4, x_d3});
	//	TRUE_(y_phi() == X_phi{-2*x_d4, x_d3});
		TRUE_(y_phi() == X_phi{-1*x_d4, x_d3});

	}
	/***/
	/**/
	TRY_("procession in-place")
	{
		T_alpha x_d4 = _fit::haplo_f(4);
		T_alpha x_d3 = _fit::haplo_f(3);
		T_alpha z_outs[2][8]{};
		auto  z_out = bond::transpack_f<void_type[2]>(8, z_outs);
		using Z_out = reiterated_t<XTAL_ALL_(z_out)>;

		auto z_psi = Z_psi::bind_f();
	//	static_assert(same_q<X_phi, decltype(z_phi.store().front())>);

		//\
		occur::resize_t<> z_req(8);
		occur::resize_t<unsigned int> z_req(8);//FIX
		occur::cursor_t<> z_ren(8);
		occur::review_t<Z_out> z_rev(z_out);

		z_psi <<= occur::math::indent_s<X_phi, 1>{x_d4};
		z_psi <<= z_req;
		
		//\
		(void) z_psi.efflux(z_rev, z_ren++);
		z_psi >>= z_ren++ >> z_rev;
		//\
		TRUE_(z_out[0] == bond::pack_t<T_alpha, T_alpha>( 1*x_d4, x_d4));
		TRUE_(z_out[0] == bond::pack_f( 1*x_d4, x_d4));
		TRUE_(z_out[1] == bond::pack_f( 2*x_d4, x_d4));
		TRUE_(z_out[2] == bond::pack_f( 3*x_d4, x_d4));
		TRUE_(z_out[3] == bond::pack_f( 4*x_d4, x_d4));
		TRUE_(z_out[4] == bond::pack_f( 5*x_d4, x_d4));
		TRUE_(z_out[5] == bond::pack_f( 6*x_d4, x_d4));
		TRUE_(z_out[6] == bond::pack_f( 7*x_d4, x_d4));
		TRUE_(z_out[7] == bond::pack_f(-8*x_d4, x_d4));

		z_psi <<= occur::math::indent_s<X_phi, 1>{x_d3};

		//\
		(void) z_psi.efflux(z_rev, z_ren++);
		z_psi >>= z_ren++ >> z_rev;
		TRUE_(z_out[0] == bond::pack_f(-3*x_d3, x_d3));
		TRUE_(z_out[1] == bond::pack_f(-2*x_d3, x_d3));
		TRUE_(z_out[2] == bond::pack_f(-1*x_d3, x_d3));
		TRUE_(z_out[3] == bond::pack_f(-0*x_d3, x_d3));
		TRUE_(z_out[4] == bond::pack_f( 1*x_d3, x_d3));
		TRUE_(z_out[5] == bond::pack_f( 2*x_d3, x_d3));
		TRUE_(z_out[6] == bond::pack_f( 3*x_d3, x_d3));
		TRUE_(z_out[7] == bond::pack_f(-4*x_d3, x_d3));

	}
	/***/
	/**/
	TRY_("procession")
	{
		T_alpha x_d4 = _fit::haplo_f(4);
		T_alpha x_d3 = _fit::haplo_f(3);
		T_alpha z_outs[2][8]{};
		auto z_out = bond::transpack_f<void_type[2]>(8, z_outs);

		auto z_phi = Z_phi::bind_f();
		static_assert(same_q<X_phi, decltype(z_phi.store().front())>);

		auto constexpr z_fit = _xtd::ranges::views::transform([] XTAL_1FN_(call) (W_phi));

		occur::resize_t<> z_req(8);
		occur::cursor_t<> z_ren(8);


		z_phi <<= occur::math::indent_s<X_phi, 1>{x_d4};
		z_phi <<= z_req;
		z_phi >>= z_ren++;
		//\
		_xtd::ranges::move(z_phi|z_fit, z_out.begin());
		_xtd::ranges::copy(z_phi|z_fit, z_out.begin());

		TRUE_(z_out[0] == bond::pack_f( 1*x_d4, x_d4));
		TRUE_(z_out[1] == bond::pack_f( 2*x_d4, x_d4));
		TRUE_(z_out[2] == bond::pack_f( 3*x_d4, x_d4));
		TRUE_(z_out[3] == bond::pack_f( 4*x_d4, x_d4));
		TRUE_(z_out[4] == bond::pack_f( 5*x_d4, x_d4));
		TRUE_(z_out[5] == bond::pack_f( 6*x_d4, x_d4));
		TRUE_(z_out[6] == bond::pack_f( 7*x_d4, x_d4));
		TRUE_(z_out[7] == bond::pack_f(-8*x_d4, x_d4));
		
		z_phi <<= occur::math::indent_s<X_phi, 1>{x_d3};
		z_phi >>= z_ren++;
		//\
		_xtd::ranges::copy(z_phi|z_fit, z_out.begin());
		_xtd::ranges::copy(z_phi|z_fit, z_out.begin());
		
		TRUE_(z_out[0] == bond::pack_f(-3*x_d3, x_d3));
		TRUE_(z_out[1] == bond::pack_f(-2*x_d3, x_d3));
		TRUE_(z_out[2] == bond::pack_f(-1*x_d3, x_d3));
		TRUE_(z_out[3] == bond::pack_f(-0*x_d3, x_d3));
		TRUE_(z_out[4] == bond::pack_f( 1*x_d3, x_d3));
		TRUE_(z_out[5] == bond::pack_f( 2*x_d3, x_d3));
		TRUE_(z_out[6] == bond::pack_f( 3*x_d3, x_d3));
		TRUE_(z_out[7] == bond::pack_f(-4*x_d3, x_d3));

	}
	/***/
	/**/
	TRY_("multiplication")
	{
		T_alpha x =  0.33, x_d4 = _fit::haplo_f(4);
		T_alpha y =  5.55;
		T_alpha z =  x*y;

		Y_phi y_phi{x, x_d4};

	//	TRUE_(Y_phi{x, x_d4} == y_phi);
		TRUE_(Y_phi{x, x_d4} == y_phi.head());
		TRUE_(X_phi{x, x_d4} == y_phi.head());
		TRUE_(check_f<6>(x, y_phi.head() (0)));

		y_phi *= y;
		TRUE_(check_f<6>(y_phi.head() (0), z - _std::round(z)));

	}
	/***/
}
TAG_("phasor trials")
{
	using namespace Eigen;
	
	using U_stored  = provision::stored<unit_type[0x1000]>;
	using U_sampled = occur::resample_t<>::template attach<>;

	using _Op = bond::fit<>;
	using T_sigma = typename _Op::sigma_type;
	using T_delta = typename _Op::delta_type;
	using T_alpha = typename _Op::alpha_type;
	
	using T_cell = T_alpha;
	using T_eigencolumns = Array<T_cell, Dynamic, 2, ColMajor>;
	using T_eigenrows    = Array<T_cell, Dynamic, 2, RowMajor>;
	using T_eigenrow     = Array<T_cell,       1, 2, RowMajor>;


	using  _phi = T_cell[2];
	using W_phi = bond::repack_t<_phi>;
	using X_phi = atom::math::phason_t<_phi>;
	
	using Y_chi = process::conveyor_t<phasor<_phi, U_sampled>>;
//	using Y_chi = process::conveyor_t<phasor<_phi>>;
	using Y_phi = phasor_t<_phi>;
	//\
	using Y_psi = phasor_t<_phi, U_sampled>;
	using Y_psi = process::lift_t<bond::repack_t<_phi>, phasor<_phi>>;
	//\
	using Y_eig = process::link_t<T_eigenrow, phasor<_phi>>;
	using Y_eig = process::lift_t<T_eigenrow, phasor<_phi>>;

	using Z_chi = processor::monomer_t<Y_chi, U_stored>;
	using Z_phi = processor::monomer_t<Y_phi, U_stored>;
	using Z_psi = processor::monomer_t<Y_psi, U_stored>;
	//\
	using Z_eig = processor::monomer_t<process::lift<bond::operate<T_eigenrow>>, Y_chi>;
//	using Z_eig = processor::monomer_t<confined_t<lift<bond::operate<_std::array<T_cell, 2>>>, phasor<_phi, U_sampled>>>;

	using _fit = bond::template fit<typename X_phi::value_type>;


	/**/
	static constexpr T_alpha x_delta  = _fit::ratio_f(7);
	
	T_sigma constexpr N_data = 0x1000;
	T_alpha   z_data[2][N_data]{};
	T_alpha  *y_data   [N_data]{z_data[0], z_data[1]};
	T_alpha **x_data = y_data;
	for (int i = 0; i < N_data; ++i) {
		z_data[0][i] =  i;
		z_data[1][i] = -1;
	}


	auto w_data  = bond::transpack_f<void_type[2]>(N_data, z_data);
	//\
	auto e_data  = Map<T_eigencolumns>(*z_data, N_data, 2).rowwise();
	auto e_data  = ConvertToEigenMatrix<2>(y_data, N_data).rowwise();
	
	auto x_phi = X_phi{}; x_phi <<=                                {_fit::ratio_f(7)};
	auto y_phi = Y_phi{}; y_phi <<= occur::math::indent_s<X_phi, 1>{_fit::ratio_f(7)}; y_phi <<= occur::resize_t<>(N_data);
	
	auto z_chi = Z_chi::bind_f(); z_chi <<= occur::math::indent_s<X_phi, 1>{_fit::ratio_f(7)}; z_chi <<= occur::resize_t<>(N_data);
	auto z_phi = Z_phi::bind_f(); z_phi <<= occur::math::indent_s<X_phi, 1>{_fit::ratio_f(7)}; z_phi <<= occur::resize_t<>(N_data);
	auto z_psi = Z_psi::bind_f(); z_psi <<= occur::math::indent_s<X_phi, 1>{_fit::ratio_f(7)}; z_psi <<= occur::resize_t<>(N_data);
//	auto z_eig = Z_eig::bind_f(); z_eig <<= occur::math::indent_s<X_phi, 1>{_fit::ratio_f(7)}; z_eig <<= occur::resize_t<>(N_data);

	occur::cursor_t<>               z_cursor(N_data);
	occur::math::indent_s<X_phi, 1> z_indent{x_delta};
	
	z_phi <<= z_indent;
	z_psi <<= z_indent;
//	z_eig <<= z_indent;

	EST_("procession (process in-place)")
	{
		auto &z_d0 = z_data[0];
		auto &z_d1 = z_data[1];
		for (int i = 0; i < N_data; ++i) {
			auto const &y = y_phi();
			z_d0[i] = y(0);
			z_d1[i] = y(1);
		}

	};
	EST_("procession (processor in-place: `ranges::copy...`)")
	{
		//\
		z_psi >>= z_cursor++ >> occur::review_f(w_data);
		z_chi >>= z_cursor++ >> occur::review_f(w_data);

	};
}


///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
