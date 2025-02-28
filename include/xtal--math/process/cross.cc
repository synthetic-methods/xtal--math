#pragma once
#include "./any.cc"





#include "./cross.hh"// testing...
XTAL_ENV_(push)
namespace xtal::process::math::_test
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

struct mix
{
	class type : public process::confine_t<type>
	{
	public:

		XTAL_DEF_(return,inline,let)
		method(auto &&...xs) const
		noexcept -> auto
		{
			return (XTAL_REF_(xs) +...+ 0);
		}

	};
};
using Px_mix = typename mix::type;


////////////////////////////////////////////////////////////////////////////////

TAG_("cross", "process")
{
	TRY_("shape with matrix")
	{
		using namespace _xtd::ranges::views;

		using U_matrix = atom::quanta_t<int[3][2]>;
		using U_remix = cross_t<U_matrix, Px_mix>;

		auto io = U_remix();
		io <<= U_matrix {{1, 2}, {3, 4}, {5, 6}};
		TRUE_(io(10, 100) == 1290);
	//	(1*10 + 2*100) + (3*10 + 4*100) + (5*10 + 6*100)

	}
}


////////////////////////////////////////////////////////////////////////////////

TAG_("cross", "processor")
{
	/**/
	TRY_("shape with matrix")
	{
		using namespace _xtd::ranges::views;

		using U_matrix = atom::quanta_t<int[3][2]>;// 3-outputs, 2-inputs
		using U_remix = cross_t<U_matrix, Px_mix>;
		using U_mixer = processor::monomer_t<U_remix, provision::stored<>>;

		//\
		auto _1 = processor::let_f(1);
		auto _1 = processor::let_f(_xtd::ranges::views::repeat(1));
		auto _n = processor::let_f(iota(0, 10));
	//
		auto io = U_mixer::bind_f(_1, _n);
		io <<= U_matrix {{1, 2}, {3, 4}, {5, 6}};
		io <<= occur::resize_t<>(3);
		io >>= occur::cursor_t<>(3);

		TRUE_(equal_f(io, _std::vector { 9, 21, 33}));
	//	(1*1 + 2*0) + (3*1 + 4*0) + (5*1 + 6*0)
	//	(1*1 + 2*1) + (3*1 + 4*1) + (5*1 + 6*1)
	//	(1*1 + 2*2) + (3*1 + 4*2) + (5*1 + 6*2)

	//	{{1,2}, {3,4}, {5,6}}.{1,0} // Total
	//	{{1,2}, {3,4}, {5,6}}.{1,1} // Total
	//	{{1,2}, {3,4}, {5,6}}.{1,2} // Total

	}
	/***/
	/**/
	TRY_("shape with matrix indent")
	{
		using namespace _xtd::ranges::views;
		//\
		using U_vector = atom::quanta_t<int[2]>;
		using U_vector = bond::pack_t<int, int>;
		using U_matrix = atom::quanta_t<U_vector[3]>;

		using U_remix = cross_t<U_matrix, Px_mix>;
		using U_mixer = processor::monomer_t<U_remix, provision::stored<>>;

		auto _1 = processor::let_f(1);
		auto _n = processor::let_f(iota(0, 10));
	//
		auto io = U_mixer::bind_f(_1, _n);
		io <<= occur::math::indent_s<U_matrix>({{1, 2}, {3, 4}, {5, 6}});
		io <<= occur::resize_t<>(3);
		io >>= occur::cursor_t<>(3);

		TRUE_(equal_f(io, _std::vector { 9, 21, 33}));

	}
	/***/
	/**/
	TRY_("shape with column indent")
	{
		using namespace _xtd::ranges::views;
		//\
		using U_matrix = atom::quanta_t<int[3][2]>;
		using U_matrix = atom::quanta_t<bond::pack_t<int, int>[3]>;

		using U_remix = cross_t<U_matrix, Px_mix>;
		using U_mixer = processor::monomer_t<U_remix, provision::stored<>>;


		auto _1 = processor::let_f(1);
		auto _n = processor::let_f(iota(0, 10));
	//
		auto io = U_mixer::bind_f(_1, _n);
		io <<= occur::math::indent_s<U_matrix, 0>({1, 2});
		io <<= occur::math::indent_s<U_matrix, 1>({3, 4});
		io <<= occur::math::indent_s<U_matrix, 2>({5, 6});
		io <<= occur::resize_t<>(3);
		io >>= occur::cursor_t<>(3);

		TRUE_(equal_f(io, _std::vector { 9, 21, 33}));

	}
	/***/
	/**/
	TRY_("shape with vector indent")
	{
		using namespace _xtd::ranges::views;
		//\
		using U_matrix = atom::quanta_t<int[3][2]>;
		using U_matrix = atom::quanta_t<bond::pack_t<int, int>[3]>;

		using U_remix = cross_t<U_matrix, Px_mix>;
		using U_mixer = processor::monomer_t<U_remix, provision::stored<>>;

		auto _1 = processor::let_f(1);
		auto _n = processor::let_f(iota(0, 10));
	//
		auto io = U_mixer::bind_f(_1, _n);
		io <<= occur::math::indent_s<U_matrix, 0, 0>(1);
		io <<= occur::math::indent_s<U_matrix, 0, 1>(2);
		io <<= occur::math::indent_s<U_matrix, 1, 0>(3);
		io <<= occur::math::indent_s<U_matrix, 1, 1>(4);
		io <<= occur::math::indent_s<U_matrix, 2, 0>(5);
		io <<= occur::math::indent_s<U_matrix, 2, 1>(6);
		io <<= occur::resize_t<>(3);
		io >>= occur::cursor_t<>(3);

		TRUE_(equal_f(io, _std::vector { 9, 21, 33}));

	}
	/***/
}

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
