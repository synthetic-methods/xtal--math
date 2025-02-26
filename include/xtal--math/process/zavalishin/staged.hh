#pragma once
#include "./any.hh"

#include "./filter.hh"




XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <auto  ...Ms>	struct  staged;
template <auto  ...Ms>	using   staged_t = process::confined_t<staged<Ms...>>;


////////////////////////////////////////////////////////////////////////////////
///\
Responds to `efflux(occur::stage_f(+0))` by resetting the filter state. \

///\todo\
Allow reset for any `dispatch`ed parameters. \

template <>
struct staged< 0>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;

	public:// ACCESS
		using S_::self;

	public:// FLOW

		template <signed N_ion>
		XTAL_DEF_(return,inline,let)
		fuse(auto &&o)
		noexcept -> signed
		{
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}
		template <signed N_ion> requires in_n<N_ion, -1>
		XTAL_DEF_(return,inline,let)
		fuse(occur::stage_q auto &&o)
		noexcept -> signed
		{
			if (o.head() == 0) {
				S_::memory(constant_t<>{});
			}
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}

	};
};


////////////////////////////////////////////////////////////////////////////////
///\
Responds to `influx(occur::stage_f(-1))` by returning `1` if the state is under threshold, \
`0` otherwise. \

///\todo\
Allow configurable threshold. \

template <>
struct staged<-1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		static_assert(filter_q<S>);
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::order_type;
		using typename S_::state_type;
		using typename S_::input_type;

	public:// ACCESS
		using S_::self;

	public:// FLOW

		template <signed N_ion>
		XTAL_DEF_(return,inline,let)
		fuse(auto &&o)
		noexcept -> signed
		{
			return S_::template fuse<N_ion>(XTAL_REF_(o));
		}
		template <signed N_ion>
		XTAL_DEF_(return,inline,let)
		fuse(flow::assessing_q<occur::stage_t<>> auto &&o)
		noexcept -> signed
		{
			return fuse<+1>(XTAL_REF_(o).tail());
		}
		template <signed N_ion> requires in_n<N_ion, +1>
		XTAL_DEF_(return,inline,let)
		fuse(occur::stage_q auto &&o)
		noexcept -> signed
		{
			auto const [states_] = S_::template memory<state_type>();
			signed x = S_::template fuse<N_ion>(XTAL_REF_(o));

			if (o.head() == -1) {
			//	TODO: Accommodate gain when calculating the product?
			//	TODO: Address clumsy dynamic-sized dot-product.

				auto const order = order_type{self()};
				input_type y{};
				for (unsigned int i{}; i < order; ++i) {
					y = term_f<1, 2>(y, states_[i]);
				}
				x &= dot_f(y) < bond::fit<input_type>::haplo_f(7 + 1);//NOTE: Squared...
			}
			return x;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
