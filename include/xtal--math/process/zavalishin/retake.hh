#pragma once
#include "./any.hh"






XTAL_ENV_(push)
namespace xtal::process::math::zavalishin
{/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

template <auto  ...Ms>	struct  retake;
template <auto  ...Ms>	using   retake_t = process::confined_t<retake<Ms...>>;


////////////////////////////////////////////////////////////////////////////////
/*!
\brief   Responds to `efflux(occur::stage_f(+0))` by resetting the filter state.
\todo    Allow reset for any `dispatch`ed parameters.
*/
template <>
struct retake< 0>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
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
		template <signed N_ion> requires in_v<N_ion, -1>
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
/*!
\brief   Responds to `influx(occur::stage_f(-1))` by returning `1` if the state is under threshold, `0` otherwise.
\node    Threshold is currently fixed at `2^-16`, or around `-48dBFS`.
\todo    Allow configurable threshold.
\todo    Accommodate gain when calculating the `dot` with state?
*/
template <>
struct retake<-1>
{
	template <class S>
	class subtype : public bond::compose_s<S>
	{
		using S_ = bond::compose_s<S>;

	public:// CONSTRUCT
		using S_::S_;
		using typename S_::order_attribute;
		using typename S_::state_type;

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
		template <signed N_ion> requires in_v<N_ion, +1>
		XTAL_DEF_(return,inline,let)
		fuse(occur::stage_q auto &&o)
		noexcept -> signed
		{
			auto constexpr dot_e = bond::fit<unstruct_t<state_type>>::haplo_f(0x11);// Epsilon (~-48dBFS) squared...
			auto const [states_] = S_::template memory<state_type>();
			static_assert(state_type::size() <= 4);

			signed x = S_::template fuse<N_ion>(XTAL_REF_(o));

			auto const order = order_attribute{self()};
			if (0 < order and o.head() == -1) {
				auto const disorder = order - one;
				XTAL_IF0
				XTAL_0IF (1 == state_type::size()) {
					x &= dot_e > dot_f(states_);
				}
				XTAL_0IF (2 == state_type::size()) {
					XTAL_IF1_(assume) (disorder == (disorder&0b01));
					switch                         (disorder&0b01) {
					case 0: x &= dot_e > dot_f(states_.self(constant_t<1>{})); break;
					case 1: x &= dot_e > dot_f(states_.self(constant_t<2>{})); break;
					}
				}
				XTAL_0IF (3 == state_type::size()) {
					XTAL_IF1_(assume) (disorder == (disorder&0b11));
					switch                         (disorder&0b11) {
					case 0: x &= dot_e > dot_f(states_.self(constant_t<1>{})); break;
					case 1: x &= dot_e > dot_f(states_.self(constant_t<2>{})); break;
					case 2: x &= dot_e > dot_f(states_.self(constant_t<3>{})); break;
					case 3:                                                    break;
					}
				}
				XTAL_0IF (4 == state_type::size()) {
					XTAL_IF1_(assume) (disorder == (disorder&0b11));
					switch                         (disorder&0b11) {
					case 0: x &= dot_e > dot_f(states_.self(constant_t<1>{})); break;
					case 1: x &= dot_e > dot_f(states_.self(constant_t<2>{})); break;
					case 2: x &= dot_e > dot_f(states_.self(constant_t<3>{})); break;
					case 3: x &= dot_e > dot_f(states_.self(constant_t<4>{})); break;
					}
				}
			}
			return x;
		}

	};
};


////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
}/////////////////////////////////////////////////////////////////////////////
XTAL_ENV_(pop)
