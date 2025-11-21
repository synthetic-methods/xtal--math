XTAL Math
=========

XTAL--Math is a DSP library that utilizes template composition as the means to build higher-order components from lower-level functional and stateful processes.

At the namespace level, the `math` namespaces specialize the `xtal` folder structure, with `xtal--math/*/` mapping to `xtal::*::math::`. Within each `math` subspace, components may be further organized by an associated author or mathematician: currently Dirichlet, Fourier, Pade, Taylor, Zavalishin.

	atom/dirichlet/
	-	symbol.hh        : Generates Dirichlet characters/Kronecker symbols.

	atom/fourier/
	-	serial.hh        : Linear convolution.
	-	series.hh        : Circular convolution & transformation.

	atom/pade/
	-	wniplex.hh       : Multiplicative complex `sign`/`abs` pair/triple with mutually reciprocal magnitude.

	process/pade/
	-	unity.hh         : Complex sinusoid with normalized period `2π`, figuratively equivalent to `1^# &`.
	-	wnity.hh         : Mutually reciprocal period-normalized sinusoids.

	process/taylor/
	-	logarithm.hh     : Various approximations of the     logarithm  `Log`.
	-	monologarithm.hh : Various approximations of the monologarithm `-Log[1 - #] &`.

	process/zavalishin/
	-	curve.hh         : (Anti)saturation curves.
	-	differ.hh        : First-order differences with amplitude normalization.
	-	filter.hh        : Arbitrary-order State Variable Filter (SVF).
	-	gate.hh          : Gate control for `filter`.
	-	vactrol.hh       : Envelope generation using `filter`.

(NOTE: Like [`xtal`](///github.com/synthetic-methods/xtal), each `*.hh` is accompanied by a test-suite `*.cc` which is not included by `CMake`/`Conan`.)

At the class level, components are parameterized both by the containing `struct` and the enclosed `method`s.
The intermediate template `subtype` provides the mechanism for decorator composition via `xtal::bond::compose`.

	template <auto ...Ms>
	struct foo
	{
		template <class S>
		class subtype
		{
			template <auto ...Ns>
			auto method(auto ...xs) -> noexcept decltype(auto)
			{
			//	...	
			}

		};
	};

The outer `Ms...` of the `struct` are used to specify compile-time configuration.
The inner `Ns...` of the `method` are enumerated and made available via _vtable_.

	template <int M_ism=1>
	struct sine
	{
		template <class S>
		class subtype
		{
			template <int N_ord=1>
			auto method(auto ...xs) -> noexcept decltype(auto)
			{
			//	Evaluate polynomial of order `N_ord`...
			//	...when `M_ism ==  2` approximate the hyperbolic function.
			//	...when `M_ism ==  1` approximate the   circular function.
			//	...when `M_ism == -1` approximate the   circular arc.
			//	...when `M_ism == -2` approximate the hyperbolic area.
			}

		};
	};

Instances of the `struct` may thus be composed with other components.

	using haversine = bond::compose<square<1>, sine<1>, halve<1>>;

For more information on the composition mechanism, see the [`xtal`](///github.com/synthetic-methods/xtal) documentation.
Detail on individual components can be found in the comments or associated Doxygen output.
The most comprehensive reference for usage is currently the test suite (sorry!).
