

Relations normalized `between -½ and ½`:

	{                Tan[x*π]/π, x * Sqrt[(1 - 2 x²)]/(1 - 4 x²)}
	{       Gudermannian[x*π]/π, x * Sqrt[(1 +   x²)]/(1 + 2 x²)}

	{             ArcTan[x*π]/π, x * Sqrt[2]/Sqrt[(1 + 8 x²) + Sqrt[(1 + 8 x²)]]}
	{InverseGudermannian[x*π]/π, x * Sqrt[2]/Sqrt[(1 - 4 x²) + Sqrt[(1 - 4 x²)]]}
	
	Sqrt[2]/Sqrt[(1 + 1/2) + Sqrt[(1 + 1/2)]]/4
	Sqrt[2]/Sqrt[(x² + 8) + Sqrt[(0)*1 + (1)*x^4 + 8 x²]]




	x Sqrt[4 - w x^2]/(2 - w x^2)
	

(Arc)Tan[x*π]/π:

	w = 32 - ½ π² - Sqrt[32 π² + ¼ π⁴]
	w = 8.6212390705554946929075956288725164247876795108982924297254982462...

	x Sqrt[4 - w x^2]/(2 - w x^2)
	Sqrt[2] x/Sqrt[(1 + w x²) + Sqrt[(1  + w x²)]]

So, for `Tan` we wrap at `1/4`, and if outside, need to invert...

	roots_f<2>((4 - w xx)*square_f(x/(2 - w xx)))]


	{x*root_f<2>(4 - wxx), (2 - wxx)}


For `ArcTan`, need to invert the argument and adjust by +/-0.5...

		using U_alpha = typename _op::alpha_type;
		using W_alpha = algebra::sector_t<U_alpha[2]>;
		auto u = wrap_f(o), n = _op::truncate_f<0,-2>(u);
		auto v = W_alpha(n == 0)*W_alpha{u, root_f<-1>(u)};
	//	


	(y + dn)*root_f<-2>(y)



ArcTan[x]:

	w = 2*(4/π)² - ½ - Sqrt[2*(4/π)² + ¼]
	w = 0.8735141470922506811011384624884910162505195635622737479437808099...

	Sqrt[2]*x/Sqrt[(1 + w x²) + Sqrt[(1  + w x²)]]
	Sqrt[2]*1/Sqrt[(w + 1 x²) + Sqrt[(x⁴ + w x²)]]
	
	
	Sqrt[2]*x/Sqrt[(1 + w x²) + Sqrt[(1  + w x²)]]

	x1 = dn*x + up*1;
	a0 = dn*1 + up*w;
	a1 = dn*w + up*k;
	y0 = 
	
	Sqrt[2]*x1/Sqrt[(a0 + a1 x²) + Sqrt[(y0 + w x²)]]



	{Arg[1 + I], Arg[-1 + I], Arg[-1 - I], Arg[1 - I]}
	[{1/4, 3/4, -3/4, -1/4}



	Sqrt[2]*1/Sqrt[(0.8735 + 1 x²) + Sqrt[(x⁴ + 0.8735 x²)]]



1.038101576831206261420099031130401897246868


Relations normalized `between -1 and 1`:

	{InverseGudermannian[x * (π/2)]/(π/2), x/sqrt(((1 - x²) + sqrt(1 - x²))/2)} between -1 and 1
	{ArcTan[x*(π/2)]/(π/2), (Sqrt[2] x)/sqrt(((1 + (Sqrt[2] x)²) + sqrt(1 + (Sqrt[2] x)²)))}
	{ArcTan[x], InverseGudermannian[x*Sqrt[-2]]/Sqrt[-2]} between -1 and 1
