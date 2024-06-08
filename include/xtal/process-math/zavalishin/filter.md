
	co := g + a1;
	HP := (x - co*s1 - s2)/(a0 + g*co);
	HP := (x - co*s1 - s2)/(1 + g*co);
	v1 := g*HP; BP := v1+s1; s1 := BP+v1;
	v2 := g*BP; LP := v2+s2; s2 := LP+v2;

2.57089e-10	1.13376e-05	0.499989

			U_alpha const g0 = a0;
			U_alpha const g1 = term_f(a1, g, g0);
			U_alpha const g2 = term_f(a2, g, g1);
			
			U_value y2 = term_f<-1>(term_f<-1>(x, s0, g0), s1, g1)/g2;
			U_value y1{0};
			U_value y0{0};



	auto const a2 = 1;
	
	auto const g0 = a0;
	auto const g1 = term_f(a1, g, g0);
	auto const g2 = term_f(a2, g, g1);
	
	auto const y2 = term_f<-1>(term_f<-1>(x, s0, g0), s1, g1)/g2;
	
	auto const y2 = (x - s0*g0 - s1*g1)/g2;
	auto const y1 = 0;
	auto const y0 = 0;
	
	auto const [y0_, y1_, y2_] = z__filter__probe_3x1x__1x_3x1x_state_2x_1g_drive_1u_curve_1v_a01_2x__$${precision}__$${saturation} (x, y0, y1, y2, s0, s1, g, u, v, a0, a1);
	y2 = y2_;





	auto const curve_add_1 = curve + 1;
	auto const curve_sub_1 = curve - 1;
	auto const     a_sub_2 =    x1 - 2;
	for (i = 0; i < ${saturation}; i += 1) {
		y1 = term_f(v1, y2, t);
		y0 = term_f(v0, y1, t);
		auto const u_y1 = u*y1;
		auto const [exp_u_y1, inv_exp_u_y1] = f__exp_1x_inv_1x__1x__approximation_1 (u_y1, ${precision});
		
		auto const p1 = curve_add_1 *     exp_u_y1;
		auto const q1 = curve_sub_1 * inv_exp_u_y1;
		
		auto const d_add_q = p1 + q1;
		auto const d_sub_q = p1 - q1;
		auto const k1 = d_sub_q +  a_sub_2;

		auto const num = term_f<-1>(u, x0, y0) - term_f(y2, y1, k1);

		y2 += (u - x0*y0 - k1*y1 - y2)/(1 + t*(t*x0 + k1 + u_y1*d_add_q));
	}
	return y0, y1, y2;
