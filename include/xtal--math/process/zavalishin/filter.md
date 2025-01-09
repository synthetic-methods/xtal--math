Fixed point:

	y2 => x - f0[y0] - f1[y1]
	y2 => x - f0[s0 + g (s1 + g y2)] - f1[s1 + g y2]

x - a0 y0 - y1 (a1$ + 2 F$1540 y1)

Newton's method:

	y[2] -= f(y[2])/f'(y[2])


	y1 = s1 + g y2;
	y0 = s0 + g y1;
	num = x - y2 - f0[y0] - f1[y1];
	nom = 1 + t(f1'(y1) + t*f0'(y0))
	y[2] += num/nom

	f[0]  = a[0]*#
	f[0]' = a[0]

	num = x - y2 - a0*y0 - f1[y1];
	nom = 1 + t(f1'(y1) + t*a0)


	f[1]  = 2 (F (#) + (½ a[1] - 1))*#
	f[1]' = 2 (F'(#) + (½ a[1] - 1))

 

(-x + a0 y0 + y1 (a1$ + 2 F$1540 y1) + y2)/(1 + 
 g (a1$ + a0 g + 4 F$1540 y1))

	co := g + a1;
	HP := (x - co*s1 - s2)/(a0 + g*co);
	HP := (x - co*s1 - s2)/(1 + g*co);
	v1 := g*HP; BP := v1+s1; s1 := BP+v1;
	v2 := g*BP; LP := v2+s2; s2 := LP+v2;

2.57089e-10	1.13376e-05	0.499989


