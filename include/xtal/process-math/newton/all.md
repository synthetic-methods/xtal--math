Root-finding for odd-functions of the form `f[x] x`:

	x *= f'[x]/(f[x] x + f'[x])
	x *= 1/(1 + x*f[x]/f'[x])

Root finding for powers of `n`:

	y *= y // Module[{u = n^-1}, Module[{v = 1 - u}, Term[v, x , u #^-n]]] &
