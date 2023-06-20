function Q = CHtoRS(Q)
	Q = cat(1, cat(2, real(Q), -imag(Q)), cat(2, imag(Q), real(Q)));
