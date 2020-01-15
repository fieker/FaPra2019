function formtoideal(f::QuadForm)
	a=f.a
	b=f.b
	c=f.c
	D=Hecke.discriminant(f)
	k,d=quadratic_field(D)
	max_order1(BigInt(D))
	zk=maximal_order(k)
	J=d*zk
	e=collect(keys(Hecke.factor(J)))
	w1=fmpz(a)
	w2=divexact(-b+d,2*a)
	I=NfOrdIdl(w1,w2)
	return I
end

function idealtoform(I::NfOrdIdl,p::BigInt)
	a1=I.gen_one
	M=representation_matrix(I.gen_two)
	b1=M[1]
	c1=M[3]
	d=fmpz(p)
	n=norm(I)
	a=divexact(a1^2,n)
	b=divexact(-2*a1*b1,n)
	c=divexact(b1^2-c1^2*d,n)
	f=QuadForm(a,b,c)
	return f
end		
