function formtoideal(f::QuadForm)
	a=f.a
	b=f.b
	c=f.c
	D=Hecke.discriminant(f)
	I=NfOrdIdl(fmpz(a),NfAbsOrdElem{AnticNumberField,nf_elem}(divexact(-b+sqrt(D),2*a)))
	#zweiter Erzeuger!!
	return I
end

function idealtoform(I::NfOrdIdl)
	a1=I.gen_one
	M=representation_matrix(I.gen_two)
	b1=M[1]
	c1=M[3]
	d=M[2]
	n=norm(I)
	a=divexact(a1^2,n)
	b=divexact(-2*a1*b1,n)
	c=divexact(b1^2-c1^2*d,n)
	f=QuadForm(a,b,c)
	return f
end		
