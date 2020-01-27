function formtoideal(f::QuadForm)
	a=f.a
	b=f.b
	c=f.c
	D=fundamental(Hecke.discriminant(f))
	k,d=quadratic_field(D)
	max_order1(BigInt(abs(D)))
	zk=maximal_order(k)
	w1=fmpz(a)
	w2=zk((-b+d)//2)
	I=NfOrdIdl(w1,w2)
	return I
end

function idealtoform(I::NfOrdIdl,p::BigInt)
	A=correctlyorderedbasis(I,p)
	a=representation_matrix(A[1])[1]
	b=representation_matrix(A[2])[1]
	c=representation_matrix(A[2])[3]
	d=-p
	n=norm(I)
	f1=a^2//n
	f2=2*a*b//n
	f3=(b^2-c^2*d)//n
	f=QuadForm(f1,f2,f3)
	return f
end		

function fundamental(D::fmpz)
	if mod(D,4)==0
		m=divexact(D,4)
			if (mod(m,4)==2 || mod(m,4)==3) && issquarefree(m)
				return m
			else 
				return false
			end
	elseif issquarefree(D)
		if mod(D,4)==1
			return D
		end
	else 
		return false
	end
end

function conjugation1(x::nf_elem,p::BigInt)
	k,a=quadratic_field(-p)
	M=representation_matrix(x)
	y=M[1]-M[3]*a
	return y
end

function conjugation2(x::NfAbsOrdElem{AnticNumberField,nf_elem},p::BigInt)
	k,a=quadratic_field(-p)
	M=representation_matrix(x)
	y=M[1]-M[3]*a
	return y
end

function iscorrectlyordered(A::Array{NfAbsOrdElem{AnticNumberField,nf_elem},1},p::BigInt)
	k,a=quadratic_field(-p)
	w1=A[1]
	w2=A[2]
	return representation_matrix(((w2*conjugation2(w1,p)-w1*conjugation2(w2,p))//a))[1]>0
end

function correctlyorderedbasis(I::NfOrdIdl,p::BigInt)
	A=basis(I)
	while iscorrectlyordered(A,p)==false
		I=NfOrdIdl(representation_matrix(A[1])[1],representation_matrix(A[1])[1]+A[2])
		A=basis(I)
	end
	return A
end
