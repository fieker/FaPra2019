function coefficientsz(x::NfAbsOrdElem{AnticNumberField,nf_elem},p::BigInt)
	k,d=quadratic_field(-p)
	a=trace(x)//2
	b=trace(((2*x-trace(x))*d//(2*d^2)))//2#
	return a,b
end

function coefficientsq(x::nf_elem,p::BigInt)
	k,d=quadratic_field(-p)
	a=divexact(trace(x),2)
	b=divexact(x-a,d)
	return a,b
end

function conjugationz(x::NfAbsOrdElem{AnticNumberField,nf_elem},p::BigInt)
	k,d=quadratic_field(-p)
	max_order1(p)
	zk=maximal_order(k)
	a,b=coefficientsz(x,p)
	y=a-b*d
	return y
end

function conjugationq(x::nf_elem,p::BigInt)
	k,d=quadratic_field(-p)
	max_order1(p)
	zk=maximal_order(k)
	a,b=coefficientsq(x,p)
	y=zk(a-b*d)
	return y
end

function iscorrectlyordered(A::Array{NfAbsOrdElem{AnticNumberField,nf_elem},1},p::BigInt)
	k,a=quadratic_field(-p)
	w1=A[1]
	w2=A[2]
	return (coefficientsq((w2*conjugationz(w1,p)-w1*conjugationz(w2,p))//a,p)[1])>0
end

function iscorrectlyordered2(A::Array{NfAbsOrdElem{AnticNumberField,nf_elem},1},p::BigInt)
	k,a=quadratic_field(-p)
	w1=A[1]
	w2=A[2]
	a=coefficientsz(w1,p)[1]
	c=coefficientsz(w2,p)[1]
	return (a*c)>0
end
	
function correctlyorderedbasis(I::NfOrdIdl,p::BigInt)
	A=basis(I)
	while iscorrectlyordered(A,p)==false
		A[2]+=A[1]
	end
	return A
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

function idealtoform(I::NfOrdIdl,p::BigInt,s::Int64)
	A=correctlyorderedbasis(I,p)
	a=coefficientsz(A[1],p)[1]
	b=coefficientsz(A[2],p)[1]
	c=coefficientsz(A[2],p)[2]
	d=-p
	#println("a=",a)
	#println("b=",b)
	#println("c=",c)
	n=norm(I)
	f1=fmpz(BigFloat(divexact(a^2,n)))#
	f2=fmpz(BigFloat(divexact(-s*2*a*b,n)))#
	f3=fmpz(BigFloat(divexact((b^2-c^2*d),n)))#
	#println("f1=",f1)
	#println("f2=",f2)
	#println("f3=",f3)
	#println("tf1=",typeof(f1))
	#println("tf2=",typeof(f2))
	#println("tf3=",typeof(f3))
	f=QuadForm(s*f1,s*f2,s*f3)
	return f
end
