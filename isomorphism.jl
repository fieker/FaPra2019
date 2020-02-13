function sign1(x::fmpq)
	return sign(numerator(x))
end

function sign1(x::fmpz)
	return sign(x)
end

function coeffs(nf_elem)
        return coeff(x,0),coeff(x,1)
end

function coeffs(x::NfAbsOrdElem{AnticNumberField,nf_elem})
        return coeff(nf(parent(x))(x),0),coeff(nf(parent(x))(x),1)
end

function conjugation(x::NfAbsOrdElem{AnticNumberField,nf_elem})
	c = coeffs(x)
    return parent(x)([c[1], -c[2]])
end

function conjugation(x::nf_elem)
	c = coeffs(x)
    return parent(x)([c[1], -c[2]])
end

@doc Markdown.doc"""
iscorrectly_ordered(A::Array{NfAbsOrdElem{AnticNumberField,nf_elem},1},p::BigInt) -> Bool
A=basis of an ideal over in O over Q[sqrt(-p)]
w1=A[1]
w2=A[2]
(norm(w2*conj(w1)-w1*conj(w2))/sqrt(-p))>0??
"""
function iscorrectly_ordered(A::Array{NfAbsOrdElem{AnticNumberField,nf_elem},1},p::BigInt)
	w1=A[1]
	w2=A[2]
	a=coeffs(w1)[1]
	c=coeffs(w2)[1]
	return sign1(a)==sign1(c)
end

@doc Markdown.doc"""
iscorrectly_ordered(A::Array{fmpz,2}) -> Bool
A=basis of an ideal over in O over Q[sqrt(-p)]
w1=A[1]
w2=A[2]
(norm(w2*conj(w1)-w1*conj(w2))/sqrt(-p))>0??
"""
function iscorrectly_ordered(A::Array{fmpz,2})
	a=A[1]
	c=A[4]
	return sign1(a)==sign1(c)
end
	
@doc Markdown.doc"""
correctlyorderedbasis(I::NfOrdIdl,p::BigInt) -> NfAbsOrdElem{AnticNumberField,nf_elem},1}
returns a correctly ordered basis of the ideal I
"""
function correctlyorderedbasis(I::NfOrdIdl,p::BigInt)
	A=basis(I)
	if iscorrectly_ordered(A,p)[1]==false
		A[2]=-1*A[2]
	end
	return A
end

@doc Markdown.doc"""
correctlyorderedbasis(A::Array{fmpz,2}) -> NfAbsOrdElem{AnticNumberField,nf_elem},1}
returns A in correctly order
"""	
function correctlyorderedbasis(A::Array{fmpz,2})
	if iscorrectly_ordered(A)[1]==false
		A[2]=-1*A[2]
		A[4]=-1*A[4]
	end
	return A
end

@doc Markdown.doc"""
fundamental(D::fmpz) -> 
returns the corresponding discriminant of a quadratic number field if D is a fundamental discriminant of a binary quadratic form
"""	
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

@doc Markdown.doc"""
formtoideal(f::QuadForm)-> NfOrdIdl
determines the corresponding ideal in O over Q[sqrt(-p)] where p is the fundamental discriminant of f
"""
function formtoideal(f::QuadForm)
	a=f.a
	b=f.b
	c=f.c
	D=Hecke.discriminant(f)
	k,d=quadratic_field(fundamental(D))
	max_order1(BigInt(abs(fundamental(D))))
	zk=maximal_order(k)
	w1=fmpz(a)
	if mod(D,4)==0
		w2=zk((-b+2*d)//2)
	else
		w2=zk((-b+d)//2)
	end
	I=NfOrdIdl(w1,w2)
	return I
end

@doc Markdown.doc"""
idealtoform(I::NfOrdIdl,p::BigInt,s::Int64) -> QuadForm
determines the corresponding binary quadratic form with discriminant p or 4*p
"""
function idealtoform(I::NfOrdIdl,p::BigInt,s::Int64)
	A=correctlyorderedbasis(I,p)
	a=coeffs(A[1])[1]
	b=coeffs(A[2])[1]
	c=coeffs(A[2])[2]
	d=-p
	n=norm(I)
	f1=FlintZZ(divexact(a^2,n))
	f2=FlintZZ(divexact(-s*2*a*b,n))
	f3=FlintZZ(divexact((b^2-c^2*d),n))
	f=QuadForm(s*f1,s*f2,s*f3)
	return f
end

@doc Markdown.doc"""
formtobasis(f::QuadForm) -> Array{fmpz,2}
determines a basis of the corresponding ideal in O over Q[sqrt(-p)] where p is the fundamental discriminant
"""
function formtobasis(f::QuadForm)
	a=f.a
	b=f.b
	c=f.c
	if mod(Hecke.discriminant(f),4)==0
		A=[a 0;fmpz(numerator(-b//2)) 1]
	else
		A=[a 0;-b//2 1//2]
	end
	return A
end

@doc Markdown.doc"""
basistoform(A::Array{fmpz,2},d::fmpz)
determines the corresponding binary quadratic form with discriminant d or 4*d
"""
function basistoform(A::Array{fmpz,2},d::fmpz)
	A=correctlyorderedbasis(A)
	p=A[1]
	a=-A[2]
	n=gcd(gcd(p^2,a^2-d),-2*a*p)#richtig?
	f1=numerator(divexact(p^2,n))
	f2=numerator(divexact(2*a*p,n))
	f3=numerator(divexact((a^2-d),n))
	f=QuadForm(f1,f2,f3)
end
