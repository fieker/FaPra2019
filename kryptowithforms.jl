######
@doc Markdown.doc"""
find(d::fmpz) -> QuadForm

returns a nontrivial reduced binary quadratic form with discriminant -d
the corresponding ideal is nonprincipal
"""
function find(d::fmpz,s::fmpz)
	i=nextprime(BigInt(s))
	while jacobi(mod(BigInt(d),i),fmpz(i))!=1 && fmpz(i)<div(root(abs(d),2),4)
		i=nextprime(i+1)
	end
	i=fmpz(i)
	I=[i 0;-lift(root(GF(i)(d),2)) 1]
	return basistoform(I,d)
end

@doc Markdown.doc"""
diffiehellman0(m::fmpz) -> fmpz
finds a random prime p>(2)^(nbits(fmpz(m))*2+100) with mod(-d,4)=3 
"""
function diffiehellman0(m::fmpz)
	i=BigInt(2)^(nbits(m)*2+100)
	d=fmpz(nextprime(rand(i:i+10^10)))
	while mod(-d,4)!=3
		d=fmpz(nextprime(rand(i:i+10^10)))
	end
	push!(Hecke.big_primes,d)
	return -d,find(-d,div(root(d,2),4*fmpz(2)^(div(nbits(d),2))))
end

@doc Markdown.doc"""
encode(m::fmpz,d::fmpz) -> QuadForm

encodes Integer m to an binary quadratic form with discriminant -d
"""
function encode(m::fmpz,d::fmpz)
	M=find(d,m*fmpz(10)^10)
	#for i=m*fmpz(10)^3:fmpz(floor(BigFloat(root(abs(d),2))/4))
	#	if Hecke.isprime(i) && jacobi(mod(d,i),i)==1
	#		M=[i 0;-lift(root(GF(i)(d),2)) 1]
	#		return basistoform(M,d)
	#	end
	#end
end

@doc Markdown.doc"""
decode(M::QuadForm) -> fmpz

decodes the decrypted binary quadratic form
"""
function decode(M::QuadForm)
	M=formtobasis(M)[1]
	m=div(M,fmpz(10)^10)
	return m
end

@doc Markdown.doc"""
diffiehellmanA(I::QuadForm,d::fmpz) -> fmpz, QuadForm

determines private and public key for the receiver 
"""
function diffiehellmanA(I::QuadForm,d::fmpz)
	a=rand(1000:fmpz(10)^6)
	A=formtobasis(formpowmodbit(I,fmpz(a)))
	while isone(A)
		a=rand(1000:fmpz(10)^6)
		A=formtobasis(formpowmodbit(I,fmpz(a)))
	end
	return a,basistoform(A,d)
end

@doc Markdown.doc"""
diffiehellmanB(I::QuadForm,d::fmpz) -> fmpz, QuadForm

determines private and public key for the transmitter
"""
function diffiehellmanB(I::QuadForm,d::fmpz)
	b=rand(1000:fmpz(10)^6)
	B=formtobasis(formpowmodbit(I,fmpz(b)))
	while isone(B)
		b=rand(1000:fmpz(10)^6)
		B=formtobasis(formpowmodbit(I,fmpz(b)))
	end
	return b,basistoform(B,d)
end

@doc Markdown.doc"""
elgamalB(b::fmpz,A::QuadForm,M::QuadForm) -> QuadForm

determines private key of the transmitter via I
|
encryts M with the public key A of the receiver 
"""
function elgamalB(b::fmpz,A::QuadForm,M::QuadForm)
	K=formpowmodbit(A,b)
	C=nucomp(K,M)
	return C
end

@doc Markdown.doc"""
elgamalA(a::fmpz,B::QuadForm,C::QuadForm) -> QuadForm

decrypts C with the public key B of the transmitter and the private key a of the receiver
"""
function elgamalA(a::fmpz,B::QuadForm,C::QuadForm)
	L=formpowmodbit(B,-a)
	M=nucomp(L,C)
	return M
end

function testmyelgamalwithforms(m::fmpz)
	d,I=diffiehellman0(m)
	#println("d=",d)
	#println("I=",I)
	#println(1)
	@time M=encode(m,d)
	#println("M=",M)
	#println(2)
	@time a,A=diffiehellmanA(I,d)
	#println("a=",a)
	#println("A=",A)
	#println(3)
	@time b,B=diffiehellmanB(I,d)
	#println("b=",b)
	#println("B=",B)
	#println(4)
	@time C=elgamalB(b,A,M)
	#println("C=",C)
	#println(5)
	@time Mneu=elgamalA(a,B,C)
	#println("Mneu=",Mneu)
	#println(6)
	@time mneu=decode(Mneu)
	#println("mneu=",mneu)
	#println(7)
	return mneu==m
end

function timeofforms(m::fmpz)
	@time testmyelgamalwithforms(m)
end 

function smalltest()
	d=fmpz(-17)
	println("d=",d)
	I=fmpz[3 0;-lift(root(GF(3)(d),2)) 1]
	I=basistoform(I,d)
	println("I=",I)
	I=positive_definite_form(formreduce(positive_definite_form(I)))
	println("I=",I)
	M=fmpz[107 0;-lift(root(GF(107)(d),2)) 1]
	M=basistoform(M,d)
	println("M=",M)
	M=positive_definite_form(formreduce(positive_definite_form(M)))
	println("M=",M)
	a=fmpz(7)
	A=formpowmod(I,a)
	println("A=",A)
	A=positive_definite_form(formreduce(positive_definite_form(A)))
	println("A=",A)
	b=fmpz(9)
	B=formpowmod(I,b)
	println("B=",B)
	B=positive_definite_form(formreduce(positive_definite_form(B)))
	println("B=",B)
	C=elgamalB(b,A,M)
	println("C=",C)
	C=positive_definite_form(formreduce(positive_definite_form(C)))
	println("C=",C)
	Mneu=elgamalA(a,B,C)
	println("Mneu=",Mneu)
	Mneu=positive_definite_form(formreduce(positive_definite_form(Mneu)))
	println("Mneu=",Mneu)
	return isequal(Mneu,M)
end

function smalltest2()
	d=fmpz(-6101)
	println("d=",d)
	m=fmpz(3)
	println("m=",m)
	I=fmpz[5 0;-lift(root(GF(5)(d),2)) 1]
	I=basistoform(I,d)
	println("I=",I)
	I=formreduce(I)
	println("I=",I)
	M=fmpz[37 0;-lift(root(GF(37)(d),2)) 1]
	M=basistoform(M,d)
	println("M=",M)
	M=formreduce(M)
	println("M=",M)
	a=fmpz(7)
	A=formpowmod(I,a)
	println("A=",A)
	A=formreduce(A)
	println("A=",A)
	b=fmpz(9)
	B=formpowmod(I,b)
	println("B=",B)
	B=formreduce(B)
	println("B=",B)
	C=elgamalB(b,A,M)
	println("C=",C)
	C=formreduce(C)
	println("C=",C)
	Mneu=elgamalA(a,B,C)
	println("Mneu=",Mneu)
	Mneu=formreduce(Mneu)
	println("Mneu=",Mneu)
	mneu=div(formtobasis(Mneu)[1],10)
	println("mneu=",mneu)
	return m==mneu
end
