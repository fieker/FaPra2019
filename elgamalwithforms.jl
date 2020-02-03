function formdiffiehellmanA(P::NfAbsOrdIdl,p::BigInt)
	A=P
	a=0
	count=0
	while count<(P.gen_one*100) && (Hecke.reduce_ideal(P*inv(A)).gen_one==1 || Hecke.reduce_ideal(A).gen_one==1)
		a=abs(rand(10000000:fmpz(10)^60))
		f=idealtoform(A,p,1)
		f=formpowmod(f,fmpz(a))
		A=Hecke.reduce_ideal(formtoideal(f))
		count+=1
	end
	return a,A
end

function formelgamalA(C::NfAbsOrdIdl,B::NfAbsOrdIdl,a::fmpz,p::BigInt)
	f=idealtoform(C,p,1)
	g=idealtoform(B,p,1)
	h=nucomp(f,formpowmod(g,fmpz(-a)))
	return Hecke.reduceideal(formtoideal(h))
end

function formelgamalB(M::NfAbsOrdIdl,A::NfAbsOrdIdl,P::NfAbsOrdIdl,p::BigInt)
	f=idealtoform(M,p,1)
	g=idealtoform(A,p,1)
	h=idealtoform(P,p,1)
	b=abs(rand(10000000:fmpz(10)^60))
	i=formpowmod(h,fmpz(b))
	B=Hecke.reduce_ideal(formtoideal(i))
	j=nucomp(f,formpowmod(g,fmpz(b)))
	C=Hecke.reduce_ideal(formtoideal(j))
	return B,C
end

function testmyelgamalwithforms(m::BigInt)
	println("m=",m)
	#@time 
	p,P=diffiehellman0(m)
	println("p=",p)
	println("P=",P)
	#@time 
	a,A=formdiffiehellmanA(P,p)
	println("A=",A)
	#@time 
	M=encode(p,m)
	println("M=",M)
	#@time 
	B,C=formelgamalB(M,A,P,p)
	println("C=",C)
	#@time 
	Mneu=formelgamalA(C,B,a,p)
	println("Mneu=",Mneu)
	#@time 
	mneu=decode(Mneu)
	println("mneu=",mneu)
	return m==mneu
end
