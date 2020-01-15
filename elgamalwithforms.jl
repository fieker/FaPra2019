function formdiffiehellmanA(P::NfAbsOrdIdl)
	A=P
	a=0
	count=0
	while count<(P.gen_one*100) && (Hecke.reduce_ideal(P*inv(A)).gen_one==1 || Hecke.reduce_ideal(A).gen_one==1)
		a=abs(rand(10000000:fmpz(10)^60))
		f=idealtoform(A)
		f=formpowmod(f,fmpz(a))
		A=Hecke.reduce_ideal(formtoideal(f))
		count+=1
	end
	return a,A
end

function formelgamalA(C::NfAbsOrdIdl,B::NfAbsOrdIdl,a::fmpz)
	f=idealtoform(C)
	g=idealtoform(B)
	h=nucomp(f,formpowmod(g,fmpz(-a))
	return Hecke.reduceideal(formtoideal(h))
end

function formelgamalB(M::NfAbsOrdIdl,A::NfAbsOrdIdl,P::NfAbsOrdIdl)
	f=idealtoform(M)
	g=idealtoform(A)
	h=idealtoform(P)
	b=abs(rand(10000000:fmpz(10)^60))
	i=formpowmod(h,fmpz(b))
	B=Hecke.reduce_ideal(formtoideal(i))
	j=nucomp(f,formpowmod(g,fmpz(b))))
	C=Hecke.reduce_ideal(formtoideal(j))
	return B,C
end

function testmyelgamalwithforms(m::BigInt)
	println("m=",m)
	@time p,P=diffiehellman0(m)
	#println("p=",p)
	#println("P=",P)
	@time a,A=formdiffiehellmanA(P)
	#println("A=",A)
	@time M=encode(p,m)
	#println("M=",M)
	@time B,C=formelgamalB(M,A,P)
	#println("C=",C)
	@time Mneu=formelgamalA(C,B,a)
	#println("Mneu=",Mneu)
	@time mneu=decode(Mneu)
	#println("mneu=",mneu)
	return m==mneu
end
