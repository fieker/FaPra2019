function encode(p::BigInt,m::BigInt)
#include find.jl
	M=findnonprincipal(BigInt(p),m*10^6,BigInt(m*10^7))
	return Hecke.reduce_ideal(M)
end

function elgamalB(M::NfAbsOrdIdl,A::NfAbsOrdIdl,P::NfAbsOrdIdl)
	b=abs(rand(Int64))
	B=Hecke.power_class(P,fmpz(b))	
	C=(M*(Hecke.power_class(A,fmpz(b))))
	return C,B
end

function elgamalA(C::NfAbsOrdIdl,B::NfAbsOrdIdl,a::Int64)
	M=C*(Hecke.power_class(B,fmpz(-a)))
	return Hecke.reduce_ideal(M)
end

function decode(M::NfAbsOrdIdl)
	m=divrem(M.gen_one,10^6)[1]
	return m
end

function testmyelgamal(i::BigInt,m::BigInt)
	#i=BigInt(10)
	println("i=",i)
	#m=BigInt(17)
	println("m=",m)
	p,P=diffiehellman0(i)
	println("p=",p)
	println("P=",P)
	a,A=diffiehellmanA(P)
	println("A=",A)
	M=encode(p,m)
	println("M=",M)
	C,B=elgamalB(M,A,P)
	println("C=",C)
	Mneu=elgamalA(C,B,a)
	println("Mneu=",Mneu)
	mneu=decode(Mneu)
	println("mneu=",mneu)
	if m==mneu
		print("ja")
	else
		print("nein")
	end
end
