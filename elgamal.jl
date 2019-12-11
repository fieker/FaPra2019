function encode(p::BigInt,m::BigInt)
#include find.jl
	M=findnonprincipal(BigInt(p),m*10^8,m*10^9)
	#println("Nachricht kodiert")
	return Hecke.reduce_ideal(M)
end

function elgamalB(M::NfAbsOrdIdl,A::NfAbsOrdIdl,P::NfAbsOrdIdl)
	b=abs(rand(10000000:fmpz(10)^60))
	B=Hecke.power_class(P,fmpz(b))	
	C=Hecke.reduce_ideal(M*(Hecke.power_class(A,fmpz(b))))
	#println("Nachricht verschlüsselt")
	return C,B
end

function elgamalA(C::NfAbsOrdIdl,B::NfAbsOrdIdl,a::fmpz)
	M=C*(Hecke.power_class(B,fmpz(-a)))
	#println("Nachricht entschlüsselt")
	
	return Hecke.reduce_ideal(M)
end

function decode(M::NfAbsOrdIdl)
	m=divrem(M.gen_one,10^8)[1]
	#println("Nachricht dekodiert")
	return m
end

function testmyelgamal(m::BigInt)
	println("m=",m)
	@time p,P=diffiehellman0(m)
	#println("p=",p)
	#println("P=",P)
	@time a,A=diffiehellmanA(P)
	#println("A=",A)
	@time M=encode(p,m)
	#println("M=",M)
	@time C,B=elgamalB(M,A,P)
	#println("C=",C)
	@time Mneu=elgamalA(C,B,a)
	#println("Mneu=",Mneu)
	@time mneu=decode(Mneu)
	#println("mneu=",mneu)
	return m==mneu
end
