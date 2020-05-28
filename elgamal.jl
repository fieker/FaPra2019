#####
@doc Markdown.doc"""
encode(p::BigInt,m::BigInt) -> NfAbsOrdIdl

encodes Integer m to an ideal in O over Q[sqrt(-p)]
"""
function encode(p::BigInt,m::BigInt)
	M=findnonprincipal(BigInt(p),m*10^8)
	#println("Nachricht kodiert")
	return Hecke.reduce_ideal(M)
end

@doc Markdown.doc"""
elgamalB(M::NfAbsOrdIdl,A::NfAbsOrdIdl,P::NfAbsOrdIdl) -> NfAbsOrdIdl, NfAbsOrdIdl

determines private key of the transmitter via P
|
encryts M with the public key A of the receiver 
"""
function elgamalB(M::NfAbsOrdIdl,A::NfAbsOrdIdl,P::NfAbsOrdIdl)
	b=abs(rand(10000000:fmpz(10)^60))
	B=Hecke.power_class(P,fmpz(b))	
	C=Hecke.reduce_ideal(M*(Hecke.power_class(A,fmpz(b))))
	#println("Nachricht verschlüsselt")
	return B,C
end

@doc Markdown.doc"""
elgamalA(C::NfAbsOrdIdl,B::NfAbsOrdIdl,a::fmpz) -> NfAbsOrdIdl

decrypts C with the public key B of the transmitter and the private key a of the receiver
"""
function elgamalA(C::NfAbsOrdIdl,B::NfAbsOrdIdl,a::fmpz)
	M=C*(Hecke.power_class(B,fmpz(-a)))
	#println("Nachricht entschlüsselt")
	return Hecke.reduce_ideal(M)
end

@doc Markdown.doc"""
decode(M::NfAbsOrdIdl) -> BigInt

decodes the decrypted ideal
"""
function decode(M::NfAbsOrdIdl)
	m=divrem(M.gen_one,10^8)[1]
	#println("Nachricht dekodiert")
	return BigInt(m)
end

function testmyelgamal(m::BigInt)
	println("m=",m)
	@time p,P=diffiehellman0(m)
	#println("p=",p)
	#println("P=",P)
	@time M=encode(p,m)
	#println("M=",M)
	@time a,A=diffiehellmanA(P)
	#println("A=",A)
	@time B,C=elgamalB(M,A,P)
	#println("C=",C)
	@time Mneu=elgamalA(C,B,a)
	#println("Mneu=",Mneu)
	@time mneu=decode(Mneu)
	#println("mneu=",mneu)
	return m==mneu
end

function timeofideals(m::BigInt)
	@time testmyelgamal(m)
end
