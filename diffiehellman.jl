######
@doc Markdown.doc"""
diffiehellman0(m::BigInt) -> BigInt, NfAbsOrdIdl
finds a random prime p>(2)^(nbits(fmpz(m))*2+100) and a nonprincipal ideal in O over Q[sqrt(-p)] 
"""
function diffiehellman0(m::BigInt)
	i=BigInt(2)^(nbits(fmpz(m))*2+100)
	p=BigInt(nextprime(rand(i:i+10^10)))
	push!(Hecke.big_primes,p)#funktioniert?
	#println("Primzahl gefunden")
	#k,a=quadratic_field(-p)
	#println("zahlkörper erstellt")
	P=findnonprincipal(BigInt(p),BigInt(floor(sqrt(p))))
	#println("ideal gefunden")
	return p,P
end

@doc Markdown.doc"""
diffiehellmanA(P::NfAbsOrdIdl) -> BigInt, NfAbsOrdIdl

determines private and public key for the receiver 
"""
function diffiehellmanA(P::NfAbsOrdIdl)
	A=P
	a=0
	while Hecke.reduce_ideal(P*inv(A)).gen_one==1 || Hecke.reduce_ideal(A).gen_one==1
		#println("nicht-trivialität geprüft")
		a=abs(rand(10000000:fmpz(10)^60))
		A=Hecke.power_class(P,fmpz(a))
		#println("Schlüssel berechnet")
	end
	#println("finaler schlüssel")
	return a,A
end

@doc Markdown.doc"""
diffiehellmanB(P::NfAbsOrdIdl) -> BigInt, NfAbsOrdIdl

determines private and public key for the transmitter
"""
function diffiehellmanB(P::NfAbsOrdIdl)
	B=P
	b=0
	while Hecke.reduce_ideal(P*inv(B)).gen_one==1 || Hecke.reduce_ideal(B).gen_one==1
		b=abs(rand(10000000:fmpz(10)^60))
		B=Hecke.power_class(P,fmpz(b))
	end
	return b,B
end

@doc Markdown.doc"""
diffiehellman1(A::NfAbsOrdIdl,b::BigInt) -> NfAbsOrdIdl

determines the key to encrypt via Diffie-Hellman
"""
function diffiehellman1(A::NfAbsOrdIdl,b::BigInt)
	K=Hecke.power_class(A,fmpz(b))
	return K
end
