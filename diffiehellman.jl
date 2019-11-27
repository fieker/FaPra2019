function diffiehellman0(m::BigInt)
#include find.jl
	i=BigInt(2)^(nbits(fmpz(m))*2+100)
	p=BigInt(nextprime(abs(rand(Int64))))
	while p<i
		p=next_prime(10*p)
	end
	k,a=quadratic_field(-p)
	P=findnonprincipal(BigInt(p),BigInt(2),BigInt(100000))
	return p,P
end

function diffiehellmanA(P::NfAbsOrdIdl)
	a=abs(rand(10000000:fmpz(10)^60))
	A=Hecke.power_class(P,fmpz(a))
	return a,A
end

function diffiehellmanB(P::NfAbsOrdIdl)
	b=abs(rand(10000000:fmpz(10)^60))
	B=Hecke.power_class(P,fmpz(b))
	return b,B
end

function diffiehellman1(A::NfAbsOrdIdl,b::Int64)
	K=Hecke.power_class(A,fmpz(b))
	return K
end
