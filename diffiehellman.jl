function diffiehellman0(m::BigInt)
#include find.jl
	i=BigInt(2)^(nbits(fmpz(m))*2+100)
	p=BigInt(nextprime(abs(rand(Int64))))
	while p<i
		p=next_prime(10*p)
	end
	k,a=quadratic_field(-p)
	P=findnonprincipal(BigInt(p),BigInt(floor(sqrt(p))),BigInt(100000*floor(sqrt(p))))
	return p,P
end

function diffiehellmanA(P::NfAbsOrdIdl)
	A=P
	a=0
	while Hecke.reduce_ideal(P*inv(A)).gen_one==1 || Hecke.reduce_ideal(A).gen_one==1
		println(1)
		a=abs(rand(10000000:fmpz(10)^60))
		A=Hecke.power_class(P,fmpz(a))
	end
	return a,A
end

function diffiehellmanB(P::NfAbsOrdIdl)
	B=P
	b=0
	while Hecke.reduce_ideal(P*inv(B)).gen_one==1 || Hecke.reduce_ideal(B).gen_one==1
		b=abs(rand(10000000:fmpz(10)^60))
		B=Hecke.power_class(P,fmpz(b))
	end
	return b,B
end

function diffiehellman1(A::NfAbsOrdIdl,b::Int64)
	K=Hecke.power_class(A,fmpz(b))
	return K
end
