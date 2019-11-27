function findnonprincipal(p::BigInt,m::BigInt,n::BigInt)
#prime_decomposition, using Primes
	k,a=quadratic_field(-p)
	zk=maximal_order(k)
	i=nextprime(m)
	while i<n
		x=prime_decomposition(zk,fmpz(i))
		for j=1:length(x)
			if isone(Hecke.reduce_ideal(x[j][1]))==false
				return x[j][1]
			end
		end
		i=nextprime(i+1)
	end
	print("alle Ideale sind Hauptideale")
end

function findnonprincipal2(p::BigInt,m::BigInt,n::BigInt)
#factor
	k,a=quadratic_field(-p)
	zk=maximal_order(k)
	i=nextprime(m)
	b=true
	while i<n
		PFZ=collect(keys(Hecke.factor(i*zk)))
		#PFZExp=collect(values(Hecke.factor(i*zk)))
		for j=1:length(PFZ)
			b=Hecke.isprincipal(Hecke.reduce_ideal(PFZ[j]))[1]
			print(b)
			if b==false
				return PFZ[j]
			end
		end
		i=i+1
	end
	print("alle Ideale sind Hauptideale")
end
				
				
