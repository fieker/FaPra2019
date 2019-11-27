function findnonprincipal(p::BigInt,m::BigInt,n::BigInt)
#prime_decomposition, using Primes
	k,a=quadratic_field(-p)
	zk=maximal_order(k)
	i=nextprime(m)
	while i<n
		x=prime_decomposition(zk,fmpz(i))
		#print(x)
		y=[]
		println(length(x))
		for j=1:length(x)
			if x[j][2]>1
				if x[j][1].is_principal!==1
					println("1v")
					#if isprincipal(x[j][1])[1]==false
					#	println("1n")					
						return x[j][1]
					#end
				end
			end
			push!(y,x[j][1].gen_one)
		end
		for j=1:length(y)
			for k=j:length(y)
				if y[j]==y[k]
					println("2v")
					if isprincipal(x[j][1])[1]==false
						println("2n")
						return x[j][1]
					end
				end
			end
		end
		#for j=1:length(x)
		#	b=isprincipal(x[j][1])
		#	if b[1]==false
		#		return x[j][1]
		#	end
		#end
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
				
				
